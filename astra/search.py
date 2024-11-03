
import os
import sys
import argparse
import pandas as pd, numpy as np
import pyhmmer
import subprocess
import collections
import shutil
from astra import initialize
from tqdm import tqdm
from platformdirs import user_config_dir
from copy import deepcopy
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
import asyncio
import logging
import time

# Pyhmmer-specific issues:
"""
- Pickle protocol not supported for sequences
- Thresholds are altered by pickling when parsing HMMs in parallel with ProcessPoolExecutor
"""
def get_results_attributes(result):
    bitscore = result.bitscore
    evalue = result.evalue
    cog = result.hmm_name
    c_evalue = result.c_evalue
    i_evalue = result.i_evalue
    query = result.sequence_id
    env_from = result.env_from
    env_to = result.env_to
    dom_bitscore = result.dom_bitscore
    return [query, cog, bitscore, evalue, c_evalue, i_evalue, env_from, env_to, bitscore]

#Store as a global so we don't have to define it multiple times
Result = collections.namedtuple("Result", ["sequence_id", "hmm_name", "bitscore", "evalue","c_evalue", "i_evalue", 
                                          "env_from", "env_to", "dom_bitscore"])

def extract_sequences(results_dataframes_dict, outdir):
    # Create tmp_ids directory within outdir
    tmp_ids_dir = os.path.join(outdir, 'tmp_ids')
    os.makedirs(tmp_ids_dir, exist_ok=True)
    
    # Create fastas directory within outdir
    fastas_dir = os.path.join(outdir, 'fastas')
    os.makedirs(fastas_dir, exist_ok=True)
    
    for genome_file, df in results_dataframes_dict.items():
        for hmm_name in df['hmm_name'].unique():
            # Extract the IDs corresponding to the current HMM
            ids_to_extract = df[df['hmm_name'] == hmm_name]['sequence_id'].tolist()
            
            # Write IDs to a temporary file
            idfile = os.path.join(tmp_ids_dir, f"{hmm_name}_ids.txt")
            with open(idfile, 'w') as f:
                f.write("\n".join(ids_to_extract))
            
            # Define the output FASTA file for hits
            hits_fasta = os.path.join(fastas_dir, f"{hmm_name}.faa")
            
            # Run pullseq command to extract sequences
            pullseq_cmd = f"cat {idfile} | pullseq -i {genome_file} -N >> {hits_fasta}"
            subprocess.run(pullseq_cmd, shell=True)

    # Remove tmp_ids directory
    shutil.rmtree(tmp_ids_dir)

def has_thresholds(x):
    flag = x.cutoffs.gathering_available() or \
           x.cutoffs.noise_available() or \
           x.cutoffs.trusted_available()

    return flag

def write_temp_result(result, temp_dir, sequence_name):
    tmp_file = os.path.join(temp_dir, f"{sequence_name}_results.tsv")
    with open(tmp_file, 'a') as f:
        f.write(f"{result.sequence_id}\t{result.hmm_name}\t{result.bitscore}\t{result.evalue}\t"
                f"{result.c_evalue}\t{result.i_evalue}\t{result.env_from}\t{result.env_to}\t"
                f"{result.dom_bitscore}\n")

def hmmsearch(protein_dict, hmms, threads, options, db_name=None):
    hmmsearch_kwargs = define_kwargs(options)
    tmp_dir = os.path.join(options['outdir'], 'tmp_results')
    os.makedirs(tmp_dir, exist_ok=True)

    def get_best_cutoff(hmm):
        if options['cascade']:
            cutoff_order = [
                hmmsearch_kwargs.get('preferred_cutoff', 'trusted'),
                'trusted', 'gathering', 'noise'
            ]
            for cutoff in cutoff_order:
                if getattr(hmm.cutoffs, f"{cutoff}_available")():
                    return cutoff
        elif 'bit_cutoffs' in hmmsearch_kwargs:
            if getattr(hmm.cutoffs, f"{hmmsearch_kwargs['bit_cutoffs']}_available")():
                return hmmsearch_kwargs['bit_cutoffs']
        return None

    for fasta_file, sequences in tqdm(protein_dict.items()):
        safe_filename = ''.join(c if c.isalnum() else '_' for c in os.path.basename(fasta_file))
        tmp_file = os.path.join(tmp_dir, f"{safe_filename}_results.tsv")

        # Write header to the temporary file
        with open(tmp_file, 'w') as f:
            f.write("sequence_id\thmm_name\tbitscore\tevalue\tc_evalue\ti_evalue\tenv_from\tenv_to\tdom_bitscore\n")

        # Group HMMs by their best available cutoff
        hmm_groups = {}
        for hmm in hmms:
            best_cutoff = get_best_cutoff(hmm)
            hmm_groups.setdefault(best_cutoff, []).append(hmm)

        # Perform searches for each group
        for cutoff, hmm_group in hmm_groups.items():
            kwargs = hmmsearch_kwargs.copy()
            if cutoff:
                kwargs['bit_cutoffs'] = cutoff
            elif 'bit_cutoffs' in kwargs:
                del kwargs['bit_cutoffs']

            if 'preferred_cutoff' in kwargs:
                del kwargs['preferred_cutoff']  # Remove this as it's not a valid pyhmmer parameter

            for hits in pyhmmer.hmmsearch(hmm_group, sequences, cpus=threads, **kwargs):
                process_hits(hits, tmp_file)

    return tmp_dir


def process_hits(hits, tmp_file):
    cog = hits.query_name.decode()
    with open(tmp_file, 'a') as f:
        for hit in hits:
            if hit.included:
                hit_name = hit.name.decode()
                full_bitscore = hit.score 
                full_evalue = hit.evalue
                for domain in hit.domains.reported:
                    f.write(f"{hit_name}\t{cog}\t{full_bitscore:.2f}\t{full_evalue:.2e}\t{domain.c_evalue:.2e}\t"
                            f"{domain.i_evalue:.2e}\t{domain.env_from}\t{domain.env_to}\t{domain.score:.2f}\n")


def extract_sequences_from_tmp(tmp_dir, protein_dict, outdir):
    fastas_dir = os.path.join(outdir, 'fastas')
    os.makedirs(fastas_dir, exist_ok=True)

    for filename in os.listdir(tmp_dir):
        if filename.endswith('_results.tsv'):
            file_path = os.path.join(tmp_dir, filename)
            df = pd.read_csv(file_path, sep='\t')
            
            for hmm_name in df['hmm_name'].unique():
                ids_to_extract = df[df['hmm_name'] == hmm_name]['sequence_id'].tolist()
                hits_fasta = os.path.join(fastas_dir, f"{hmm_name}.faa")

                with open(hits_fasta, 'w') as out_f:
                    for fasta_file, sequences in protein_dict.items():
                        for seq in sequences:
                            if seq.name.decode() in ids_to_extract:
                                out_f.write(f">{seq.name.decode()}\n{seq.sequence}\n")

def cleanup_temp_files(temp_dir):
    shutil.rmtree(temp_dir)
    print(f"Temporary files removed from {temp_dir}")



def parse_single_hmm(hmm_path):
    #Single-file parser for parallelization
    with pyhmmer.plan7.HMMFile(hmm_path) as hmm_file:
        return hmm_file.read()

# Modified to include explicit loop reference
async def parse_single_hmm_async(hmm_path, sem):
    #print(f"Processing {hmm_path}")  # Debug: Check if function is called
    async with sem:
        #print(f"Acquired semaphore for {hmm_path}")  # Debug: Check if semaphore is acquired
        loop = asyncio.get_event_loop()
        #print(f"Got event loop for {hmm_path}")  # Debug: Check if event loop is obtained
        result = await loop.run_in_executor(None, parse_single_hmm, hmm_path)
        #print(f"Executor completed for {hmm_path}")  # Debug: Check if executor has completed
        return result

def parse_hmms(hmm_in):
    #Checks first whether HMMs are provided as a single file or as a directory.

    hmms = []  # Initialize an empty list to store parsed HMMs
    print("Parsing HMMs...")
    # Check if hmm_in is a directory or a single file
    if os.path.isdir(hmm_in):
        if not os.listdir(hmm_in):
            print("hmm_in directory is empty.")
            logging.info('hmm_in directory is empty.')
            sys.exit(1)

        num_files = len(os.listdir(hmm_in))
        if num_files == 1:
            #Only one HMM file in input directory
            #Get full path to file
            hmm_path = os.path.join(hmm_in, os.listdir(hmm_in)[0])
            with pyhmmer.plan7.HMMFile(hmm_path) as hmm_file:
                #Works in case of single-model or multi-model HMM file
                hmms = list(hmm_file)

        else:
            hmm_files = list(filter(lambda x: x.endswith(('.hmm', '.HMM')), os.listdir(hmm_in)))
            hmm_paths = [os.path.join(hmm_in, hmm_file) for hmm_file in hmm_files]
            
            #I have tried!! Every possible method! To parallelize this!
            #It does not work. SINGLE THREADED IT IS!
            hmms = list(tqdm(map(parse_single_hmm, hmm_paths)))

            """
            #WITNESS THE RESULT OF MY FOLLY!
            #GAZE UPON MY MISDEEDS AND DESPAIR!
            loop = asyncio.get_event_loop()
            sem = asyncio.Semaphore(threads)  # Explicit loop reference
            
            async def gather_tasks():
                tasks = [parse_single_hmm_async(hmm_path, sem) for hmm_path in hmm_paths]
                return await asyncio.gather(*tasks)

            hmms = loop.run_until_complete(gather_tasks())

            """

    elif os.path.isfile(hmm_in):
        if os.path.getsize(hmm_in) == 0:
            print("hmm_in file is empty.")
            logging.info('hmm_in file is empty.')
            sys.exit(1)
        # Parse the single HMM file; handles multi-model files
        with pyhmmer.plan7.HMMFile(hmm_in) as hmm_file:
            hmms = list(hmm_file)
    else:
        print("Invalid HMM input.")
        logging.info("Invalid HMM input.")
        print("If you used pre-installed HMMs, check hmm_databases.json")
        logging.info("If you used pre-installed HMMs, check hmm_databases.json")
        print("Which is located in the databases directory.")
        logging.info("Which is located in the databases directory.")

        print("Thing that threw the error: {}".format(hmm_in))
        sys.exit(1)

    print("HMMs parsed.")

    return list(hmms)

def process_fasta(fasta_file):
    # Function to handle each file for parallelism
    with pyhmmer.easel.SequenceFile(fasta_file, digital=True, alphabet=pyhmmer.easel.Alphabet.amino()) as seq_file:
        sequences = seq_file.read_block()
    return fasta_file, sequences

def parse_protein_input(prot_in, threads):
    print("Parsing protein input sequences...")
    protein_dict = {}  # Initialize an empty dictionary to store parsed proteins
    
    # Check if prot_in is a directory or a single file
    if os.path.isdir(prot_in):
        if not os.listdir(prot_in):
            print("prot_in directory is empty.")
            logging.info("prot_in directory is empty.")
            sys.exit(1)

        # Initialize an empty dictionary to hold protein sequences
        protein_dict = {}

        #I formatted this as a loop because I was trying to parallelize it
        #but the sequence object for pyHMMER doesn't have pickle protocol support
        #Anyway I left it as a weird map with a tqdm you're welcome enjoy
        results = list(map(process_fasta, 
            tqdm(
                [os.path.join(prot_in, x) for x in os.listdir(prot_in)]
                )
            ))

        # Populate the protein_dict
        for fasta_path, sequences in results:
            protein_dict[fasta_path] = sequences
    elif os.path.isfile(prot_in):
        if os.path.getsize(prot_in) == 0:
            print("prot_in file is empty.")
            logging.info("prot_in file is empty.")
            sys.exit(1)
        # Parse the single protein FASTA file
        with pyhmmer.easel.SequenceFile(prot_in, digital=True) as seq_file:
            sequences = seq_file.read_block()
        protein_dict[prot_in] = sequences
    else:
        print("Invalid input for prot_in.")
        logging.info("Invalid input for prot_in.")
        sys.exit(1)
    
    return protein_dict

def get_installed_hmm_paths(hmm_names):
    #Loads relevant information from astra DB json file
    parsed_json = initialize.load_config()
    hmm_paths = []
    for db in parsed_json['db_urls']:
        if db['name'] in hmm_names and db['installed']:
            hmm_paths.append(os.path.join(db['installation_dir'], db['name']))
    return hmm_paths

def define_kwargs(options):
    kwargs = {}
    
    if options['cascade']:
        kwargs['bit_cutoffs'] = 'cascade'
        if options['cut_tc']:
            kwargs['preferred_cutoff'] = 'trusted'
        elif options['cut_ga']:
            kwargs['preferred_cutoff'] = 'gathering'
        elif options['cut_nc']:
            kwargs['preferred_cutoff'] = 'noise'
        else:
            kwargs['preferred_cutoff'] = 'trusted'  # Default to trusted if no specific cutoff is specified
    elif options['cut_ga']:
        kwargs['bit_cutoffs'] = 'gathering'
    elif options['cut_nc']:
        kwargs['bit_cutoffs'] = 'noise'
    elif options['cut_tc']:
        kwargs['bit_cutoffs'] = 'trusted'


    #Numerical threshold parameters
    if options['bitscore'] is not None:
        #Make sure it's the right format, or castable as such!
        if not isinstance(options['bitscore'], float):
            try:
                kwargs['T'] = float(options['bitscore'])
            except ValueError:
                print("Error: bitscore threshold must be a float or castable as a float.")
                logging.info("Error: bitscore threshold must be a float or castable as a float.")

    if options['domE'] is not None:
        #Make sure it's the right format, or castable as such!
        if not isinstance(options['domE'], float):
            try:
                kwargs['domE'] = float(options['domE'])
            except ValueError:
                print("Error: domE must be a float or castable to float.")
                logging.info("Error: domE must be a float or castable to float.")

    if options['domT'] is not None:
        #Make sure it's the right format, or castable as such!
        if not isinstance(options['domT'], float):
            try:
                kwargs['domT'] = float(options['domT'])
            except ValueError:
                print("Error: domT must be a float or castable to float.")
                logging.info("Error: domT must be a float or castable to float.")

    if options['incE'] is not None:
        #Make sure it's the right format, or castable as such!
        if not isinstance(options['incE'], float):
            try:
                kwargs['incE'] = float(options['incE'])
            except ValueError:
                print("Error: domT must be a float or castable to float.")
                logging.error("Error: domT must be a float or castable to float.")

    if options['incT'] is not None:
        #Make sure it's the right format, or castable as such!
        if not isinstance(options['incT'], float):
            try:
                kwargs['incT'] = float(options['incT'])
            except ValueError:
                print("Error: incT must be a float or castable to float.")
                logging.error("Error: incT must be a float or castable to float.")

    if options['incdomE'] is not None:
        #Make sure it's the right format, or castable as such!
        if not isinstance(options['incdomE'], float):
            try:
                kwargs['incdomE'] = float(options['incdomE'])
            except ValueError:
                print("Error: incdomE must be a float or castable to float.")
                logging.error("Error: incdomE must be a float or castable to float.")

    if options['incdomT'] is not None:
        #Make sure it's the right format, or castable as such!
        if not isinstance(options['incdomT'], float):
            try:
                kwargs['incdomT'] = float(options['incdomT'])
            except ValueError:
                print("Error: incdomT must be a float or castable to float.")
                logging.error("Error: incdomT must be a float or castable to float.")

    if options['evalue'] is not None:
        #Make sure it's the right format, or castable as such!
        if not isinstance(options['incdomT'], float):
            try:
                kwargs['E'] = float(options['evalue'])
            except ValueError:
                print("Error: evalue must be a float or castable to float.")
                logging.error("Error: evalue must be a float or castable to float.")

    return kwargs

def sanitize_filename(filename):
    return ''.join(c for c in filename if c.isalnum() or c in '._- ')

def create_temp_directory(outdir):
    temp_dir = os.path.join(outdir, 'tmp')
    os.makedirs(temp_dir, exist_ok=True)
    return temp_dir

def combine_results(tmp_dir, output_file):
    all_results = []
    for filename in os.listdir(tmp_dir):
        if filename.endswith('_results.tsv'):
            file_path = os.path.join(tmp_dir, filename)
            df = pd.read_csv(file_path, sep='\t')
            print(f"Sample from {filename}:")
            print(df.head().to_string())  # This will print the first few rows
            all_results.append(df)

    if all_results:
        combined_df = pd.concat(all_results, ignore_index=True).sort_values(by=['sequence_id', 'env_from'], ascending=[True, True])
        combined_df.to_csv(output_file, sep='\t', index=False)
        print(f"Combined results written to {output_file}")
        print("Sample from combined results:")
        print(combined_df.head().to_string())  # This will print the first few rows of the combined results
    else:
        print("No results found to combine.")

def main(args):
    t1 = time.time()
    hmm_in = args.hmm_in
    prot_in = args.prot_in 
    outdir = args.outdir
    log_file_path = os.path.join(outdir, 'astra_search_log.txt')

    if not os.path.exists(outdir):
        os.makedirs(outdir)
        if args.write_seqs:
            os.makedirs(os.path.join(outdir, 'fastas'))

    logging.basicConfig(filename=log_file_path, level=logging.INFO,
                        format='%(asctime)s %(levelname)s: %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')

    hmmsearch_options = {
        "cascade": args.cascade,
        "cut_ga": args.cut_ga,
        "cut_nc": args.cut_nc,
        "cut_tc": args.cut_tc,
        "evalue": args.evalue,
        "bitscore": args.bitscore,
        "domE": args.domE,
        "domT": args.domT,
        "incE": args.incE,
        "incT": args.incT,
        "incdomE": args.incdomE,
        "incdomT": args.incdomT,
        "outdir": outdir,
        "meta": args.meta
    }

    if hmm_in is None and args.installed_hmms is None:
        error_out = "Either a user-provided or pre-installed HMM database must be specified."
        print(error_out)
        logging.error(error_out)
        sys.exit(1)

    protein_dict = parse_protein_input(prot_in, args.threads)

    if hmm_in is not None:
        print("Searching with user-provided HMM(s)...")
        logging.info("Searching with user-provided HMM(s)...")
        user_hmms = parse_hmms(args.hmm_in)
        tmp_dir = hmmsearch(protein_dict, user_hmms, args.threads, hmmsearch_options)
        combine_results(tmp_dir, os.path.join(outdir, 'user_hmms_hits_df.tsv'))
        if args.write_seqs:
            extract_sequences_from_tmp(tmp_dir, protein_dict, outdir)
        del user_hmms

    if args.installed_hmms is not None:
        installed_hmm_names = args.installed_hmms.split(',') if ',' in args.installed_hmms else [args.installed_hmms]
        print(f"Searching with pre-installed HMMs: {', '.join(installed_hmm_names)}")
        logging.info(f"Searching with pre-installed HMMs: {', '.join(installed_hmm_names)}")

        parsed_json = initialize.load_config()

        if 'all_prot' in installed_hmm_names:
            installed_hmm_names = [db['name'] for db in parsed_json['db_urls'] if db['molecule_type'] == 'protein' and db['installed']]

        for hmm_db in installed_hmm_names:
            installed_hmm_in = next((item for item in parsed_json['db_urls'] if item["name"] == hmm_db), None)
            if installed_hmm_in is not None:
                installation_dir = installed_hmm_in['installation_dir']
                db_hmms = parse_hmms(installation_dir)
                tmp_dir = hmmsearch(protein_dict, db_hmms, args.threads, hmmsearch_options, hmm_db)
                combine_results(tmp_dir, os.path.join(outdir, f'{hmm_db}_hits_df.tsv'))
                if args.write_seqs:
                    extract_sequences_from_tmp(tmp_dir, protein_dict, outdir)
            else:
                print(f"No installation_dir specified for db {hmm_db}")
                logging.info(f"No installation_dir specified for db {hmm_db}")

    # Clean up temporary directory
    cleanup_temp_files(os.path.join(outdir, 'tmp_results'))

    time_printout = f"Process took {time.time()-t1} seconds."
    print(time_printout)
    logging.info(time_printout)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="ASTRA search tool")
    # Add existing arguments here
    args = parser.parse_args()
    main(args)


#TODO:
"""
- Multithread extract_sequences
"""
