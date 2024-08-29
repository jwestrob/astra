
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

def hmmsearch(protein_dict, hmms, threads, options, individual_results_dir=None, db_name=None):
    hmmsearch_kwargs = define_kwargs(options)

    def get_best_cutoff(hmm):
        if options['cascade']:
            cutoff_order = [
                hmmsearch_kwargs.get('preferred_cutoff', 'trusted'),
                'trusted',
                'gathering',
                'noise'
            ]
            for cutoff in cutoff_order:
                if getattr(hmm.cutoffs, f"{cutoff}_available")():
                    return cutoff
        elif 'bit_cutoffs' in hmmsearch_kwargs:
            if getattr(hmm.cutoffs, f"{hmmsearch_kwargs['bit_cutoffs']}_available")():
                return hmmsearch_kwargs['bit_cutoffs']
        return None

    results_dataframes = {}

    for fasta_file, sequences in tqdm(protein_dict.items()):
        results = []

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
                process_hits(hits, results)

        result_df = pd.DataFrame(
            list(map(get_results_attributes, results)), 
            columns=["sequence_id", "hmm_name", "bitscore", "evalue", "c_evalue", 
                     "i_evalue", "env_from", "env_to", "dom_bitscore"]
        )

        if individual_results_dir:
            basename_fasta = os.path.basename(fasta_file)
            output_filename = f"{basename_fasta}_{db_name}_results.tsv" if db_name else f"{basename_fasta}_results.tsv"
            result_df.to_csv(os.path.join(individual_results_dir, output_filename), sep='\t', index=False)
        elif options.get('meta', False):
            basename_fasta = os.path.basename(fasta_file)
            output_filename = f"{basename_fasta}_{db_name}_results.tsv" if db_name else f"{basename_fasta}_results.tsv"
            result_df.to_csv(os.path.join(options['outdir'], output_filename), sep='\t', index=False)
        else:
            results_dataframes[fasta_file] = result_df

        # Clear results to free up memory
        results.clear()

    return results_dataframes if not individual_results_dir and not options.get('meta', False) else None

def process_hits(hits, results):
    cog = hits.query_name.decode()
    for hit in hits:
        if hit.included:
            hit_name = hit.name.decode()
            full_bitscore = hit.score 
            full_evalue = hit.evalue
            for domain in hit.domains.reported:
                results.append(Result(hit_name, cog, full_bitscore, full_evalue, domain.c_evalue, 
                              domain.i_evalue, domain.env_from, domain.env_to, domain.score))

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


def combine_results(individual_results_dir, output_file):
    all_results = []
    for filename in os.listdir(individual_results_dir):
        if filename.endswith('_results.tsv'):
            file_path = os.path.join(individual_results_dir, filename)
            df = pd.read_csv(file_path, sep='\t')
            all_results.append(df)
    
    combined_df = pd.concat(all_results, ignore_index=True)
    combined_df.to_csv(output_file, sep='\t', index=False)
    print(f"Combined results written to {output_file}")

def main(args):
    t1 = time.time()
    # Required arguments
    hmm_in = args.hmm_in
    prot_in = args.prot_in 

    #boolean; indicates input is metagenomic files
    global meta
    meta = args.meta

    #Set this as global; we don't want to have to pass it
    global outdir 
    outdir = args.outdir
    log_file_path = os.path.join(outdir, 'astra_search_log.txt')

    # Check if the output directory already exists
    if not os.path.exists(outdir):
        os.makedirs(outdir)
        if args.write_seqs:
            os.makedirs(os.path.join(outdir, 'fastas'))  # Also create a 'fastas' folder within the output directory
    
    # Create individual_results directory if --individual_results is specified
    if args.individual_results:
        individual_results_dir = os.path.join(outdir, 'individual_results')
        os.makedirs(individual_results_dir, exist_ok=True)
    else:
        individual_results_dir = None

    logging.basicConfig(filename=log_file_path, level=logging.INFO,
                        format='%(asctime)s %(levelname)s: %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')
    installed_hmms = args.installed_hmms

    # Optional arguments
    evalue = args.evalue
    bitscore = args.bitscore

    # Boolean flags
    cut_ga = args.cut_ga
    cut_nc = args.cut_nc
    cut_tc = args.cut_tc

    write_seqs = args.write_seqs

    #again i am too lazy to pass this parameter in a function call SUE ME
    global threads
    threads = args.threads

    #initialize default options
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
    }


    if hmm_in is None and installed_hmms is None:
        error_out = "Either a user-provided or pre-installed HMM database must be specified. You know better."
        print(error_out)
        logging.error(error_out)
        sys.exit(1)

    # Check if more than one of --evalue, --bitscore, --cut_nc, --cut_tc, and --cut_ga are specified
    specified_flags = [args.cut_nc, args.cut_tc, args.cut_ga]
    if sum(specified_flags) > 1:
        print("Error: You can only specify one of --bitscore, --cut_nc, --cut_tc, and --cut_ga.")
        logging.info("Error: You can only specify one of --bitscore, --cut_nc, --cut_tc, and --cut_ga.")
        print("If you specify a bitscore threshold and a pre-defined cutoff (e.g. --cut_ga) the pre-defined cutoff will be used")
        logging.info("If you specify a bitscore threshold and a pre-defined cutoff (e.g. --cut_ga) the pre-defined cutoff will be used")
        print("where available, otherwise the specified bitscore threshold will be used.")
        logging.info("where available, otherwise the specified bitscore threshold will be used.")





    #Check protein input and parse
    protein_dict = parse_protein_input(prot_in, threads)
    
    if hmm_in is not None:
        print("Searching with user-provided HMM(s)...")
        logging.info("Searching with user-provided HMM(s)...")
        #Check HMM input and parse
        user_hmms = parse_hmms(args.hmm_in)
        #Obtain dictionary containing results dataframes for each input FASTA
        results_dataframes_dict = hmmsearch(protein_dict, user_hmms, threads, hmmsearch_options, individual_results_dir)
        
        if args.write_seqs:
            extract_sequences(results_dataframes_dict, outdir)

        if individual_results_dir:
            combine_results(individual_results_dir, os.path.join(outdir, 'user_hmms_hits_df.tsv'))
        else:
            all_results_df = pd.concat([results_dataframes_dict[key] for key in results_dataframes_dict.keys()])
            all_results_df.to_csv(os.path.join(outdir,'user_hmms_hits_df.tsv'), sep='\t', index=False)
        del user_hmms
        
    if installed_hmms is not None:
        #check HMM input and parse

        # Step 1: Get paths for installed HMM databases
        installed_hmm_paths = []
        if ',' in installed_hmms:
            installed_hmm_names = installed_hmms.split(',')
        else:
            installed_hmm_names = [installed_hmms]  # Single element list

        #Two conditions in case there isn't a , in the installed_hmm_names
        if ',' in installed_hmm_names:
            print("Searching with pre-installed HMMs: ", ', '.join(installed_hmm_names))
            logging.info("Searching with pre-installed HMMs: ", ', '.join(installed_hmm_names))
        else:
            print("Searching with pre-installed HMMs: {}".format(installed_hmm_names[0]))
            logging.info("Searching with pre-installed HMMs: {}".format(installed_hmm_names[0]))

        #Load JSON with database and procedural information
        parsed_json = initialize.load_config()

        if 'all_prot' in installed_hmm_names:

            # Replace 'all_prot' with all installed protein HMM database names
            installed_hmm_names = []
            installed_hmm_paths = []
            for db in parsed_json['db_urls']:
                if db['molecule_type'] == 'protein' and db['installed']:
                    installed_hmm_names.append(db['name'])

        for hmm_db in installed_hmm_names:

            installed_hmm_in = next((item for item in parsed_json['db_urls'] if item["name"] == hmm_db), None)
            if installed_hmm_in is not None:
                installation_dir = installed_hmm_in['installation_dir']
                db_hmms = parse_hmms(installation_dir)
            else:
                #No installation_dir specified; print this and move on
                print("No installation_dir specified for db " + hmm_db)
                logging.info("No installation_dir specified for db " + hmm_db)
                continue
            #if we're in meta mode, we don't want to keep all that shit in memory
            #and the hmmsearch function will write a file for each DB and each protein file
            #because they're huge
            if not meta:
                results_dataframes_dict = hmmsearch(protein_dict, db_hmms, threads, hmmsearch_options, hmm_db)
            else:
                hmmsearch(protein_dict, db_hmms, threads, hmmsearch_options, hmm_db)

            if args.write_seqs:
                extract_sequences(results_dataframes_dict, outdir)

            if not meta:
                db_results_df = pd.concat([results_dataframes_dict[key] for key in results_dataframes_dict.keys()])
                db_results_df.to_csv(os.path.join(outdir,hmm_db + '_hits_df.tsv'), sep='\t', index=False)
    time_printout  = "Process took {} seconds.".format(time.time()-t1)
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