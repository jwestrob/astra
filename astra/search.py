
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
from concurrent.futures import ProcessPoolExecutor, as_completed


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
    return x.cutoffs.gathering_available() is not None or \
           x.cutoffs.noise_available() is not None or \
           x.cutoffs.trusted_available() is not None

def hmmsearch(protein_dict, hmms, threads, options):
    #Runs HMMscan on all provided FASTA files using 'threads' threads
    #uses default parameters unless specified


    results_dataframes = {}  # Initialize an empty dictionary to store results as DataFrames

    #Construct search options; need to make sure names are consistent with pyHMMER
    #and we only specify one bitscore threshold
    hmmsearch_kwargs = define_kwargs(options)
    print(hmmsearch_kwargs)
    if hmmsearch_kwargs['bit_cutoffs'] != None:
        print(hmmsearch_kwargs)
        #pyHMMER rightly throws an error when you try to use thresholds that don't exist in the model.
        #Let's separate these out because often a single set of HMMs will contain thresholded
        #as well as unthresholded models.
        print("Separating thresholded and non-thresholded HMMs...")

        with ProcessPoolExecutor(threads) as executor:
            # Create a Boolean mask indicating which HMMs have thresholds
            has_thresholds_mask = list(executor.map(has_thresholds, hmms))

            # Convert the mask to a NumPy array for efficient indexing
            has_thresholds_mask_np = np.array(has_thresholds_mask)

            # Filter HMMs with thresholds using the mask
            hmms_with_thresholds = np.array(hmms)[has_thresholds_mask_np].tolist()

            # Filter HMMs without thresholds using the inverse of the mask
            hmms_without_thresholds = np.array(hmms)[~has_thresholds_mask_np].tolist()



        hmmsearch_kwargs_nothreshold = hmmsearch_kwargs
        #Remove this; having a None value for bit_cutoffs will annoy pyHMMER
        del hmmsearch_kwargs_nothreshold['bit_cutoffs']

        #Remove evalue and bitscore thresholds that might be otherwise imposed when running with predefined cutoff scores
        #And remove other thresholds because pyhmmer gets mad if you specify nonetype
        bit_cutoff = hmmsearch_kwargs['bit_cutoffs']
        hmmsearch_kwargs = {'bit_cutoffs':bit_cutoff}
    else:
        hmms_without_thresholds = hmms
        hmms_with_thresholds = None

    print("Searching...")
    for fasta_file, sequences in tqdm(protein_dict.items()):
        results = []
        
        if hmms_with_thresholds is not None:
            print("scanning with thresholded HMMs...")
            print(hmmsearch_kwargs)
            # Run the thresholded HMMs
            for hits in pyhmmer.hmmsearch(hmms_with_thresholds, sequences, cpus=threads, **hmmsearch_kwargs):
                cog = hits.query_name.decode()
                for hit in hits:
                    if hit.included:
                        hit_name = hit.name.decode()
                        full_bitscore = hit.score 
                        full_evalue = hit.evalue
                        for domain in hit.domains.reported:
                            results.append(Result(hit_name, cog, full_bitscore, full_evalue, domain.c_evalue, 
                                  domain.i_evalue, domain.env_from, domain.env_to, domain.score))

        if hmms_without_thresholds is not None:
            print("scanning with unthresholded HMMs...")
            #Run the unthresholded HMMs, making sure to specify bit_cutoffs=None
            for hits in pyhmmer.hmmsearch(hmms_without_thresholds, sequences, cpus=threads, **hmmsearch_kwargs_nothreshold):
                cog = hits.query_name.decode()
                for hit in hits:
                    if hit.included:
                        hit_name = hit.name.decode()
                        full_bitscore = hit.score 
                        full_evalue = hit.evalue
                        for domain in hit.domains.reported:
                            results.append(Result(hit_name, cog, full_bitscore, full_evalue, domain.c_evalue, 
                                  domain.i_evalue, domain.env_from, domain.env_to, domain.score))
                    
        # Convert the results to a DataFrame
        #Is it necessary to cast it as a list?
        result_df = pd.DataFrame(list(map(get_results_attributes, results)), columns=["sequence_id", "hmm_name", "bitscore", "evalue","c_evalue", "i_evalue", "env_from", "env_to", "dom_bitscore"])
        
        if meta == False:
            # Store the DataFrame in the dictionary
            results_dataframes[fasta_file] = result_df
        else:
            try:
                result_df.to_csv(os.path.join(outdir, fasta_file + '_' + hmm_db + 'results.tsv'), sep='\t', index=False)
            except Error as e:
                print(Error)
                result_df.to_csv(fasta_file + '_' + hmm_db + '.results.tsv', sep='\t', index=False)
    
    return results_dataframes

def parse_single_hmm(hmm_path):
    #Single-file parser for parallelization
    with pyhmmer.plan7.HMMFile(hmm_path) as hmm_file:
        return list(hmm_file)

def parse_hmms(hmm_in):
    #Checks first whether HMMs are provided as a single file or as a directory.

    hmms = []  # Initialize an empty list to store parsed HMMs
    
    # Check if hmm_in is a directory or a single file
    if os.path.isdir(hmm_in):
        if not os.listdir(hmm_in):
            print("hmm_in directory is empty.")
            sys.exit(1)
        
        # Generate list of HMM files
        hmm_files = list(filter(lambda x: x.endswith(('.hmm', '.HMM')), os.listdir(hmm_in)))
        hmm_paths = [os.path.join(hmm_in, hmm_file) for hmm_file in hmm_files]

        # Initialize progress bar
        pbar = tqdm(total=len(hmm_files), desc="Parsing HMMs")

        # Parse each HMM file in the directory
        with ProcessPoolExecutor(threads) as executor:
            future_to_hmm = {executor.submit(parse_single_hmm, hmm_path): hmm_path for hmm_path in hmm_paths}
            for future in as_completed(future_to_hmm):
                hmms.extend(future.result())
                pbar.update(1)

        pbar.close()
    elif os.path.isfile(hmm_in):
        if os.path.getsize(hmm_in) == 0:
            print("hmm_in file is empty.")
            sys.exit(1)
        # Parse the single HMM file
        with pyhmmer.plan7.HMMFile(hmm_in) as hmm_file:
            hmms = list(hmm_file)
    else:
        print("Invalid input for hmm_in.")
        sys.exit(1)
    
    return hmms

def process_fasta(fasta_file):
    # Function to handle each file for parallelism
    with pyhmmer.easel.SequenceFile(fasta_file, digital=True) as seq_file:
        sequences = seq_file.read_block()
    return fasta_file, sequences

def parse_protein_input(prot_in, threads):
    print("Parsing protein input sequences...")
    protein_dict = {}  # Initialize an empty dictionary to store parsed proteins
    
    # Check if prot_in is a directory or a single file
    if os.path.isdir(prot_in):
        if not os.listdir(prot_in):
            print("prot_in directory is empty.")
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
            sys.exit(1)
        # Parse the single protein FASTA file
        with pyhmmer.easel.SequenceFile(prot_in, digital=True) as seq_file:
            sequences = seq_file.read_block()
        protein_dict[prot_in] = sequences
    else:
        print("Invalid input for prot_in.")
        sys.exit(1)
    
    return protein_dict

def get_installed_hmm_paths(hmm_names):
    #Loads relevant information from astra DB json file
    parsed_json = initialize.load_json()
    hmm_paths = []
    for db in parsed_json['db_urls']:
        if db['name'] in hmm_names and db['installed']:
            hmm_paths.append(os.path.join(db['installation_dir'], db['name']))
    return hmm_paths

def define_kwargs(options):
    kwargs = {}

    #Calibrated threshold parameters
    if options['cut_ga']:
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

    if options['domE'] is not None:
        #Make sure it's the right format, or castable as such!
        if not isinstance(options['domE'], float):
            try:
                kwargs['domE'] = float(options['domE'])
            except ValueError:
                print("Error: domE must be a float or castable to float.")

    if options['domT'] is not None:
        #Make sure it's the right format, or castable as such!
        if not isinstance(options['domT'], float):
            try:
                kwargs['domT'] = float(options['domT'])
            except ValueError:
                print("Error: domT must be a float or castable to float.")

    if options['incE'] is not None:
        #Make sure it's the right format, or castable as such!
        if not isinstance(options['incE'], float):
            try:
                kwargs['incE'] = float(options['incE'])
            except ValueError:
                print("Error: domT must be a float or castable to float.")

    if options['incT'] is not None:
        #Make sure it's the right format, or castable as such!
        if not isinstance(options['incT'], float):
            try:
                kwargs['incT'] = float(options['incT'])
            except ValueError:
                print("Error: domT must be a float or castable to float.")

    if options['incdomE'] is not None:
        #Make sure it's the right format, or castable as such!
        if not isinstance(options['incdomE'], float):
            try:
                kwargs['incdomE'] = float(options['incdomE'])
            except ValueError:
                print("Error: domT must be a float or castable to float.")

    if options['incdomT'] is not None:
        #Make sure it's the right format, or castable as such!
        if not isinstance(options['incdomT'], float):
            try:
                kwargs['incdomT'] = float(options['incdomT'])
            except ValueError:
                print("Error: domT must be a float or castable to float.")

    if options['evalue'] is not None:
        #Make sure it's the right format, or castable as such!
        if not isinstance(options['incdomT'], float):
            try:
                kwargs['E'] = float(options['evalue'])
            except ValueError:
                print("Error: domT must be a float or castable to float.")

    return kwargs


def main(args):
    # Required arguments
    hmm_in = args.hmm_in
    prot_in = args.prot_in 

    #boolean; indicates input is metagenomic files
    global meta
    meta = args.meta

    #Set this as global; we don't want to have to pass it
    global outdir 
    outdir = args.outdir

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
    "cut_ga":cut_ga,
    "cut_nc":cut_nc,
    "cut_tc":cut_tc,
    "evalue":evalue,
    "bitscore":bitscore,
    "domE":args.domE,
    "domT":args.domT,
    "incE":args.incE,
    "incT":args.incT,
    "incdomE":args.incdomE,
    "incdomT":args.incdomT,
    }

    if hmm_in is None and installed_hmms is None:
        print("Either a user-provided or pre-installed HMM database must be specified. You know better.")
        sys.exit(1)

    # Check if more than one of --evalue, --bitscore, --cut_nc, --cut_tc, and --cut_ga are specified
    specified_flags = [args.cut_nc, args.cut_tc, args.cut_ga]
    if sum(specified_flags) > 1:
        print("Error: You can only specify one of --bitscore, --cut_nc, --cut_tc, and --cut_ga.")
        print("If you specify a bitscore threshold and a pre-defined cutoff (e.g. --cut_ga) the pre-defined cutoff will be used")
        print("where available, otherwise the specified bitscore threshold will be used.")

    # Check if the output directory already exists
    if not os.path.exists(outdir):
        os.makedirs(outdir)
        if write_seqs:
            os.makedirs(os.path.join(outdir, 'fastas'))  # Also create a 'fastas' folder within the output directory




    #Check protein input and parse
    protein_dict = parse_protein_input(prot_in, threads)

    if hmm_in is not None:
        print("Searching with user-provided HMM(s)...")
        #Check HMM input and parse
        user_hmms = parse_hmms(args.hmm_in)
        
        #Obtain dictionary containing results dataframes for each input FASTA
        results_dataframes_dict = hmmsearch(protein_dict, user_hmms, threads, hmmsearch_options)
        
        if args.write_seqs:
            extract_sequences(results_dataframes_dict, outdir)

        all_results_df = pd.concat([results_dataframes_dict[key] for key in results_dataframes_dict.keys()])
        all_results_df.to_csv(os.path.join(outdir,'all_hits_df.tsv'), sep='\t', index=False)

    if installed_hmms is not None:
        #check HMM input and parse

        # Step 1: Get paths for installed HMM databases
        installed_hmm_paths = []
        if ',' in installed_hmms:
            installed_hmm_names = installed_hmms.split(',')
        else:
            installed_hmm_names = [installed_hmms]  # Single element list
        print("Searching with pre-installed HMMs: ", ', '.join(installed_hmm_names))
        #Load JSON with database and procedural information
        parsed_json = initialize.load_json()

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
                continue

            #if we're in meta mode, we don't want to keep all that shit in memory
            #and the hmmsearch function will write a file for each DB and each protein file
            #because they're huge
            if not meta:
                results_dataframes_dict = hmmsearch(protein_dict, db_hmms, threads, hmmsearch_options)
            else:
                hmmsearch(protein_dict, db_hmms, threads, hmmsearch_options)

            if args.write_seqs:
                extract_sequences(results_dataframes_dict, outdir)

            if not meta:
                db_results_df = pd.concat([results_dataframes_dict[key] for key in results_dataframes_dict.keys()])
                db_results_df.to_csv(os.path.join(outdir,hmm_db + '_hits_df.tsv'), sep='\t', index=False)

if __name__ == "__main__":
    main()


#TODO:
"""
- implement thresholding; pass arguments to 'hmmsearch' function
- Implement custom thresholding for databases
"""