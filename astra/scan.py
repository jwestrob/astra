
import argparse
import asyncio
import collections
import logging
import os
import shutil
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from copy import deepcopy

import numpy as np
import pandas as pd
import pyhmmer
from tqdm import tqdm

from astra import initialize
from astra.search import(
    Result,
    get_results_attributes,
    extract_sequences,
    has_thresholds,
    parse_hmms,
    parse_protein_input,
)

# Pyhmmer-specific issues:
"""
- Pickle protocol not supported for sequences
- Thresholds are altered by pickling when parsing HMMs in parallel with ProcessPoolExecutor
"""

def hmmscan(protein_dict, hmms, threads, options, db_name = None):
    #Runs HMMscan on all provided FASTA files using 'threads' threads
    #uses default parameters unless specified


    results_dataframes = {}  # Initialize an empty dictionary to store results as DataFrames

    #Construct scan options; need to make sure names are consistent with pyHMMER
    #and we only specify one bitscore threshold
    hmmscan_kwargs = define_kwargs(options)

    #Do we need to check whether we have a mixture of thresholded and non-thresholded models?
    if 'bit_cutoffs' in hmmscan_kwargs and not db_name in ['PFAM', 'FOAM']:
        #pyHMMER rightly throws an error when you try to use thresholds that don't exist in the model.
        #Let's separate these out because often a single set of HMMs will contain thresholded
        #as well as unthresholded models.
        print("Separating thresholded and non-thresholded HMMs...")
        logging.info("Separating thresholded and non-thresholded HMMs...")

        with ProcessPoolExecutor(threads) as executor:
            # Create a Boolean mask indicating which HMMs have thresholds
            has_thresholds_mask = list(executor.map(has_thresholds, hmms))

        # Convert the mask to a NumPy array for efficient indexing
        has_thresholds_mask_np = np.array(has_thresholds_mask)

        # Filter HMMs with thresholds using the mask
        hmms_with_thresholds = np.array(hmms)[has_thresholds_mask_np].tolist()

        # Filter HMMs without thresholds using the inverse of the mask
        hmms_without_thresholds = np.array(hmms)[~has_thresholds_mask_np].tolist()

        if len(hmms_with_thresholds) == 0:
            print("Bitscore cutoffs were specified, but specified HMMs do not contain these thresholds.")
            print("Defaulting to other specified threshold parameters (if none were specified, none will be applied)...")
            logging.info("Bitscore cutoffs were specified, but specified HMMs do not contain these thresholds.")
            logging.info("Defaulting to other specified threshold parameters (if none were specified, none will be applied)...")
            hmms_with_thresholds = None

        if len(hmms_without_thresholds) == 0:
            hmms_without_thresholds = None



        #This was a bit tricky. pyHMMER doesn't like NoneType for bit_cutoffs
        #And I don't want to specify any other thresholds for models with cutoffs
        #So we have to isolate that parameter, and remove it from the kwargs 
        #Used in a scan for non-thresholded models
        bit_cutoff = hmmscan_kwargs['bit_cutoffs']

        del hmmscan_kwargs['bit_cutoffs']
    elif db_name in ['PFAM', 'FOAM'] and 'bit_cutoffs' in hmmscan_kwargs:

        #These dbs have thresholds for every HMM
        hmms_with_thresholds = hmms
        hmms_without_thresholds = None
        bit_cutoff = hmmscan_kwargs['bit_cutoffs']
    else:
        #All unthresholded
        hmms_with_thresholds = None
        hmms_without_thresholds = hmms
        

    print("Scanning...")

    for fasta_file, sequences in tqdm(protein_dict.items()):
        results = []
        
        if hmms_with_thresholds is not None:
            #print("Scanning with {} thresholded HMMs...".format(len(hmms_with_thresholds)))
            # Run the thresholded HMMs
            for hits in pyhmmer.hmmscan(sequences, hmms_with_thresholds, cpus=threads, bit_cutoffs=bit_cutoff):
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
            #print("Scanning with {} unthresholded HMMs...".format(len(hmms_without_thresholds)))
            #print(hmmscan_kwargs)
            #Run the unthresholded HMMs, making sure to specify bit_cutoffs=None
            for hits in pyhmmer.hmmscan(hmms_without_thresholds, sequences, cpus=threads, **hmmscan_kwargs):
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
            #If meta is true, we don't want to hold all the results in RAM. We want to write an output file for every DB-metagenome scan.
            basename_fasta = os.path.basename(fasta_file)
            try:
                #Make sure the outdir exists and db_name is specified
                result_df.to_csv(os.path.join(outdir, basename_fasta + '_' + db_name + '_results.tsv'), sep='\t', index=False)
            except:
                #Hey man idk, maybe it doesn't? Maybe you called scan as a function from a python script?
                #If so, write output files to the current working directory instead.
                if db_name is None:
                    #Is there no db_name and meta is specified?
                    result_df.to_csv(os.path.join(outdir, basename_fasta + '_results.tsv'), sep='\t', index=False)
                else:
                    result_df.to_csv(fasta_file + '_' + db_name + '.results.tsv', sep='\t', index=False)
    if meta == False:
        return results_dataframes
    else:
        return


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
    log_file_path = os.path.join(outdir, 'astra_scan_log.txt')
    write_seqs = args.write_seqs

    # Check if the output directory already exists
    if not os.path.exists(outdir):
        os.makedirs(outdir)
        if write_seqs:
            os.makedirs(os.path.join(outdir, 'fastas'), exist_ok = True)  # Also create a 'fastas' folder within the output directory


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


    #again i am too lazy to pass this parameter in a function call SUE ME
    global threads
    threads = args.threads

    #initialize default options
    hmmscan_options = {
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
        print("Scanning with user-provided HMM(s)...")
        logging.info("Scanning with user-provided HMM(s)...")
        #Check HMM input and parse
        user_hmms = parse_hmms(args.hmm_in)
        #Obtain dictionary containing results dataframes for each input FASTA
        results_dataframes_dict = hmmscan(protein_dict, user_hmms, threads, hmmscan_options)
        
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
        if ',' in installed_hmm_names:
            print("Scanning with pre-installed HMMs: ", ', '.join(installed_hmm_names))
            logging.info("Scanning with pre-installed HMMs: ", ', '.join(installed_hmm_names))
        else:
            print("Scanning with pre-installed HMMs: {}".format(installed_hmm_names[0]))
            logging.info("Scanning with pre-installed HMMs: {}".format(installed_hmm_names[0]))
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
                logging.info("No installation_dir specified for db " + hmm_db)
                continue

            #if we're in meta mode, we don't want to keep all that shit in memory
            #and the hmmscan function will write a file for each DB and each protein file
            #because they're huge
            if not meta:
                results_dataframes_dict = hmmscan(protein_dict, db_hmms, threads, hmmscan_options, hmm_db)
            else:
                hmmscan(protein_dict, db_hmms, threads, hmmscan_options, hmm_db)

            if args.write_seqs:
                extract_sequences(results_dataframes_dict, outdir)

            if not meta:
                db_results_df = pd.concat([results_dataframes_dict[key] for key in results_dataframes_dict.keys()])
                db_results_df.to_csv(os.path.join(outdir,hmm_db + '_hits_df.tsv'), sep='\t', index=False)
    time_printout  = "Process took {} seconds.".format(time.time()-t1)
    print(time_printout)
    logging.info(time_printout)

if __name__ == "__main__":
    main()


#TODO:
"""
- Multithread extract_sequences
"""