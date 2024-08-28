import os
import sys
import argparse
import pandas as pd
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

Result = collections.namedtuple("Result", ["sequence_id", "hmm_name", "bitscore", "evalue",
                                            "env_from", "env_to"])

def get_results_attributes(result):
    bitscore = result.bitscore
    evalue = result.evalue
    cog = result.hmm_name
    query = result.sequence_id
    env_from = result.env_from
    env_to = result.env_to
    return [query, cog, bitscore, evalue, env_from, env_to]

def nucsearch(nuc_dict, hmms, threads, options, db_name=None):
    results_dataframes = {}
    nucsearch_kwargs = define_kwargs(options)

    for fasta_file, sequences in tqdm(nuc_dict.items()):
        results = []
        
        for hits in pyhmmer.nhmmer(hmms, sequences, cpus=threads, **nucsearch_kwargs):
            cog = hits.query_name.decode()
            for hit in hits:
                if hit.included:
                    hit_name = hit.name.decode()
                    full_bitscore = hit.score 
                    full_evalue = hit.evalue
                    results.append(Result(hit_name, cog, full_bitscore, full_evalue, 
                                          hit.env_from, hit.env_to))

        result_df = pd.DataFrame(list(map(get_results_attributes, results)), 
                                 columns=["sequence_id", "hmm_name", "bitscore", "evalue", 
                                          "env_from", "env_to"])
        
        if not options['meta']:
            results_dataframes[fasta_file] = result_df
        else:
            basename_fasta = os.path.basename(fasta_file)
            output_filename = f"{basename_fasta}_{db_name}_results.tsv" if db_name else f"{basename_fasta}_results.tsv"
            result_df.to_csv(os.path.join(options['outdir'], output_filename), sep='\t', index=False)

    return results_dataframes if not options['meta'] else None

def parse_nuc_input(nuc_in, threads):
    print("Parsing nucleotide input sequences...")
    nuc_dict = {}
    
    if os.path.isdir(nuc_in):
        if not os.listdir(nuc_in):
            print("nuc_in directory is empty.")
            logging.info("nuc_in directory is empty.")
            sys.exit(1)

        for fasta_file in tqdm([os.path.join(nuc_in, x) for x in os.listdir(nuc_in)]):
            with pyhmmer.easel.SequenceFile(fasta_file, digital=True, alphabet=pyhmmer.easel.Alphabet.dna()) as seq_file:
                sequences = seq_file.read_block()
            nuc_dict[fasta_file] = sequences
    elif os.path.isfile(nuc_in):
        if os.path.getsize(nuc_in) == 0:
            print("nuc_in file is empty.")
            logging.info("nuc_in file is empty.")
            sys.exit(1)
        with pyhmmer.easel.SequenceFile(nuc_in, digital=True, alphabet=pyhmmer.easel.Alphabet.dna()) as seq_file:
            sequences = seq_file.read_block()
        nuc_dict[nuc_in] = sequences
    else:
        print("Invalid input for nuc_in.")
        logging.info("Invalid input for nuc_in.")
        sys.exit(1)
    
    return nuc_dict

def parse_hmms(hmm_in):
    hmms = []
    print("Parsing HMMs...")
    if os.path.isdir(hmm_in):
        if not os.listdir(hmm_in):
            print("hmm_in directory is empty.")
            logging.info('hmm_in directory is empty.')
            sys.exit(1)

        hmm_files = list(filter(lambda x: x.endswith(('.hmm', '.HMM')), os.listdir(hmm_in)))
        hmm_paths = [os.path.join(hmm_in, hmm_file) for hmm_file in hmm_files]
        
        hmms = list(tqdm(map(lambda x: pyhmmer.plan7.HMMFile(x).read(), hmm_paths)))

    elif os.path.isfile(hmm_in):
        if os.path.getsize(hmm_in) == 0:
            print("hmm_in file is empty.")
            logging.info('hmm_in file is empty.')
            sys.exit(1)
        with pyhmmer.plan7.HMMFile(hmm_in) as hmm_file:
            hmms = list(hmm_file)
    else:
        print("Invalid HMM input.")
        logging.info("Invalid HMM input.")
        sys.exit(1)

    print("HMMs parsed.")
    return hmms

def define_kwargs(options):
    kwargs = {}
    
    if options['cut_ga']:
        kwargs['cut_ga'] = True
    elif options['cut_nc']:
        kwargs['cut_nc'] = True
    elif options['cut_tc']:
        kwargs['cut_tc'] = True

    if options['bitscore'] is not None:
        kwargs['T'] = float(options['bitscore'])
    if options['evalue'] is not None:
        kwargs['E'] = float(options['evalue'])

    return kwargs

def main(args):
    t1 = time.time()
    nuc_in = args.nuc_in
    hmm_in = args.hmm_in
    outdir = args.outdir
    installed_hmms = args.installed_hmms

    log_file_path = os.path.join(outdir, 'astra_nucsearch_log.txt')

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    logging.basicConfig(filename=log_file_path, level=logging.INFO,
                        format='%(asctime)s %(levelname)s: %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')

    options = {
        "meta": args.meta,
        "cut_ga": args.cut_ga,
        "cut_nc": args.cut_nc,
        "cut_tc": args.cut_tc,
        "evalue": args.evalue,
        "bitscore": args.bitscore,
        "outdir": outdir,
    }

    if hmm_in is None and installed_hmms is None:
        print("Either a user-provided or pre-installed HMM database must be specified.")
        logging.error("No HMM database specified.")
        sys.exit(1)

    nuc_dict = parse_nuc_input(nuc_in, args.threads)
    
    if hmm_in is not None:
        print("Searching with user-provided HMM(s)...")
        logging.info("Searching with user-provided HMM(s)...")
        user_hmms = parse_hmms(args.hmm_in)
        results_dataframes_dict = nucsearch(nuc_dict, user_hmms, args.threads, options)
        
        if not options['meta']:
            all_results_df = pd.concat([results_dataframes_dict[key] for key in results_dataframes_dict.keys()])
            all_results_df.to_csv(os.path.join(outdir,'user_hmms_hits_df.tsv'), sep='\t', index=False)
    
    if installed_hmms is not None:
        parsed_json = initialize.load_config()
        installed_hmm_names = installed_hmms.split(',') if ',' in installed_hmms else [installed_hmms]

        print(f"Searching with pre-installed HMMs: {', '.join(installed_hmm_names)}")
        logging.info(f"Searching with pre-installed HMMs: {', '.join(installed_hmm_names)}")

        for hmm_db in installed_hmm_names:
            installed_hmm_in = next((item for item in parsed_json['db_urls'] if item["name"] == hmm_db), None)
            if installed_hmm_in is not None:
                installation_dir = installed_hmm_in['installation_dir']
                db_hmms = parse_hmms(installation_dir)
            else:
                print(f"No installation_dir specified for db {hmm_db}")
                logging.info(f"No installation_dir specified for db {hmm_db}")
                continue

            if not options['meta']:
                results_dataframes_dict = nucsearch(nuc_dict, db_hmms, args.threads, options, hmm_db)
                db_results_df = pd.concat([results_dataframes_dict[key] for key in results_dataframes_dict.keys()])
                db_results_df.to_csv(os.path.join(outdir, f"{hmm_db}_hits_df.tsv"), sep='\t', index=False)
            else:
                nucsearch(nuc_dict, db_hmms, args.threads, options, hmm_db)

    time_printout = f"Process took {time.time()-t1} seconds."
    print(time_printout)
    logging.info(time_printout)

if __name__ == "__main__":
    main()
