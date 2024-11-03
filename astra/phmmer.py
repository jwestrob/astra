import collections
import os
import sys
import time

import pandas as pd
import pyhmmer
import pyhmmer.plan7

Result = collections.namedtuple("Result", ["sequence_id", "evalue", 
                                      "env_from", "env_to", "bitscore"])
def get_results_attributes(result):
    bitscore = result.bitscore
    c_evalue = result.evalue
    query = result.sequence_id
    env_from = result.env_from
    env_to = result.env_to
    return [query, c_evalue, env_from, env_to, bitscore]


def phmmer(query, sequence_db, threads):

	search_protocol = pyhmmer.hmmer.phmmer(query, sequence_db, cpus=threads)

	search_results = list(search_protocol)

	results = []
	for hit in search_results[0]:
	    if hit.included:
	        hit_name = hit.name.decode()
	        results.append(Result(hit_name, hit.evalue, 
	                hit.domains[0].env_from, hit.domains[0].env_to, hit.score))
	return results

def setup_logging(outdir: str) -> None:
    """Setup logging configuration"""
    log_file = os.path.join(outdir, 'astra_phmmer.log')
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )

def parse_args() -> argparse.Namespace:
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description="ASTRA phmmer search tool")
    
    parser.add_argument('-q', '--query', 
                        required=True,
                        help='Query sequence file in FASTA format')
    
    parser.add_argument('-d', '--database',
                        required=True,
                        help='Target sequence database in FASTA format')
    
    parser.add_argument('-o', '--outdir',
                        default='phmmer_results',
                        help='Output directory (default: phmmer_results)')
    
    parser.add_argument('-t', '--threads',
                        type=int,
                        default=1,
                        help='Number of CPU threads to use (default: 1)')
    
    parser.add_argument('-E', '--evalue',
                        type=float,
                        default=10.0,
                        help='E-value threshold for reporting hits (default: 10.0)')
    
    parser.add_argument('-T', '--score',
                        type=float,
                        help='Bit score threshold for reporting hits')
    
    parser.add_argument('--domE',
                        type=float,
                        help='Domain E-value threshold')
    
    parser.add_argument('--domT',
                        type=float,
                        help='Domain bit score threshold')
    
    parser.add_argument('--write-fasta',
                        action='store_true',
                        help='Write matching sequences to FASTA files')
    
    return parser.parse_args()

def main():
    args = parse_args()
    
    # Create output directory
    os.makedirs(args.outdir, exist_ok=True)
    
    # Setup logging
    setup_logging(args.outdir)
    logging.info(f"Starting phmmer search with parameters: {vars(args)}")
    
    t_start = time.time()
    
    # Parse sequences
    logging.info("Reading query sequences...")
    with pyhmmer.easel.SequenceFile(args.query, digital=True) as query_reader:
        query = query_reader.read_block()
    
    logging.info("Reading target database...")
    with pyhmmer.easel.SequenceFile(args.database, digital=True) as seq_file:
        sequence_db = seq_file.read_block()
    
    # Run search
    logging.info("Running phmmer search...")
    results = phmmer(query, sequence_db, args.threads)
    
    # Process results
    logging.info("Processing results...")
    results_df = pd.DataFrame(
        list(map(get_results_attributes, results)),
        columns=["sequence_id", "evalue", "env_from", "env_to", "bitscore"]
    ).sort_values(by=['sequence_id', 'evalue'])
    
    # Write results
    out_file = os.path.join(args.outdir, 'phmmer_results.tsv')
    results_df.to_csv(out_file, sep='\t', index=False)
    logging.info(f"Results written to {out_file}")
    
    if args.write_fasta:
        logging.info("Writing matching sequences to FASTA...")
        fasta_dir = os.path.join(args.outdir, 'matching_sequences')
        os.makedirs(fasta_dir, exist_ok=True)
        
        with open(os.path.join(fasta_dir, 'matches.fasta'), 'w') as f:
            for seq in sequence_db:
                if seq.name.decode() in results_df['sequence_id'].values:
                    f.write(f">{seq.name.decode()}\n{seq.sequence}\n")
    
    t_end = time.time()
    logging.info(f"phmmer search completed in {t_end - t_start:.2f} seconds")

if __name__ == "__main__":
    main()
