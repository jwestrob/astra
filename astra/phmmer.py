import os, sys, pandas as pd
import pyhmmer
import pyhmmer.plan7
import time 
import collections

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

def main(args):
	query_file = args.query_seqs
	database_file = args.subject_seqs
	threads = args.threads

	outdir = args.outdir

	# Check if the output directory already exists
	if not os.path.exists(outdir):
		os.makedirs(outdir)

	#Parse query sequences
	with pyhmmer.easel.SequenceFile(query_file, digital=True) as query_reader:
		query = query_reader.read_block()

	#Parse subject sequences
	with pyhmmer.easel.SequenceFile(database_file, digital=True) as seq_file:
		sequence_db = seq_file.read_block()

	results = phmmer(query, sequence_db, threads)

	results_df = pd.DataFrame(
					list(map(get_results_attributes, results)), 
					columns=["sequence_id", "evalue", "env_from", "env_to", "bitscore"]).sort_values(by='sequence_id')

	results_df.to_csv(os.path.join(outdir, 'phmmer_out.tsv'), sep='\t', index=False)


if __name__ == "__main__":
	main(args)