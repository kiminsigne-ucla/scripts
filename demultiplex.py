'''
This script takes as input a fastq file of reads that have the sequencing index as the last 6bp and a tab
separated text file containing sample names in the first column and second column with the index (no header).
It separates the input fastq file into separate fastq files based on their index. This script only works with
single-end read files.
'''

import argparse
from itertools import islice
import Levenshtein

def reverse_complement(seq):
	"""
	Return the reverse complement of a nucleotide string
	"""
	complement = {'A': 'T', 'T':'A', 'C':'G', 'G':'C'}
	
	rc = ''.join([complement[nt] for nt in seq[::-1]])
	return rc

if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser.add_argument('reads_file', help='.fastq file')
	parser.add_argument('index_reads_file', help='.fastq file, where the index is the first 6bp of the read')
	parser.add_argument('index_file', help="""tab separated text file where first column is sample name 
												and second column is index""")
	parser.add_argument('index_length', type=int, help="length of index, assuming it's the first n bases")
	parser.add_argument('-rev', action='store_true', help="""If sequences in index_file are reverse complement 
															relative to index reads file""")
	args = parser.parse_args()		

	# read in index file
	index_file = open(args.index_file)
	indices = {}
	for line in index_file:
		name, index = line.strip().split('\t')
		if args.rev:
			index = reverse_complement(index)
		indices[index] = name

	samples = indices.values()
	index_to_handle = {sample : open(sample+'.fastq', 'w') for sample in samples}
	bad_index = open('bad_index.fastq', 'w')

	reads_file = open(args.reads_file)
	index_reads_file = open(args.index_reads_file)

	count = 0

	while True:
		reads_info = list(islice(reads_file, 4))
		index_info = list(islice(index_reads_file, 4))
		if not reads_info:
			break

		if count % 10000000 == 0:
			print count, '...'

		index = index_info[1].strip()[:args.index_length]

		# default
		filehandle = bad_index

		# if index is within 1bp of known index
		for x in indices:
			if Levenshtein.distance(x, index) <= 1:
				index_match = x
				filehandle = index_to_handle[ indices[index_match] ]
			
		filehandle.write(reads_info[0])
		filehandle.write(reads_info[1])
		filehandle.write(reads_info[2])
		filehandle.write(reads_info[3])

		count += 1

	for x in index_to_handle:
		index_to_handle[x].close()


