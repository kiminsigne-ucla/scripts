'''
This script takes as input a fastq file of reads that have the sequencing index as the last 6bp and a tab
separated text file containing sample names in the first column and second column with the index (no header).
It separates the input fastq file into separate fastq files based on their index. This script only works with
paired-end reads files.
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
	parser.add_argument('read1_file', help='.fastq file of R1')
	parser.add_argument('read2_file', help='.fastq file of R2')
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
	index_to_handle = {}
	for sample in samples:
		index_to_handle[sample+'_read1'] = open(sample+'_read1.fastq', 'w')
		index_to_handle[sample+'_read2'] = open(sample+'_read2.fastq', 'w')

	
	bad_index_1 = open('bad_index_1.fastq', 'w')
	bad_index_3 = open('bad_index_3.fastq', 'w')

	read1_file = open(args.read1_file)
	read2_file = open(args.read2_file)
	index_reads_file = open(args.index_reads_file)

	count = 0

	while True:
		read1_info = list(islice(read1_file, 4))
		read2_info = list(islice(read2_file, 4))
		index_info = list(islice(index_reads_file, 4))
		if not read1_info:
			break

		if count % 10000000 == 0:
			print count, '...'

		index = index_info[1].strip()[:args.index_length]

		# default
		filehandle1 = bad_index_1
		filehandle2 = bad_index_3

		# if index is within 1bp of known index
		for x in indices:
			if Levenshtein.distance(x, index) <= 1:
				index_match = x
				filehandle1 = index_to_handle[indices[index_match] + '_read1']
				filehandle2 = index_to_handle[indices[index_match] + '_read2']
			
		filehandle1.write(read1_info[0])
		filehandle1.write(read1_info[1])
		filehandle1.write(read1_info[2])
		filehandle1.write(read1_info[3])
		filehandle2.write(read2_info[0])
		filehandle2.write(read2_info[1])
		filehandle2.write(read2_info[2])
		filehandle2.write(read2_info[3])

		count += 1

	for x in index_to_handle:
		index_to_handle[x].close()


