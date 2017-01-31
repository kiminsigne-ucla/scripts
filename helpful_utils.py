"""
Collection of useful functions

@author Kimberly Insigne
kiminsigne@gmail.com

"""

from itertools import islice
import string


def fasta_reader(filename):
	"""
	Input: str, name of file
	Output: dictionary, key = header, value = sequence
	"""

	seqs = {}
	with open(filename) as infile:
		entry = list(islice(infile, 2))
		while entry:
			# grab two lines from file at a time, strip \n
			header, seq = map(string.strip, entry)
			# strip '>' from header, first character
			seqs[header[1:]] = seq.upper()
			entry = list(islice(infile, 2))

		return seqs


def reverse_complement(seq):
	"""
	Return the reverse complement of a nucleotide string
	"""
	complement = {'A': 'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
	
	rc = ''.join([complement[nt] for nt in seq[::-1]])
	return rc


def mutate_1_bp(seq):

	nts = set(['A', 'C', 'G', 'T'])
	mutants = []

	for i in range(len(seq)):
		mod = list(seq)
		diff = nts - set(mod[i])
		for nt in diff:
			mod[i] = nt
			mutants.append(''.join(mod))

	return mutants
