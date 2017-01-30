"""
Collection of useful functions

@author Kimberly Insigne

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
