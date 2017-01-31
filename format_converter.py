import argparse

def fasta_reader(filename):

	infile = open(filename)
	seqs = {}

	for line in infile:
		name = line.strip()[1:] # remove leading '>'
		seq = next(infile).strip()
		seqs[name] = seq

	return seqs

def tab_reader(filename):
	infile = open(filename)
	seqs = {}

	for line in infile:
		name, seq = line.strip().split('\t')
		seqs[name] = seq
		seqs[name] = seq

	return seqs

def csv_reader(filename):
	infile = open(filename)
	seqs = {}

	for line in infile:
		name, seq = line.strip().split(',')
		seqs[name] = seq

	return seqs

def fasta_writer(filename, seqs):
	outfile = open(filename, 'w')
	
	for x in seqs:
		outfile.write('>'+x+'\n'+seqs[x]+'\n')

	outfile.close()

def tab_writer(filename, seqs):
	outfile = open(filename, 'w')
	
	for x in seqs:
		outfile.write(x+'\t'+seqs[x]+'\n')

	outfile.close()

def csv_writer(filename, seqs):
	outfile = open(filename, 'w')
	
	for x in seqs:
		outfile.write(x+','+seqs[x]+'\n')

	outfile.close()

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Switch between either csv,tab or fasta formats')
	parser.add_argument('filename', help='Name of input file')
	parser.add_argument('input_type', help='Type of input file')
	parser.add_argument('output_type', help='Type of output file')
	parser.add_argument('output_name', help='Name of output file')

	args = parser.parse_args()

	input_filename = args.filename
	output_type = args.output_type
	input_type = args.input_type
	output_name = args.output_name

	if input_type == 'fasta':
		seqs = fasta_reader(input_filename)
	elif input_type == 'csv':
		seqs = csv_reader(input_filename)
	elif input_type == 'tab':
		seqs = tab_reader(input_filename)
	else:
		raise ValueError('Input type must be csv, tab or fasta')

	if output_type == 'fasta':
		fasta_writer(output_name, seqs)
	elif output_type == 'csv':
		csv_writer(output_name, seqs)
	elif output_type == 'tab':
		tab_writer(output_name, seqs)
	else:
		raise ValueError('Output type must be csv, tab or fasta')

