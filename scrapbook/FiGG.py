"""Iterate over 2 Fastq records as string tuples.

Modified from Bio.SeqIO.QualityIO import FastqGeneralIterator to return
two reads at a time.

"""

from sys import argv

handle = argv[1]

def FiGG(handle):
	handle_readline = handle.readline
	while True:
		line = handle_readline()
		if line == "" : return
		if line[0] == "@":
			break
	while True:
		title_lines = []
		seq_strings = []
		quality_strings = []
		count = 0
		while count < 2 :
			if line[0] != "@":
				raise ValueError("Records in Fastq files should start with '@' character")
			title_line = line[1:].rstrip()
			title_lines.append(title_line)
			seq_string = handle_readline().rstrip()
			seq_strings.append(seq_string)
			while True:
				line = handle_readline()
				if not line:
					raise ValueError("End of file without quality information.")
				if line[0] == "+":
					second_title = line[1:].rstrip()
					if second_title and second_title != title_line:
						raise ValueError("Sequence and quality captions differ.")
					break
				seq_string += line.rstrip() #removes trailing newlines
			if " " in seq_string or "\t" in seq_string:
				raise ValueError("Whitespace is not allowed in the sequence.")
			seq_len = len(seq_string)
			quality_string = handle_readline().rstrip()
			quality_strings.append(quality_string)
			while True:
				line = handle_readline()
				if not line : break #end of file
				if line[0] == "@":
					if len(quality_string) >= seq_len:
						break
				quality_string += line.rstrip()
			if seq_len != len(quality_string):
				raise ValueError("Lengths of sequence and quality values differs "
								 " for %s (%i and %i)." \
								 % (title_line, seq_len, len(quality_string)))
			count +=1
		yield (title_lines, seq_strings, quality_strings)
		if not line : return #StopIteration at end of file
	assert False, "Should not reach this line"
	
count = 0

for duple in FiGG(open(handle)):
	count +=1
	print count, duple[0], duple[1]
	if count ==10 :
		break