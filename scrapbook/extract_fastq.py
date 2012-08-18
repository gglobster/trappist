import sys
from sys import argv
from Bio.SeqIO.QualityIO import FastqGeneralIterator

source_file = argv[1]
number = int(argv[2])

count = 0
for title, seq, qual in FastqGeneralIterator(open(source_file)) :
	if count < number:
		print seq
		count +=1
	else: 
		sys.exit()
