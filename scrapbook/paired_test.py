### Test for paired ends ###

import re
from sys import argv
from Bio.SeqIO.QualityIO import FastqGeneralIterator

source_file = argv[1]

print "-- Testing for paired ends --"
counter = 0
count1 = 0
count2 = 0
for title, seq, qual in FastqGeneralIterator(open(source_file)) :
	if int(title[-1]) is 1:
		count1 +=1
	elif int(title[-1]) is 2:
		count2 +=1
	counter +=1
	if counter%100000==0: 
		print len(title), title
		print "\t"+str(count1), "fwd reads", "and", str(count2), "reverse reads out of", str(counter), "reads processed"

print "Total:", counter, "reads in dataset"
