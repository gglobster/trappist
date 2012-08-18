### Blast Illumina data against a reference DB
### Filter out positive matches to separate file

import datetime
import sys,subprocess
from sys import argv
from Bio.SeqIO.QualityIO import FastqGeneralIterator

source_file = argv[1]
nickname = argv[2]
target_file = argv[3]
bridge_file = argv[4]

print "--- Converting to multifasta ---"
print datetime.datetime.now()

read_count = 0
multifasta = open(target_file, 'w')
tracker = open(bridge_file, 'w')

for title, seq, qual in FastqGeneralIterator(open(source_file)) :
	read_count +=1
	id_string = nickname+"_"+str(read_count)
	multifasta.write(">"+id_string+"\n"+seq+"\n")
	tracker.write(title+"\t"+id_string+"\n")
	if read_count%100000==0:
		print read_count, "reads processed"

print "--- Finished conversion ---"
print datetime.datetime.now()
print read_count, "reads written"

    