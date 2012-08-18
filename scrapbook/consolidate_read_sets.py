### Consolidate sets of accepted reads for assembly ###

import numpy
from sys import argv
from Bio.SeqIO.QualityIO import FastqGeneralIterator

bin_names = ["PA_L", "PA_S", "PB_L", "PB_S"]

source_dir = "data/trimmed/"
mask_dir = "data/tracking/mask/"
destin_dir = "data/assembly/"

source_files = ["PA_1_trimmed_85.txt", 
				"PA_1_trimmed_50.txt",
				"PB_1_trimmed_85.txt",
				"PB_1_trimmed_50.txt"
				]

for name in bin_names:
	source_file = source_dir+source_files[bin_names.index(name)]
	out_file = open(destin_dir+name+"_reads_ok.txt", 'w')
	print "--", name, "--"
	counter = 0
	tru_count = 0
	mask = numpy.load(mask_dir+name+"_mask.npy")
	for title, seq, qual in FastqGeneralIterator(open(source_file)) :
		if mask[counter] == True:
			out_file.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
			tru_count +=1
		else:
			pass
		counter +=1
 		if counter%100000==0: 
 			print "\t"+str(tru_count), "reads selected of", str(counter), "reads processed"
	out_file.close()
print "DONE!"
