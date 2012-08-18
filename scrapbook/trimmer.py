### Trim reads

import pprint
from sys import argv
from Bio.SeqIO.QualityIO import FastqGeneralIterator

def trim_illumina(source_file, bin_files, bin_sizes, clip):
	"""Trim reads based on Illumina's PHRED quality scores.
	
	Trimmed reads are binned by size.
	Reads are loaded and processed as plain strings to avoid extra overhead from SeqIO.
	Extended from http://news.open-bio.org/news/2010/04/illumina-q2-trim-fastq/
	
	"""
	# Prep bin files
	untrimmed_out = open(bin_files['untrimmed'], "w")
	#small_bin_out = open(bin_files['small_bin'], "w")
	large_bin_out = open(bin_files['large_bin'], "w")
	rejected_out = open(bin_files['rejected'], "w")
	# Initialize counters
	total_count = 0
	untrimmed_count = 0
	#small_bin_count = 0
	large_bin_count = 0
	rejected_count = 0
	# Cycle through reads
	for title, seq, qual in FastqGeneralIterator(open(source_file)) :
		#Find the location of the first "B" (PHRED quality 2)
		trim = qual.find("B")
		if trim == -1:
			trim = bin_sizes['length']-clip[1]
			untrimmed_out.write("@%s\n%s\n+\n%s\n" % (title, seq[clip[0]:trim], qual[clip[0]:trim]))
			untrimmed_count +=1
		elif trim < bin_sizes['mid']: # use 'min' if using two trimmed bin sizes
			rejected_out.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
			rejected_count +=1
# 		elif trim < bin_sizes['mid']:
# 			trim = bin_sizes['min']
# 			small_bin_out.write("@%s\n%s\n+\n%s\n" % (title, seq[clip[0]:trim], qual[clip[0]:trim]))
# 			small_bin_count +=1
		elif trim < bin_sizes['length']:
			trim = bin_sizes['mid']
			large_bin_out.write("@%s\n%s\n+\n%s\n" % (title, seq[clip[0]:trim], qual[clip[0]:trim]))
			large_bin_count +=1
		total_count +=1
		if total_count%500000==0:
			print total_count, "reads processed"
	# Close bin files
	untrimmed_out.close()
	#small_bin_out.close()
	large_bin_out.close()
	rejected_out.close()
	# Summarize count results
	counts_dict = {'untrimmed': untrimmed_count, 
				   'rejected': rejected_count, 
				   #'small': small_bin_count, 
				   'large': large_bin_count, 
				   'total': total_count
				   }
	return counts_dict

source_file = argv[1]
nickname = argv[2]

bin_sizes = {'length': 101, 
			 'min': 35, 
			 'mid': 55
			 }

clip = [5, 11]

untrimmed = nickname+"_trimmed_"+str(bin_sizes['length']-clip[0]-clip[1])+".txt"
#small_bin = nickname+"_trimmed_"+str(bin_sizes['min']-clip[0])+".txt"
large_bin = nickname+"_trimmed_"+str(bin_sizes['mid']-clip[0])+".txt"
rejected = nickname+"_rejected.txt"

bin_files = {'untrimmed': untrimmed, 
			 #'small_bin': small_bin, 
			 'large_bin': large_bin, 
			 'rejected': rejected
			 }

# Run trimmer
counts = trim_illumina(source_file, bin_files, bin_sizes, clip)
printout = pprint.PrettyPrinter(indent=4)
printout.pprint(counts)