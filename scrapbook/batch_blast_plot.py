### Map & plot batch Blast results against a reference ###

import os
import re
import numpy
import datetime
from sys import argv
import matplotlib.pyplot as pylot

def read_array(filename, dtype, separator='\t'):
	# From Numpy cookbook
    """ Read a file with an arbitrary number of columns.
        The type of data in each column is arbitrary
        It will be cast to the given dtype at runtime
    """
    cast = numpy.cast
    data = [[] for dummy in xrange(len(dtype))]
    for line in open(filename, 'r'):
        fields = line.strip().split(separator)
        for i, number in enumerate(fields):
            data[i].append(number)
    for i in xrange(len(dtype)):
        data[i] = cast[dtype[i]](data[i])
    return numpy.rec.array(data, dtype=dtype)


# define inputs
blast_out_dir = argv[1]
bin_names = (argv[2], argv[3], argv[4])
ref_name = argv[5]
phi_name = argv[6]
match_dir = argv[7]

data_descr = numpy.dtype([('query', 'S16'), 
					  	  ('dbhit', 'S32'), 
					  	  ('idp', 'float'), 
					  	  ('mlen', 'uint8'),
					  	  ('mms', 'uint8'),
					  	  ('gaps', 'uint8'), 
					 	  ('q_start', 'uint8'), 
						  ('q_end', 'uint8'),
						  ('r_start', 'uint8'), 
						  ('r_end', 'uint8'),
						  ('evalue', 'S5'),
						  ('bitscore', 'float'), 
						  ])

### Main script ### 

print "--- Consolidating data series ---"
print datetime.datetime.now()

# compile data series
all_series = []
for name in bin_names:
	index = 1
	print "\tconsolidating", name, "data"
	bin_arrays =[]
	while os.path.isfile(blast_out_dir+name+"/"+name+"_"+str(index)+"_blast.out"):
		infile = blast_out_dir+name+"/"+name+"_"+str(index)+"_blast.out"
		rec_array = read_array(infile, data_descr)
		bin_arrays.append(rec_array)
		index +=1
	print "\t\t"+str(len(bin_arrays)), "arrays in", name
	series = numpy.hstack(bin_arrays)
	print "\t\t"+str(len(series)), "records in the series"
	all_series.append(series)

print len(all_series), "data series consolidated"

print "--- Parsing and binning matches, positions on the reference ---"
print datetime.datetime.now()

# parse and bin
bin_index = 0
binned_pos = []
for series in all_series:
	positions = []
	match_read = []
	for row in series:	
		# collect match read info while we're at it
		# use regex to extract query index 
 		query_pattern = re.compile(r'\w*\_(\d*)')
 		query_match = query_pattern.match(row[0])
 		query_index = int(query_match.group(1))
		match_read.append(query_index)
		# use regex to extract ref coords 
		#ref_pattern = re.compile(r'\w*\_(\d*)\_(\d*)') #if we want both
		ref_pattern = re.compile(r'\w*\_\d*\_\d*\_(\d*)') # but here we don't
		ref_match = ref_pattern.match(row[1])
		#ref_coord1 = int(ref_match.group(1))
		#ref_coord2 = int(ref_match.group(2))
		ref_pos = int(ref_match.group(1))
		pos_scaled = ref_pos/100				# adjust to db segment length
		positions.append(pos_scaled)
	# uniquify the match read array
	unique_matches = numpy.unique(match_read)
	print "\t"+str(len(unique_matches)), "unique matching reads in", bin_names[bin_index]
	# write to file for future reference
	numpy.save(match_dir+bin_names[bin_index]+"_match.npy", unique_matches)
	# now count ocurrences per position
	pos_np = numpy.array(positions)
	binned = numpy.bincount(pos_np)	
	binned_pos.append(binned)
	bin_index +=1

# determine what order to use so the plot will look good
series_index = 0
averages = []
for binned in binned_pos:
	pos_count_average = numpy.average(binned)
	#print "\t"+bin_names[series_index], ":", len(all_series[series_index]), "matches"
	averages.append((pos_count_average, series_index))
	series_index +=1
averages.sort()
averages.reverse()
order_indices = []
for pair in averages:
	order_indices.append(pair[1])
	
print "--- Plotting data series ---"
print datetime.datetime.now()

# plot the lot
pylot.autoscale(enable=True, axis='both', tight=True)
pylot.xlabel('Position on the chromosome (/100)')
pylot.ylabel('Number of matches (includes multiples)')
pylot.title(phi_name+' matches to '+ref_name)
pylot.grid(True)
for index in order_indices:
	pylot.plot(binned_pos[index], label=bin_names[index]+" ("+str(numpy.sum(binned_pos[index]))+")")
pylot.legend(loc=1)
pylot.show()
print "\tdone, see plot"
print datetime.datetime.now()

