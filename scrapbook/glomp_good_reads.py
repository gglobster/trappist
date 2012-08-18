### Use matching reads lists to glomp good reads ###

import re
import numpy
import datetime
from sys import argv

def read_array(filename, dtype, separator='\t'):
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


base_dir = argv[1]

match_dir = base_dir+"match/"
q2a_dir = base_dir+"Q2A/"
mask_dir = base_dir+"mask/"

genome_names = [("Pstutz", "PstutzA", "PstutzLMG"), ("Pfluo1", "Pfluo5", "PfluoSBW")]
bin_names = [("PAtL", "PAtS", "PAR"), ("PBtL", "PBtS", "PBR")]

print "\n--- Fusing lists of chromosome-matching reads per bin ---"
print datetime.datetime.now()

index = 0
c_mreads_by_bin = {}
while index < len(bin_names):
	for name in bin_names[index]:
		print "\n--", name, "--"
		# collect match files from each genome
		bin_data_arrays = []
		for genome in genome_names[index]:
			data_array = numpy.load(match_dir+genome+"/"+name+"_match.npy")
			bin_data_arrays.append(data_array)
			print "\t"+genome, len(data_array), "unique matching reads"
		# take the intersection of the match sets
		inter_first = numpy.intersect1d(bin_data_arrays[0], bin_data_arrays[1])
		inter_then = numpy.intersect1d(inter_first, bin_data_arrays[2])
		print len(inter_then), "present in all reference genomes"
		# take the union of the match sets
		uni_first = numpy.union1d(bin_data_arrays[0], bin_data_arrays[1])
		uni_then = numpy.union1d(uni_first, bin_data_arrays[2])
		print len(uni_then), "matching reads all together"		
		c_mreads_by_bin[name] = uni_then
	index +=1

print "\n--- Fusing lists of phage-matching reads per bin ---"
print datetime.datetime.now()

index = 0
p_mreads_by_bin = {}
while index < len(bin_names):
	for name in bin_names[index]:
		print "\n--", name, "--"
		# collect match files from each genome
		data_array = numpy.load(match_dir+"Pphages/"+name+"_match.npy")
		print "\tPphages", len(data_array), "unique matching reads"	
		p_mreads_by_bin[name] = data_array
	index +=1

print "\n--- Comparing match lists to rescue phage material ---"
print datetime.datetime.now()

all_bin_names = ["PAtL", "PAtS", "PBtL", "PBtS"]
# From now on we're not including the rejected bins PAR and PBR

rescue_list = {}
rejection_list = {}
for name in all_bin_names:
	print "\n--", name, "--"
	rescued = []
	rejected = []
	c_list = c_mreads_by_bin[name]
	p_list = p_mreads_by_bin[name]
	for element in c_list:			# optimize by using numpy array masking
		if element in p_list:
			rescued.append(element)
		else:
			rejected.append(element)
	print "\t"+str(len(rescued)), "to rescue"
	print "\t"+str(len(rejected)), "to reject"
	rescue_list[name] = rescued
	rejection_list[name] = rejected

print "\n--- Separating out accepted vs. rejected read titles ---"
print datetime.datetime.now()

accept_titles = {}
reject_titles = {}
q2a_by_bin = {}
for name in all_bin_names:
	print "\n--", name, "--"
	# read in Q2A files
	q2a_file = q2a_dir+name+"_q2a.txt"
	dtype = numpy.dtype([('title', 'S50'), ('bincode', 'S15')])
	pair_array = read_array(q2a_file, dtype, separator='\t')
	# create masking array 
	mask = numpy.zeros(len(pair_array), bool)
	for item in rejection_list[name]:
		mask[item-1] = True 	# True means must reject!
	bin_reject_titles = pair_array[mask]
	print "\t"+str(len(bin_reject_titles)), "read titles to reject"
	reject_titles[name] = bin_reject_titles
	# now take the inverse set
	inv_mask = numpy.invert(mask)
	bin_accept_titles = pair_array[inv_mask]
	print "\t"+str(len(bin_accept_titles)), "read titles to accept"
	accept_titles[name] = bin_accept_titles
	# save inverted mask to file (where True means keep, False means reject)
	numpy.save(mask_dir+name+"_mask.npy", inv_mask)
# write to files for future reference
#numpy.save(q2a_dir+"titles_accept.npy", [(accept_titles)])
#numpy.save(q2a_dir+"titles_reject.npy", [(reject_titles)])
 
print datetime.datetime.now()
print "DONE!"