#!/usr/bin/env python
# a stacked bar plot 
import numpy 
import matplotlib.pyplot as plt
from sys import argv

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
    
statsfile = argv[1]
set_name = argv[2]

stats_dtypes = numpy.dtype([('ID', 'uint8'), 
							('length', 'uint16'), 
							('out', 'uint8'),
							('in', 'uint8'),
							('long_cov', 'float'),
							('short1_cov', 'float'),
							('short1_0cov', 'float'),
							('short2_cov', 'float'),
							('short2_0cov', 'float'),
							('long_nb', 'uint8'),
							('short1_nb', 'uint32'),
							('short2_nb', 'uint8')])		
					
ass_stats = read_array(statsfile, stats_dtypes)

lengths = ass_stats['length'] 
coverage = ass_stats['short1_cov']

width = 100       	 

p1 = plt.bar(lengths, coverage, width)

plt.xlabel('Length of contigs')
plt.ylabel('Read coverage')
plt.title(set_name)
plt.grid(True)
plt.autoscale(enable=False, axis='both', tight=False)

plt.show()
