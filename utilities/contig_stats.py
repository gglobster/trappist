## script to compile basic stats about sets of contigs

import re
from sys import argv
from libs.common import load_multifasta, from_dir
import matplotlib.pyplot as plt
import numpy as np

data_dir = "data/"+argv[1]

filenames = from_dir(data_dir, re.compile(r'.*\.fas.*'))

ctg_ns = []
n50s = []

for filename in filenames:
    # load contigs from file
    contig_list = load_multifasta(data_dir+"/"+filename)
    # count contigs
    ctg_count = len(contig_list)
    if ctg_count < 200:
        ctg_ns.append(ctg_count)
    else:
        ctg_ns.append(200)

    # sort contig list by size
    contig_list.sort(key=len)
    contig_list.reverse()

    # count full sequence length
    full_seq_length = 0
    for contig in contig_list:
        full_seq_length += len(contig.seq)

    # calculate N50
    calc_n50 = 0
    low_n50 = 0
    for contig in contig_list:
        calc_n50 += len(contig.seq)
        low_n50 = len(contig.seq)
        if calc_n50 > full_seq_length/2:
            if low_n50 > 300000:
                n50s.append(300)
            else:
                n50s.append(low_n50/1000)
            break

    print filename,
    print ctg_count,
    print full_seq_length,
    print low_n50

# plot metrics summary

ind = np.arange(len(ctg_ns))    # the x locations for the groups
width = 0.2                     # the width of the bars

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax2 = ax1.twinx()

ax1.bar(ind+width, ctg_ns, width, color='red')
ax2.bar(ind+width*2, n50s, width, color='blue')

ax1.set_title('Basic assembly statistics')
ax1.set_ylabel('Number of contigs')
ax2.set_ylabel('N50 contig size (kb)')

plt.show()



