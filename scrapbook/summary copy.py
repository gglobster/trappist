#!/usr/bin/env python
# a stacked bar plot 
import numpy as np
import matplotlib.pyplot as plt


N = 6
accepted = (995431, 3000994, 0, 1380608, 2767587, 0)
rejected = (421853, 717483, 1094219, 520684, 584721, 885753)

ind = (1, 2, 3, 4.5, 5.5, 6.5)   
width = 0.8       	 

p1 = plt.bar(ind, accepted, width, color='g')
p2 = plt.bar(ind, rejected, width, color='r', bottom=accepted)

plt.ylabel('Number of reads')
plt.title('Proportion of accepted vs. rejected reads')
plt.xticks((1.4, 2.4, 3.4, 4.9, 5.9, 6.9), 
		   ('PAtL(85)', 'PAtS(50)', 'PAR', 'PBtL(85)', 'PBtS(50)', 'PBR') )
plt.legend( (p1[0], p2[0]), ('Accepted', 'Rejected') )
plt.grid(True)
plt.autoscale(enable=False, axis='both', tight=None)

plt.show()
