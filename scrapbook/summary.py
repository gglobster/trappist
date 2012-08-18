#!/usr/bin/env python
# a stacked bar plot 
import numpy as np
import matplotlib.pyplot as plt


N = 6
accepted = (637598, 5955272, 0, 1135420, 6030630, 0)
rejected = (343230, 2320966, 3202894, 540966, 2079084, 2492682)

ind = (1, 2, 3, 4.5, 5.5, 6.5)   
width = 0.8       	 

p1 = plt.bar(ind, accepted, width, color='g')
p2 = plt.bar(ind, rejected, width, color='r', bottom=accepted)

plt.ylabel('Number of reads')
plt.title('Proportion of accepted vs. rejected reads')
plt.xticks((1.4, 2.4, 3.4, 4.9, 5.9, 6.9), 
		   ('PA_L (75)', 'PA_S (45)', 'PA_R', 'PB_L(75)', 'PB_S(45)', 'PB_R') )
plt.legend( (p1[0], p2[0]), ('Accepted', 'Rejected') )
plt.grid(True)
plt.autoscale(enable=False, axis='both', tight=None)

plt.show()
