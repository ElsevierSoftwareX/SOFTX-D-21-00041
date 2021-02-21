#!/usr/bin/env python
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt

exname = 'basic_example'

# get time data
Tdata = []
with open(exname+".time",'r') as fl:
    lines = fl.readlines()
for line in lines:
    Tdata.append(float(line.strip().split()[-1]))

# get space data
Xdata = []
with open(exname+".coord",'r') as fl:
    lines = fl.readlines()
for line in lines:
    Xdata.append(float(line.strip().split()[0]))

# get field data
Vdata = []
with open(exname+"-DataFiles/cohesion_0.out",'r') as fl:
    lines = fl.readlines()
for line in lines:
    Vdata.append([float(i) for i in line.strip().split()])
Vdata = np.array(Vdata)

# plot
fig = plt.figure()
ax = fig.add_subplot(111)
XV, TV = np.meshgrid(Xdata, Tdata)
pc = ax.pcolor(XV,TV,Vdata)
cbar = fig.colorbar(pc)
cbar.set_label('cohesion')
ax.set_xlabel('x')
ax.set_ylabel('t')
plt.show()
