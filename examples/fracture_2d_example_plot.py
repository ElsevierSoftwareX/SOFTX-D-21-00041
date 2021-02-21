#!/usr/bin/env python
from __future__ import print_function
import sys
import numpy as np
import matplotlib.pyplot as plt

exname = 'fracture_2d_example'

def plot(fldname,ax):

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
    with open(exname+"-DataFiles/"+fldname+".out",'r') as fl:
        lines = fl.readlines()
    for line in lines:
        Vdata.append([float(i) for i in line.strip().split()])
    Vdata = np.array(Vdata)

    # plot
    XV, TV = np.meshgrid(Xdata, Tdata)
    pc = ax.pcolor(XV,TV,Vdata)
    return pc

# -----------------------------------------------------------------------------
if __name__ == "__main__":

    if len(sys.argv) not in [2]:
        sys.exit('Missing argument! usage: ./fracture_2d_example_plot.py '
                 + 'fieldname '
                 + '(options: cohesion_0 cohesion_1 top_disp_0 top_disp_1 '
                 + 'bot_disp_0 bot_disp_1 tau_max)')

    fldname=str(sys.argv[1])
        
    fig = plt.figure()
    ax = fig.add_subplot(111)
    pc = plot(fldname,ax)
    cbar = fig.colorbar(pc)
    cbar.set_label(fldname)
    ax.set_xlabel('x')
    ax.set_ylabel('t')
    plt.show()
