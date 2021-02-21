#!/usr/bin/env python
from __future__ import print_function
import sys
import numpy as np
import matplotlib.pyplot as plt

exname = 'fracture_3d_example'

def plot(fldname,time):

    # get time data
    Tdata = []
    with open(exname+".time",'r') as fl:
        lines = fl.readlines()
    for line in lines:
        Tdata.append(float(line.strip().split()[-1]))
    Tdata = np.array(Tdata)
        
    tidx = np.argmin(abs(Tdata-time))
    print("time: {}".format(Tdata[tidx]))
        
    # get space data
    Xdata = []
    Zdata = []
    with open(exname+".coord",'r') as fl:
        lines = fl.readlines()
    for line in lines:
        c = line.strip().split()
        Xdata.append(float(c[0]))
        Zdata.append(float(c[2]))
    Xdata = np.array(Xdata)
    Zdata = np.array(Zdata)
        
    # get field data
    Vdata = []
    with open(exname+"-DataFiles/"+fldname+".out",'r') as fl:
        lines = fl.readlines()
    for i,line in enumerate(lines):
        if i==tidx:
            Vdata.append([float(i) for i in line.strip().split()])
    Vdata = np.array(Vdata[0])

    print(Xdata.shape,Zdata.shape,Vdata.shape)
    
    # structure data for plot
    p = Xdata.argsort()
    Xdata = Xdata[p]
    Zdata = Zdata[p]
    Vdata = Vdata[p]

    Xr = Xdata[::-1]
    Xmin = Xr[-1]
    lng = len(Xr) - np.where(Xr==Xmin)[0][0]
    Xdata.shape = [-1,lng]
    Zdata.shape = [-1,lng]
    Vdata.shape = [-1,lng]

    soi = np.argsort(Zdata, axis=1)
    sti = np.indices(Zdata.shape)
    Xdata = Xdata[sti[0], soi]
    Zdata = Zdata[sti[0], soi]
    Vdata = Vdata[sti[0], soi]
    
    # plot
    fig = plt.figure()
    ax = fig.add_subplot(111)
    pc = ax.pcolor(Xdata,Zdata,Vdata)
    cbar = fig.colorbar(pc)
    cbar.set_label(fldname)
    ax.set_xlabel('x')
    ax.set_ylabel('t')


# -----------------------------------------------------------------------------
if __name__ == "__main__":

    if len(sys.argv) not in [3]:
        sys.exit('Missing argument! usage: ./fracture_2d_example_plot.py '
                 + 'fieldname time'
                 + '(options for fieldname: cohesion_0 cohesion_1 cohesion_2 '
                 + 'top_disp_0 top_disp_1 top_disp_2)')

    plot(str(sys.argv[1]),float(sys.argv[2]))
    plt.show()
