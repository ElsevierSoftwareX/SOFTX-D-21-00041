#!/usr/bin/env python
from __future__ import print_function
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from fracture_2d_example_plot import plot as f2dplot
exname = 'fracture_2d_example'

params = {#'backend': 'ps',
        'text.latex.preamble': [
            #r'\usepackage{gensymb}',
            #r'\usepackage[utf8x]{inputenc}',
            #r'\usepackage[T1]{fontenc}',
            r'\usepackage{amsmath}',
            r'\usepackage{bm}',
            r'\usepackage{amssymb}',
            #r'\usepackage[cm]{sfmath}',
            #r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
            #r'\sansmath'  
        ],
        'font.size': 7, # was 10
        'figure.titlesize' : 8,
        'axes.labelsize': 6, # fontsize for x and y labels (was 10)
        'axes.titlesize': 8,
        'legend.fontsize': 7,#8, # was 10
        'legend.frameon':False,
        'legend.columnspacing' : 1,
        'legend.numpoints': 1,
        'legend.scatterpoints': 1,
        'legend.handlelength': 1,
        'lines.linewidth' : 0.5,
        'axes.linewidth' : 0.5,
        'xtick.labelsize': 6,
        'ytick.labelsize': 6,
        'xtick.major.size' : 2,
        'ytick.major.size' : 2,
        'xtick.minor.size' : 2,
        'ytick.minor.size' : 2,
        'xtick.major.width': 0.5,
        'ytick.major.width': 0.5,
        'xtick.minor.width': 0.5,
        'ytick.minor.width': 0.5,
        #'axes.prop_cycle': cycler(color=matlab_colors),
        #'xtick.major.pad'  : 1,
        #'ytick.major.pad'  : 1,
        #'xtick.minor.pad'  : 1,
        #'ytick.minor.pad'  : 1,
        #'text.usetex': True,
        #'figure.figsize': [fig_width,fig_height],
        #'image.cmap': u'rocket',
        #'axes.labelcolor': '.15',
        #'mathtext.fontset':'stixsans',
        #'font.family': 'serif',
        #'font.serif': [u'Computer Modern'],#u'Times New Roman'],
        #'font.family': 'sans-serif',
        #'font.sans-serif': [u'Arial',
        #                    u'DejaVu Sans',
        #                    u'Liberation Sans',
        #                    u'Bitstream Vera Sans',
        #                    u'sans-serif'],
}

# -----------------------------------------------------------------------------
if __name__ == "__main__":
    
    #mpl.style.use('classic') 
    mpl.rcParams.update(params)

    fig = plt.figure()
    fig.set_size_inches(3.375,3.0)

    lft=0.16
    rgt=0.88
    top= 0.97
    bot=0.12
    v2=0.85
    h1=0.35
    dd=0.016

    cbarshort = 0.1
    
    gs1 = gridspec.GridSpec(1,1)
    gs1.update(left=lft, bottom=h1+dd, right=v2-dd, top=top,
               hspace=0.08, wspace=0.3)
    ax1 = plt.subplot(gs1[0,0])

    gs2 = gridspec.GridSpec(1,1)
    gs2.update(left=lft, bottom=bot, right=v2-dd, top=h1-dd,
               hspace=0.08, wspace=0.3)
    ax2 = plt.subplot(gs2[0,0])

    gsc1 = gridspec.GridSpec(1,1)
    gsc1.update(left=v2+dd, bottom=h1+dd+cbarshort,
                right=rgt, top=top-cbarshort,
               hspace=0.08, wspace=0.3)
    cax1 = plt.subplot(gsc1[0,0])


    pc = f2dplot('cohesion_0',ax1)
    
    cbar = fig.colorbar(pc,cax=cax1)
    cax1.annotate(r'shear cohesion $\tau$ [MPa]',
                  xy=(6.5,0.82), xycoords='axes fraction', rotation=90, size=6)


    # plot initial setup
    X = []
    with open(exname+".coord",'r') as fl:
        lines = fl.readlines()
    for line in lines:
        X.append(float(line.strip().split()[0]))
    
    with open(exname+"-DataFiles/cohesion_0.out",'r') as fl:
        line = fl.readline()
    coh0 = [float(i) for i in line.strip().split()]

    with open(exname+"-DataFiles/cohesion_1.out",'r') as fl:
        line = fl.readline()
    coh1 = [float(i) for i in line.strip().split()]

    with open(exname+"-DataFiles/tau_max.out",'r') as fl:
        line = fl.readline()
    taumax = [float(i) for i in line.strip().split()]

    ax2.plot(X, taumax, '-', label=r'$\tau_{max}$')
    ax2.plot(X, coh0, '-', label=r'$\tau_0$')
    ax2.plot(X, coh1, '-', label=r'$\tau_1$')

    ax2.set_xlim([0,1])
    ax2.set_ylim([0,4e6])
    
    # normalize
    ax1.set_yticklabels(['{:.2f}'.format(1e3*i) for i in ax1.get_yticks()])
    cax1.set_yticklabels(['{:.1f}'.format(float(i.get_text())*1e-6)
                          for i in cax1.get_yticklabels()])
    ax2.set_yticklabels(['{:.1f}'.format(1e-6*i) for i in ax2.get_yticks()])
    
    ax2.set_xlabel(r'space, $x$ [m]')
    ax2.set_ylabel(r'stress, $\tau$ [MPa]')
    ax1.set_ylabel(r'time, $t$ [ms]')
    ax1.set_xticklabels([])
    
    ax2.legend(loc='center left', bbox_to_anchor=(1., 0.5))
    
    fig.savefig('fracture_2d_example.png',dpi=300)
    
