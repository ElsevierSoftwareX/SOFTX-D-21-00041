#!/usr/bin/env python

###
# @file   kernel_plotter.py
#
# @author David S. Kammer <dkammer@ethz.ch>
# @author Gabriele Albertini <ga288@cornell.edu>
# @author Chun-Yu Ke <ck659@cornell.edu>
#
# @date creation: Tue Jun 30 2020
# @date last modification: Tue Jun 30 2020
#
# @brief  python script to plot precomputed kernels from files
#
#
# Copyright (C) 2021 ETH Zurich (David S. Kammer)
#
# This file is part of uguca.
#
# uguca is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# uguca is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public License
# along with uguca.  If not, see <https://www.gnu.org/licenses/>.
####

import sys
import numpy as np
import matplotlib.pyplot as plt

from kernel_io import *

argv=sys.argv
if len(argv)==2:
    nu=argv[1]
    pstress=''
elif len(argv)==3:
    nu=argv[1]
    pstress='_pstress'
else:
    raise RuntimeError('usage: ./kernel_plotter.py nu [pstress]')

fig,axes=plt.subplots(2,sharex=True)

fnames=["nu{}{}_{}.txt".format(nu,pstress,k) for k in ['h00','h01','h11']]
fnames.append("h22.txt")

for ftxt in fnames:
    try:
        ft,tt=load_kernel_txt(ftxt)
    except:
        raise RuntimeError("could not load kernel {}".format(ftxt))
    else:
        axes[0].plot(tt,ft,
                     label="$H_{"+ftxt.split('h')[1][:2]+"}(T)$")
        axes[1].plot(tt[:-1],np.diff(ft)/np.diff(tt),
                     label="$H_{"+ftxt.split('h')[1][:2]+"}(T)$")
    
axes[0].set_ylabel('convolution kernels $H(T)$')
axes[1].set_ylabel('$\partial_T H(T)$')
axes[1].set_xlabel('$T$')

plt.tight_layout()    
plt.show()
