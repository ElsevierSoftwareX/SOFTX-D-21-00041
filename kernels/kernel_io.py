###
# @file   kernel_io.py
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
import numpy as np

# read kernel files
def load_kernel_txt(fname):
    with open(fname) as f:
        lines = f.readlines()

    dt = float(lines[0])
    ftxt=lines[1].split(',')[:-1]
    f = [float(i) for i in ftxt]
    t = np.arange(len(f))*dt

    return f,t

# save kernel files
def save_to_txt(hinv,dt,fname):
    with open(fname, 'w') as f:
        f.write("{:.12e}\n".format(float(dt)))
        for h in hinv:
            f.write("{:.12e},".format(float(h)))
