#!/usr/bin/env python

###
# @file   laplace_inverse.py
#
# @author David S. Kammer <dkammer@ethz.ch>
# @author Gabriele Albertini <ga288@cornell.edu>
# @author Chun-Yu Ke <ck659@cornell.edu>
#
# @date creation: Tue Jun 30 2020
# @date last modification: Tue Jun 30 2020
#
# @brief  python script to generate precomputed kernel files
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

from mpmath  import *
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool

from kernel_io import *

#--------------------------------
# wave speed ratio Cd/Cs

def CdCs(nu,pstress):
    if pstress==False:
        return sqrt(2.0*(1.0-nu)/(1.0-2.0*nu));
    else:
        return sqrt(2.0/(1.0-nu));
#--------------------------------
# define kernels
#--------------------------------
# h00 kernel
# equiv. to Breitenfeld and Geubelle 1997 eq(7) H11
  
def h00(s,cdcs):
    return ((cdcs*cdcs + cdcs*sqrt(s*s+cdcs*cdcs)*sqrt(s*s+1.0) +
	     cdcs*s*sqrt(s*s+cdcs*cdcs) - s*sqrt(1.0+s*s) ) /
            ((sqrt(1.0+s*s)+s) * (s*s+1.0+cdcs*cdcs) ))

#--------------------------------
# h01 kernel
# equiv. to Breitenfeld and Geubelle 1997 eq(7) H12

def h01(s,cdcs):
  return cdcs*(s*s/(cdcs-sqrt(s*s+cdcs*cdcs)*sqrt(1.0+s*s))+1.0)

#--------------------------------
# h11 kernel
# equiv. to Breitenfeld and Geubelle 1997 eq(7) H22

def h11(s,cdcs):
    return (cdcs*s*(sqrt(1.0+s*s)*(-cdcs*cdcs)/(s+sqrt(s*s+cdcs*cdcs))+cdcs)
            /(sqrt(1.0+s*s)*sqrt(s*s+cdcs*cdcs)-cdcs))

#--------------------------------
# h22 kernel
# equiv to eq(8) in Geubelle and Breitenfeld (1997) H33 
# H33(t) = J1(t)/t
# where J1 is the Bessel function

def h22inv(t):
    return besselj(1,t)/t;

# inversion functions for H00, H01 and H11
def h00inv(cdcs,tt,dt,do_save,nu,psss):
    print('invert H00...')
    sol=[invertlaplace(lambda s: h00(s,cdcs),t) for t in tt]
    if do_save:
        save_to_txt(sol,dt,"nu{}{}_h00.txt".format(nu,psss))
    return sol
        
def h01inv(cdcs,tt,dt,do_save,nu,psss):
    print('invert H01...')
    sol=[invertlaplace(lambda s: h01(s,cdcs),t) for t in tt]
    if do_save:
        save_to_txt(sol,dt,"nu{}{}_h01.txt".format(nu,psss))
    return sol
        
def h11inv(cdcs,tt,dt,do_save,nu,psss):
    print('invert H11...')
    sol=[invertlaplace(lambda s: h11(s,cdcs),t) for t in tt]
    if do_save:
        save_to_txt(sol,dt,"nu{}{}_h11.txt".format(nu,psss))
    return sol

    
#--------------------------------
# do the inversion and generate files of H00, H01, H11
def invert(nu=0.25,tcut=100,dt=0.05,plane_stress=False,do_save=True,do_plot=True):
    if nu<0.0 or nu>0.5:
        return RuntimeError("bounds of nu 0.0,0.5")
    if tcut>300:
        return RuntimeError("tcut <= 300")

    # print info for user
    print('nu      = {}'.format(nu))
    print('pstress = {}'.format(plane_stress))
    print('dt      = {}'.format(dt))
    print('tcut    = {}'.format(tcut))

    # for saving file
    if plane_stress:
        psss="_pstress"
    else: psss=""
    
    cdcs=CdCs(nu,plane_stress)

    tt=mp.arange(0,tcut+dt,dt)
    tt[0]+=1e-12

    pool = Pool()
    h00res = pool.apply_async(h00inv, [cdcs,tt,dt,do_save,nu,psss])
    h01res = pool.apply_async(h01inv, [cdcs,tt,dt,do_save,nu,psss])
    h11res = pool.apply_async(h11inv, [cdcs,tt,dt,do_save,nu,psss])

    h00inv_sol = h00res.get()
    h01inv_sol = h01res.get()
    h11inv_sol = h11res.get()
    
    if do_plot:
        plt.plot(tt,h00inv_sol,'.-',label="H00")
        plt.plot(tt,h01inv_sol,'.-',label="H01")
        plt.plot(tt,h11inv_sol,'.-',label="H11")
        plt.xlabel("$t$"); plt.ylabel("$H(t)$")
        plt.legend(loc='best',frameon=False)
        plt.savefig("laplace_inverse.pdf")

#--------------------------------
# generate file for H22 from analytical solution
# note: independent on nu and plane_stress
def analytical_inverse(tcut=100,dt=0.05,do_save=True,do_plot=True):
    print('dt      = {}'.format(dt))
    print('tcut    = {}'.format(tcut))
    tt=mp.arange(0,tcut+dt,dt)
    tt[0]+=1e-12
    print('invert H22...')
    h22inv_anl=[h22inv(t) for t in tt]
    if do_save:
        save_to_txt(h22inv_anl,dt,"h22.txt")
    if do_plot:
        plt.plot(tt,h22inv_anl,'.-',label="H22")
        plt.xlabel("$t$"); plt.ylabel("$H(t)$")
        plt.legend(loc='best',frameon=False)
        plt.savefig("analytical_laplace_inverse.pdf")
        
if __name__=="__main__":

    import sys
    mp.dps = 100; mp.pretty = True
    
    if len(sys.argv)<2 or len(sys.argv)>4:
        raise RuntimeError(
            "\nUsage: ./laplace_inverse.py <nu> [pstress/pstrain] [tcut]\n       ./laplace_inverse.py h22 [tcut]")
    
    elif sys.argv[1]=="h22":
        tcut=100
        if len(sys.argv) == 3:
            tcut=int(sys.argv[2])
            mp.dps=tcut
        analytical_inverse(tcut=tcut)
        
    else:
        nu=float(sys.argv[1])

        pstress = False
        if len(sys.argv) > 2:
            if str(sys.argv[2]) == 'pstress':
                pstress = True
            elif str(sys.argv[2]) == 'pstrain':
                pstress = False
            else:
                raise RuntimeError("did not generate kernel. did you mean 'pstress'?")

        tcut=100
        if len(sys.argv) > 3:
            tcut=int(sys.argv[3])
            mp.dps=tcut
            print('mp.dps = {}'.format(mp.dps))
            
        invert(nu=nu,tcut=tcut,plane_stress=pstress)
