#!/usr/bin/env python
from __future__ import print_function
from __future__ import division
from glob import glob
import sys
import struct
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import *
import matplotlib.tri as tri
from plotting_utils import *


def inspect_results(full_path, bname):
  domain_xfactor = float(bname.split('_fx')[1].split('_')[0])
  domain_zfactor = float(bname.split('_fz')[1].split('_')[0])

  length_x_rpt = 30e3
  length_y_rpt = 15e3
  length_x = domain_xfactor * length_x_rpt
  length_y = domain_zfactor * length_y_rpt

  stations = [[0.0  ,  6e3],
              [7.5e3, 0.0 ]]
  for station in stations:
    x = station[0]
    y = station[1]
    compare_station(full_path, bname, x, y)
  compare_cplot(full_path,bname)

def compare_station(full_path, bname, x_interest, y_interest):
  spec = bname.split('_')
  domain_xfactor = float(bname.split('_fx')[1].split('_')[0])
  domain_zfactor = float(bname.split('_fz')[1].split('_')[0])
  length_x_rpt = 30e3
  length_y_rpt = 15e3
  length_x = domain_xfactor * length_x_rpt
  length_y = domain_zfactor * length_y_rpt
  nb_nodes_x = int(bname.split('_N')[1].split('_')[0])
  nb_nodes_y = int(nb_nodes_x/length_x*length_y)
  nb_nodes = nb_nodes_x * nb_nodes_y
  dx = length_x / nb_nodes_x
  dy = length_y / nb_nodes_y
  x = np.arange(nb_nodes_x) * dx - length_x / 2
  y = -np.arange(nb_nodes_y) * dy + length_y / 2 
  idx_x = np.argmin(np.abs(x - x_interest))
  idx_y = np.argmin(np.abs(y - y_interest))
  idx = (idx_x - 1) * nb_nodes_y + idx_y
  print('Warning: (%.2e, %.2e) != (%e, %e)' %
        (x_interest, y_interest, x[idx_x], y[idx_y]))

  t = np.fromfile('%s.time' % full_path, sep=' ')
  t = t[1:-1:2]
  nt = len(t)
  delta     = read_data('%s-DataFiles/top_disp_0.out' % full_path, nb_nodes_x, nb_nodes_y, nt, idx_x,idx_y) * 2
  delta_dot = read_data('%s-DataFiles/top_velo_0.out' % full_path, nb_nodes_x, nb_nodes_y, nt, idx_x,idx_y) * 2
  cohesion  = read_data('%s-DataFiles/cohesion_0.out' % full_path, nb_nodes_x, nb_nodes_y, nt, idx_x,idx_y)

  params = {
    'text.latex.preamble': [
        r'\usepackage{siunitx}',
        r'\usepackage{sfmath}',
        r'\sisetup{detect-family = true}',
        r'\usepackage{amsmath}'],
    'font.family': 'sans-serif',
    'font.sans-serif': ['Helvetica']
  }
  plt.rcParams.update(params)
  fig = plt.figure()
  fig.subplots_adjust(top=0.97, bottom=0.1, left=0.1, right=0.97, wspace=0.27, hspace=0.35)
  ax1 = fig.add_subplot(3, 1, 1)
  ax1.plot(t, delta, '-k', label='uguca', linewidth=1, zorder=10)
  ax1.set_xlabel(r'$t$ (sec)')
  ax1.set_ylabel(r'$\delta$ (m)')
  ax2 = fig.add_subplot(3, 1, 2)
  ax2.plot(t, delta_dot, '-k', label='uguca', linewidth=1, zorder=10)
  ax2.set_xlabel(r'$t$ (sec)')
  ax2.set_ylabel(r'$\dot{\delta}$ (m/s)')
  ax3 = fig.add_subplot(3, 1, 3)
  ax3.plot(t, cohesion * 1e-6, '-k', label='uguca', linewidth=1, zorder=10)
  ax3.set_xlabel(r'$t$ (sec)')
  ax3.set_ylabel(r'$\tau$ (MPa)')

  stname='faultst%03ddp%03d'%(x_interest / 100, y_interest / 100)
  if x_interest<0:
    stname='faultst%04ddp%03d'%(x_interest / 100, y_interest / 100)
  data_dir = './ref/'+stname+'/*'
  refs = glob(data_dir)
  for ref in refs:
    author = ref.split('/')[-1].split('.txt')[0]
    if author in dont_exclude_author:
      pass
    else:
      continue
    data = read_scec(ref)
    ax1.plot(data['time'], data['delta_x'], '--', label=author, linewidth=1)
    ax2.plot(data['time'], data['delta_dot_x'], '--', label=author, linewidth=1)
    ax3.plot(data['time'], data['tau_x'], '--', label=author, linewidth=1)

  ax2.legend(loc='upper right', ncol=1, fontsize=7)
  # plt.show()
  plt.tight_layout()
  print('saving %s_%s.pdf' % (bname, stname))
  plt.savefig('%s_%s.pdf' % (bname, stname))



def compare_cplot(full_path,bname):
  spec = bname.split('_')
  nb_nodes_x = int(bname.split('_N')[1].split('_')[0])

  domain_xfactor = float(bname.split('_fx')[1].split('_')[0])
  domain_zfactor = float(bname.split('_fz')[1].split('_')[0])
  
  length_x_rpt = 30e3
  length_y_rpt = 15e3
  length_x = domain_xfactor * length_x_rpt
  length_y = domain_zfactor * length_y_rpt
  nb_nodes_y = int(nb_nodes_x/length_x*length_y)
  nb_nodes = nb_nodes_x * nb_nodes_y
  dx = length_x / nb_nodes_x
  dy = length_y / nb_nodes_y
  x = np.arange(nb_nodes_x) * dx - length_x / 2
  y = np.arange(nb_nodes_y) * dy - length_y / 2
  t = np.fromfile('%s.time' % full_path, sep=' ')
  t = t[1:-1:2]
  nt = len(t)
  xx,yy = np.meshgrid(x,y)
  xx=xx.T
  yy=yy.T
  tau0 = read_data_cplot('%s-DataFiles/cohesion_0.out' % full_path,nb_nodes_x,nb_nodes_y,nt)
  
  rpt_arrival=np.zeros(xx.shape)
  for j in range(nb_nodes_x):
    for k in range(nb_nodes_y):
      i = np.argmax(tau0[:,j,k])
      rpt_arrival[j,k]=t[i]

  fig,ax = plt.subplots(1)
  times=np.arange(0,8,0.5)
  
  stname='cplot'
  data_dir = './ref/'+stname+'/*'
  refs = glob(data_dir)
  i=0
  for ref in refs:
    author = ref.split('/')[-1].split('.txt')[0]
    if author in dont_exclude_author:
      c='bgrcmy'[i]
      i+=1
      pass
    else:
      continue
    X,Y,T = read_scec_cplot(ref)#,nb_nodes_x,nb_nodes_y,nt)
    ax.contour(X,Y,T, colors=c, linestyles='--', linewidth=1,levels=times)
    ax.plot([],[],color=c, linestyle='--',label=author)
  ax.contour(xx,yy,rpt_arrival,colors='k',levels=times)
  ax.plot([],[],color='k', linestyle='-',label='uguca')
  #ax.pcolor(xx.T,yy.T,tau0[200])

  ax.legend(loc='upper right', ncol=1, fontsize=7)

  ax.set_ylim(-length_y_rpt/2,length_y_rpt/2)
  ax.set_xlim(-length_x_rpt/2,length_x_rpt/2)
  ax.set_aspect(1)
  ax.set_xlabel(r'$x$ (m)')
  ax.set_ylabel(r'$y$ (m)')
  plt.tight_layout()
  plt.savefig('%s_%s.pdf' % (bname, stname))

  


if __name__ == '__main__':
  usage = """./TPV3_inspect_results.py <path_to_bname>"""
  if len(sys.argv) != 2:
      sys.exit(usage)
  try:
    full_path = sys.argv[1]
    bname = sys.argv[1].split('/')[-1]
  except:
    sys.exit(usage)
  
  inspect_results(full_path, bname)
