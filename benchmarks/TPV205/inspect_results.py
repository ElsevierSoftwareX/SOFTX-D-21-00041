#!/usr/bin/env python
from __future__ import print_function, division, absolute_import
from glob import glob
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import *
import matplotlib.tri as tri
from plotting_utils import *

dont_exclude_author=np.array(['noda','kaneko','kaneko.2','duan','duan.2'])

dont_exclude_author_c=np.array(['noda','duan','duan.2'])

def inspect_results(full_path, bname):
  stations = glob('./ref/faultst*/')
  if not stations:
    stations = ['faultst000dp030/', 'faultst000dp075/',
                'faultst045dp000/', 'faultst045dp075/',
                'faultst120dp000/', 'faultst120dp075/']
  for station in stations:
    location = station[:-1].split('st')[1].split('dp')
    x = int(location[0]) * 100
    y = int(location[1]) * 100
    compare_station(full_path, bname, x, y)
  compare_cplot(full_path,bname)

def compare_station(full_path, bname, x_interest, y_interest):
  spec = bname.split('_')
  nb_nodes_x = int(bname.split('_N')[1].split('_')[0])
  nb_nodes_y = nb_nodes_x;
  domain_factor = float(bname.split('_s')[1].split('_')[0])
  nb_nodes = nb_nodes_x * nb_nodes_y
  length_x_rpt = 30e3
  length_y_rpt = 30e3
  length_x = domain_factor * length_x_rpt
  length_y = domain_factor * length_y_rpt
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
  plt.tight_layout()
  print('savefig {}_{}.pdf'.format(bname,stname))
  plt.savefig('{}_{}.pdf'.format(bname, stname))



def read_scec(path):
  raw = np.genfromtxt(path, skip_header=1)
  raw = np.delete(raw, 0, 0) # don't know why it didn't discard column headers
  raw = np.array(raw)
  data = dict()
  data['time']        = raw[1::, 0]
  data['delta_x']     = raw[1::, 1]
  data['delta_dot_x'] = raw[1::, 2]
  data['tau_x']       = raw[1::, 3]
  data['delta_y']     = raw[1::, 4]
  data['delta_dot_y'] = raw[1::, 5]
  data['tau_y']       = raw[1::, 6]
  return data


def compare_cplot(full_path,bname):
  spec = bname.split('_')
  nb_nodes_x = int(bname.split('_N')[1].split('_')[0])
  nb_nodes_y = nb_nodes_x;
  domain_factor = float(bname.split('_s')[1].split('_')[0])
  nb_nodes = nb_nodes_x * nb_nodes_y
  length_x_rpt = 30e3
  length_y_rpt = 30e3
  length_x = domain_factor * length_x_rpt
  length_y = domain_factor * length_y_rpt
  dx = length_x / nb_nodes_x
  dy = length_y / nb_nodes_y
  x = np.arange(nb_nodes_x) * dx - length_x / 2
  y = np.arange(nb_nodes_y) * dy - length_y / 2 
  t = np.fromfile('%s.time' % full_path, sep=' ')
  t = t[1:-1:2]
  nt = len(t)
  xx,yy = np.meshgrid(x,y)

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
  ax.contour(xx.T,yy.T,rpt_arrival,colors='k',levels=times)
  ax.plot([],[],color='k', linestyle='-',label='uguca')
  #ax.pcolor(xx.T,yy.T,tau0[200])

  ax.legend(loc='upper right', ncol=1, fontsize=7)

  ax.set_ylim(-length_y_rpt/2,0)
  ax.set_xlim(-length_x_rpt/2,length_x_rpt/2)
  ax.set_aspect(1)
  ax.set_xlabel(r'$x$ (m)')
  ax.set_ylabel(r'$y$ (m)')
  plt.tight_layout()
  print('savefig {}_{}.pdf'.format(bname, stname))
  plt.savefig('{}_{}.pdf'.format(bname, stname))
    

def read_scec_cplot(path):
  raw = np.genfromtxt(path, skip_header=1)
  raw = np.delete(raw, 0, 0) # don't know why it didn't discard column headers
  raw = np.array(raw)
  X = raw[0::, 0]# distance along strike
  Y = raw[0::, 1]# distance along dip
  Z = raw[0::, 2]# rpture time (s)
  ny = len(X[X==X[0]])
  nx = len(X[Y==Y[0]])
  n = len(X)
  if n!=ny*nx:
    raise RuntimeError("nx*ny!=n")
  if True:
    if X[0]==X[1]:
      X= X.reshape(nx,ny)
      Y=-Y.reshape(nx,ny)
      Z= Z.reshape(nx,ny)
    else:
      x=np.linspace(np.min(X),np.max(X),100)
      y=np.linspace(np.min(Y),np.max(Y),100)
      triang = tri.Triangulation(X, Y)
      interpolator = tri.LinearTriInterpolator(triang, Z)
      Xi, Yi = np.meshgrid(x, y)
      Zi = interpolator(Xi, Yi)
      X=Xi
      Y=-Yi
      Z=Zi
    return X,Y,Z


if __name__ == '__main__':
  usage = """./TPV205_inspect_results.py <path_to_bname>"""
  if len(sys.argv) != 2:
      sys.exit(usage)
  try:
    full_path = sys.argv[1]
    bname = sys.argv[1].split('/')[-1]
  except:
    sys.exit(usage)
  
  inspect_results(full_path, bname)
