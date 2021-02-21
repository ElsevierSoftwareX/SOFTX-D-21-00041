#!/usr/bin/env python3
from __future__ import print_function
from __future__ import division
from glob import glob
import sys
import struct
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from plotting_utils import *

dont_exclude_author = np.array(['barall', 'dalguer2', 'dunham', 'kaneko', 'liu'])

dont_exclude_author_c=np.array(['barall', 'dalguer2', 'dunham', 'kaneko'])


def inspect_results(full_path, bname):
  params = {
    'text.latex.preamble': r'\usepackage{siunitx} \sisetup{detect-all} \usepackage{sfmath} \sisetup{detect-family=true} \newcommand{\minus}{\scalebox{0.75}[1.0]{$-$}} \usepackage{amsmath}',
    'font.family': 'sans-serif',
    'font.sans-serif': ['Helvetica'],
    'font.size': 7,
    'axes.labelsize': 7,
    'axes.labelpad': 1,
    'axes.titlesize': 7,
    'axes.titlepad': 3,
    'axes.linewidth': 0.5,
    'lines.markeredgewidth': 0.0,
    'lines.markersize': 5,
    'legend.fontsize': 6,
    'legend.fancybox': False,
    'legend.framealpha': 0,
    'legend.numpoints': 1,
    'legend.borderpad': 0.1,
    'legend.labelspacing': -0.0,
    'legend.columnspacing': 0.0,
    'legend.borderaxespad': 0.1,
    'legend.markerscale': 1.0,
    'legend.handlelength': 1.4,
    'legend.handleheight': 1,
    'legend.handletextpad': 0.5,
    'xtick.labelsize': 7,
    'ytick.labelsize': 7,
    'xtick.major.size': 2,
    'ytick.major.size': 2,
    'xtick.minor.size': 2,
    'ytick.minor.size': 2,
    'xtick.major.width': 0.3,
    'ytick.major.width': 0.3,
    'xtick.minor.width': 0.3,
    'ytick.minor.width': 0.3,
    'xtick.major.pad': 1.1,
    'ytick.major.pad': 1.1,
    'xtick.minor.pad': 1.1,
    'ytick.minor.pad': 1.1,
    'lines.linewidth': 1,
    'text.usetex': True,
    # 'text.latex.preview': True,
    'figure.figsize': [3.5, 3.5]  # [3.74, 4.53]
    }  # 7.2 for full-width, 3.5 for half-width
  plt.rcParams.update(params)
  fig1 = plt.figure()
  fig1.subplots_adjust(top=0.957, bottom=0.074, left=0.121,
                     right=0.97, wspace=0.400, hspace=0.06)
  gs = gridspec.GridSpec(4, 3, height_ratios=[2, 0.5, 1, 1], width_ratios=[1, 1, 1])
  ax0 = fig1.add_subplot(gs[0, :])
  ax1a = fig1.add_subplot(gs[2, 0])
  ax1b = fig1.add_subplot(gs[2, 1])
  ax1c = fig1.add_subplot(gs[2, 2])
  ax2a = fig1.add_subplot(gs[3, 0])
  ax2b = fig1.add_subplot(gs[3, 1])
  ax2c = fig1.add_subplot(gs[3, 2])

  compare_cplot(full_path, bname, ax0, fig1)
  
  stations = ['faultst000dp030/', 'faultst090dp075/', 'faultst120dp030/']
  labels = ['A', 'B', 'C']
  ax1s = [ax1a, ax1b, ax1c]
  ax2s = [ax2a, ax2b, ax2c]
  legend_ons = [True, False, False]
  for station, lbl, ax1, ax2, legend_on in zip(stations, labels, ax1s, ax2s, legend_ons):
    location = station.split('faultst')[1].split('/')[-2]
    location = location.split('dp')
    x = int(location[0]) * 100
    y = int(location[1]) * 100
    compare_station(full_path, bname, x, y, ax1, ax2, legend_on)
    rotate_yticklabels(ax1)
    rotate_yticklabels(ax2)
    add_station(ax0, x, y, lbl)
    ax1.set_title(r'$\mathbf{%s}$' % lbl)
  
  x_offset = -0.070
  y_offset = -0.000
  axes = [ax0, ax1a, ax1b, ax1c, ax2a, ax2b, ax2c]
  for i in range(len(axes)):
      pa = axes[i].get_position()
      fig1.text(pa.xmin + x_offset, pa.ymax + y_offset,
               r'\textbf{%s}' % chr(97 + i), fontsize=8, ha='center', va='center')

  # plt.show()
  fig1.savefig('TPV101_result_summary.png', dpi=600)


def rotate_yticklabels(ax, minor_visible=False):
    for label in ax.get_yticklabels():
        label.set_ha('right')
        label.set_va('center')
        label.set_rotation(90)
    ax.tick_params(which='major', length=4)
    for label in ax.yaxis.get_minorticklabels():
        label.set_visible(minor_visible)
        if minor_visible:
            label.set_ha('right')
            label.set_va('center')
            label.set_rotation(90)

def set_visible_xticklabels(ax, visible):
    for tl in ax.get_xticklabels():
        tl.set_visible(visible)

def add_station(ax, x, y, lbl):
  x = x * 1e-3
  y = y * 1e-3
  ax.plot(x, y, 'k.')
  ax.annotate(lbl, [x + 0.2, y - 0.2], color='k')

def compare_cplot(full_path, bname, ax, fig):
  spec = bname.split('_')
  nb_nodes_x = int(bname.split('_Nx')[1].split('_')[0])
  nb_nodes_z = int(bname.split('_Nz')[1].split('_')[0])
  domain_factor = float(bname.split('_s')[1].split('_')[0])
  nb_nodes = nb_nodes_x * nb_nodes_z
  length_x_rpt = 36e3
  length_z_rpt = 18e3
  length_x = domain_factor * length_x_rpt
  length_z = domain_factor * length_z_rpt
  dx = length_x / nb_nodes_x
  dz = length_z / nb_nodes_z
  x = np.arange(nb_nodes_x) * dx - length_x / 2
  z = -7.5e3 + length_z / 2 - np.arange(nb_nodes_z) * dz
  t = np.fromfile('%s.time' % full_path, sep=' ')
  t = t[1:-1:2]
  nt = len(t)
  xx, zz = np.meshgrid(x, z, indexing='ij')
  tau0 = read_data_cplot('%s-DataFiles/cohesion_0.out' % full_path, nb_nodes_x, nb_nodes_z, nt)
  rpt_arrival = np.zeros(xx.shape)
  for j in range(nb_nodes_x):
    for k in range(nb_nodes_z):
      i = np.argmax(tau0[:, j, k])
      rpt_arrival[j, k] = t[i]
  ax.set_title(r'$\mathbf{rupture\ front}$')

  times = np.arange(0, 8, 0.5)
  h = ax.contourf(xx * 1e-3, -zz*1e-3, rpt_arrival, levels=times, cmap=plt.cm.viridis, extend='max')
  h2 = ax.contour(xx * 1e-3, -zz*1e-3, rpt_arrival, levels=times, cmap=plt.cm.gray)
  ax.set_ylim(-1, 16)
  ax.set_xlim(-16, 16)
  ax.set_aspect(1)
  ax.set_xlabel(r'$x$ (km)')
  ax.set_ylabel(r'$y$ (km)')
  ax.invert_yaxis()
  divider = make_axes_locatable(ax)
  spacer = divider.append_axes('right', size='10%', pad=0.05)
  spacer.set_visible(False)
  cax = inset_axes(ax, width="5%", height="100%", loc='lower left',
                     bbox_to_anchor=(1.01, 0., 1, 1),
                     bbox_transform=ax.transAxes,
                     borderpad=0)
  cb = fig.colorbar(h, cax=cax)
  cb.set_label(r'time (sec)', rotation=270, labelpad=8)
  cb.add_lines(h2)
  for level in h2.collections:
    for kp, path in reversed(list(enumerate(level.get_paths()))):
        verts = path.vertices  # (N,2)-shape array of contour line coordinates
        diameter = np.max(verts.max(axis=0) - verts.min(axis=0))
        if diameter < 1:  # threshold to be refined for your actual dimensions!
            del(level.get_paths()[kp])  # no remove() for Path objects:(

def compare_station(full_path, bname, x_interest, z_interest, ax1, ax2, legend_on=False):
  spec = bname.split('_')
  nb_nodes_x = int(spec[1][2::])
  nb_nodes_z = int(spec[2][2::])
  domain_factor = float(spec[3][1::])
  nb_nodes = nb_nodes_x * nb_nodes_z
  length_x_rpt = 36e3
  length_z_rpt = 18e3
  length_x = domain_factor * length_x_rpt
  length_z = domain_factor * length_z_rpt
  dx = length_x / nb_nodes_x
  dz = length_z / nb_nodes_z
  x = np.arange(nb_nodes_x) * dx - length_x / 2
  z = 7.5e3 + length_z / 2 - np.arange(nb_nodes_z) * dz
  idx_x = np.argmin(np.abs(x - x_interest))
  idx_z = np.argmin(np.abs(z - z_interest))
  idx = (idx_x - 1) * nb_nodes_z + idx_z
  print('Warning: (%.2e, %.2e) != (%e, %e)' %
        (x_interest, z_interest, x[idx_x], z[idx_z]))

  t = np.fromfile('%s.time' % full_path, sep=' ')
  t = t[1:-1:2]
  nt = len(t)
  delta = read_data('%s-DataFiles/top_disp_0.out' % full_path, nb_nodes_x, nb_nodes_z, nt, idx_x, idx_z) * 2
  delta_dot = read_data('%s-DataFiles/top_velo_0.out' % full_path, nb_nodes_x, nb_nodes_z, nt, idx_x, idx_z) * 2
  cohesion = read_data('%s-DataFiles/cohesion_0.out' % full_path, nb_nodes_x, nb_nodes_z, nt, idx_x, idx_z)

  ax1.plot(t, cohesion * 1e-6, '-k', label='uguca', linewidth=1, zorder=10)
  ax1.set_ylabel(r'$\tau$ (MPa)')
  # ax2.plot(t, delta, '-k', label='uguca', linewidth=1, zorder=10)
  # ax2.set_ylabel(r'$\Delta u$ (m)')
  ax2.plot(t, delta_dot, '-k', label='uguca', linewidth=1, zorder=10)
  ax2.set_ylabel(r'$\Delta \dot{u}$ (m/s)')
  ax2.set_xlabel(r'time (s)')
  set_visible_xticklabels(ax1, False)
  xticks = range(0, round(t[-1]+1), 5)
  ax1.set_xticks(xticks)
  ax2.set_xticks(xticks)

  data_dir = './ref/faultst%03ddp%03d/*' % (x_interest / 100, z_interest / 100)
  refs = glob(data_dir)
  for ref in refs:
    author = ref.split('/')[-1].split('.')[0]
    data = read_scec(ref)
    ax1.plot(data['time'], data['tau_x'], '--', label=author, linewidth=1)
    # ax2.plot(data['time'], data['delta_x'], '--', label=author, linewidth=1)
    ax2.plot(data['time'], data['delta_dot_x'], '--', label=author, linewidth=1)
  if legend_on:
    ax1.legend(loc='upper right', ncol=1, fontsize=6)


def read_scec(path):
  raw = np.genfromtxt(path, skip_header=1)
  raw = np.delete(raw, 0, 0)  # don't know why it didn't discard column headers
  raw = np.array(raw)
  data = dict()
  data['time'] = raw[1::, 0]
  data['delta_x'] = raw[1::, 1]
  data['delta_dot_x'] = raw[1::, 2]
  data['tau_x'] = raw[1::, 3]
  data['delta_y'] = raw[1::, 4]
  data['delta_dot_y'] = raw[1::, 5]
  data['tau_y'] = raw[1::, 6]
  data['sigma'] = raw[1::, 7]
  data['log_theta'] = raw[1::, 8]
  return data


def read_scec_cplot(path):
  raw = np.genfromtxt(path, skip_header=1)
  raw = np.delete(raw, 0, 0)  # don't know why it didn't discard column headers
  raw = np.array(raw)
  X = raw[0::, 0]  # distance along strike
  Y = raw[0::, 1]  # distance along dip
  Z = raw[0::, 2]  # rpture time (s)
  ny = len(X[X == X[0]])
  nx = len(X[Y == Y[0]])
  n = len(X)
  if n != ny*nx:
    raise RuntimeError("nx*ny!=n")
  if True:
    if X[0] == X[1]:
      X = X.reshape(nx, ny)
      Y = -Y.reshape(nx, ny)
      Z = Z.reshape(nx, ny)
    else:
      x = np.linspace(np.min(X), np.max(X), 100)
      y = np.linspace(np.min(Y), np.max(Y), 100)
      triang = tri.Triangulation(X, Y)
      interpolator = tri.LinearTriInterpolator(triang, Z)
      Xi, Yi = np.meshgrid(x, y)
      Zi = interpolator(Xi, Yi)
      X = Xi
      Y = -Yi
      Z = Zi
    return X, Y, Z

if __name__ == '__main__':
  usage = """./inspect_results.py <path_to_bname>"""
  if len(sys.argv) != 2:
      sys.exit(usage)
  try:
    full_path = sys.argv[1]
    bname = sys.argv[1].split('/')[-1]
  except:
    sys.exit(usage)

  inspect_results(full_path, bname)
