from __future__ import print_function, division, absolute_import
import numpy as np
import struct

def read_full_data(path, nx, ny, nt):
  # Only reads binary file (float32)
  file = open(path, 'r')
  data=np.fromfile(file,dtype=np.float32)
  data=data[:nt*nx*ny]
  data=data.reshape((nt,nx,ny))
  file.close()
  return data

def read_data(path, nx, ny, nt, idxx,idxy):
  # Only reads binary file (float32)
  data= read_full_data(path,nx,ny,nt)
  return data[:,idxx,idxy]

def read_data_cplot(path, nx, ny, nt):
  return read_full_data(path,nx,ny,nt)

