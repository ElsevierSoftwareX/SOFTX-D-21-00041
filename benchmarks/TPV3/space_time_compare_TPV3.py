#!/usr/bin/env python
from __future__ import print_function
import matplotlib.pyplot as plt
import numpy as np
import sys
from mpl_toolkits.mplot3d import axis3d
from matplotlib.colors import LogNorm

#from postprocess import get_input_data

from ifasha.datamanager import FieldId
from ifasha.datamanager import FieldCollectionAnalysis
from ifasha.datamanager import DataManagerAnalysis

#from simid_bname_data import simid_to_bname

def space_time(bname, group, fldid, **kwargs):

    wdir = kwargs.get('wdir','data/')

    # number of points along space
    nb_x_elements = kwargs.get('nb_x_elements',400)
    start_fct = kwargs.get('start_fct',0.)
    end_fct = kwargs.get('end_fct',1.)

    start_time = kwargs.get('start_time',0)
    end_time = kwargs.get('end_time',None) # None if until end of sim
    nb_t_points = kwargs.get('nb_t_points',400)

    zmax = kwargs.get('zmax',None)
    zmin = kwargs.get('zmin',None)

    # in cases bnames are actually sim-ids
    # bname = simid_to_bname(bname,**kwargs)

    # get input data of simulation
    #input_data = get_input_data(bname,group,**kwargs)
    
    # add plot on ax of figure if already provided
    ax = kwargs.get('ax', None)
    new_figure = True if ax is None else False
    if kwargs.get('no_fig'): new_figure = False 
    if new_figure:
        fig = plt.figure()
        ax = fig.add_subplot(111)  
    
    # ---------------------------------------------------------------------
    # load the simulation data
    dma = DataManagerAnalysis(bname,wdir)
    data = dma(group)


    z_coord_idx = kwargs.get('z_coord_idx',2)

    try:
        data.get_full_field(FieldId('coord',2))
    except:
        is_3d_in_2d = False
    else:
        is_3d_in_2d = True
        z_coord_tmp = kwargs.get('z_coord',0.0)
        z_coords = np.array(sorted(set(data.get_field_at_t_index(FieldId('coord',z_coord_idx),0)[0])))
        z_coord = z_coords[np.argmin(np.abs(z_coords-z_coord_tmp))]
        print('is 3d in 2d, z_coord_tmp = ',z_coord_tmp,'z_coord = ',z_coord)

    # ---------------------------------------------------------------------
    start = data.get_index_of_closest_time(FieldId('time'),start_time)
    if end_time is None:
        last = data.get_t_index('last')
    else:
        last = data.get_index_of_closest_time(FieldId('time'),end_time)

    delta_t = max(1,(last - start) / nb_t_points)
    print('start',start,'last',last,'delta_t',delta_t)

    # position
    x_coord_idx = 0
    if z_coord_idx==0:
        x_coord_idx=2;

    print('x_coord_idx', x_coord_idx)

    pos_fldid = FieldId('coord',x_coord_idx)

    
    if is_3d_in_2d:
        print('z_coord',z_coord)
        idx = data.get_index_of_closest_position(FieldId('coord',z_coord_idx),z_coord)
        zpos = data.get_field_at_node_index(FieldId('coord',z_coord_idx),idx)[0][0]
        idcs = data.get_indices_of_nodes_on_line(FieldId('coord',z_coord_idx),zpos)
        nb_elements = len(idcs)
    else:
        nb_elements = len(data.get_field_at_t_index(pos_fldid,0)[0])
        
    st_x_elements = int(start_fct * nb_elements)
    e_x_elements = int(end_fct * nb_elements)
    int_x_elements = max(1,int((e_x_elements - st_x_elements) / nb_x_elements))
    print('xs',st_x_elements,'xe',e_x_elements,'delta_x',int_x_elements)

    if is_3d_in_2d:
        X,T,Z = data.get_sliced_x_sliced_t_plot_at_node_index(
            pos_fldid,
            np.arange(st_x_elements,e_x_elements,int_x_elements),
            FieldId("time"),
            np.arange(start,last,delta_t),
            fldid,
            idcs)
    else:
        X,T,Z = data.get_sliced_x_sliced_t_plot(
            pos_fldid,
            np.arange(st_x_elements,e_x_elements,int_x_elements),
            FieldId("time"),
            np.arange(start,last,delta_t),
            fldid)

    print(X.shape,T.shape,Z.shape)
    
    print('maxX',np.max(X))
    print('minX',np.min(X))
    print('maxT',np.max(T))
    print('minT',np.min(T))
    print('maxZ',np.max(Z))
    print('minZ',np.min(Z))

    if kwargs.get('no_fig'):
        return X,T,Z

    if zmin is not None and zmax is not None:
        fg1 = ax.pcolor(X,T,Z,vmin=zmin,vmax=zmax)
    elif zmin is not None:
        fg1 = ax.pcolor(X,T,Z,vmin=zmin)
    elif zmax is not None:
        fg1 = ax.pcolor(X,T,Z,vmax=zmax)
    else:
        fg1 = ax.pcolor(X,T,Z)

    #Z = np.gradient(Z)
    #Z = np.abs(Z)
    #Z = Z[1]
    #fg1 = ax.pcolor(X,T,Z, norm=LogNorm(vmin=1, vmax=Z.max()))
    if new_figure:
        ax.set_title(bname+" "+fldid.get_string())

        plt.colorbar(fg1)

    return fg1


plt.rcParams['image.cmap'] = 'jet'#'gist_ncar'#'jet''nipy_spectral'#

bname = sys.argv[1]

fct = float(bname.split('s')[1].split('_')[0])

print(fct)

fctx=fct
fcty=fct*1.5


normal_load = -120e6;
  
f_c = 0.677;
f_r = 0.525;
dc = 0.4;

tau_c = -f_c * normal_load;
tau_r = -f_r * normal_load;

space_time(bname,'interface',FieldId('cohesion',0),z_coord_idx=2,
           z_coord=7500*fcty,
           zmin=tau_r,zmax=tau_c)

plt.xlim(15000*fctx,15000*fctx+15000)
plt.ylim(0,6)
plt.xticks(np.arange(15000*fctx,15000*fctx+15000+1,5000),np.arange(0,15000+1,5000))
plt.xlabel('Fault along inplane (m)')
plt.ylabel('Time (sec)')
plt.tight_layout()

plt.savefig(bname+"_space_time_x.png",dpi=300)

space_time(bname,'interface',FieldId('cohesion',0),z_coord_idx=0,
           z_coord=15000*fctx,
           zmin=tau_r,zmax=tau_c)

plt.xlim(7500*fcty,7500*fcty+7500)
plt.ylim(0,6)
plt.xticks(np.arange(7500*fcty,7500*fcty+7500+1,1000),np.arange(0,7500+1,1000))
plt.xlabel('Fault along inplane (m)')
plt.ylabel('Time (sec)')
plt.tight_layout()

plt.savefig(bname+"_space_time_y.png",dpi=300)



