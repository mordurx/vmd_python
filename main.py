#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 15:42:13 2019

@author: eniac
"""
from trajectory import Trajectory
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.style
import matplotlib as mpl

'''
paths files
'''
query='index 5799 5800 5801'
path="/media/eniac/mdd1/paper_membranas/analisis/coef_diff_agua/"

dcds=['mytraj.dcd','mytra22j.dcd','tox_unwrap1micro.dcd','unwrapfix_tox_micro.dcd','unwrapped75ns_agua.dcd']
dcd1=path+str(dcds[4])
psfs=['snx_water_ions.psf','snx.rep.cen.psf','solvate.psf']
psf1=path+str(psfs[2])

filename='agua_agua75ns'+query+'.dat'
filepathname=path+str(filename)

#llamo a mi clase tratectoria rep0
Trajectory1=Trajectory(dcd1,psf1,first=100000,last=150000,stride=1,waitfor=-1)
print(Trajectory1.num_frames())
com_Traj_bigdcd=Trajectory1.center_mass(query,filepathname)
time_line=np.linspace(0, Trajectory1.num_frames(), num=Trajectory1.num_frames())
com_Traj_bigdcd=np.transpose(com_Traj_bigdcd)
plt.plot(time_line, com_Traj_bigdcd[0], label='x')
plt.plot(time_line, com_Traj_bigdcd[1], label='y')
plt.plot(time_line, com_Traj_bigdcd[2], label='z')
plt.legend()
plt.show()
#com_Traj_bigdcd=np.transpose(com_Traj_bigdcd)
Trajectory1.close()
# ####### traj vmd 
# Trajectory2=Trajectory(dcd2,psf2,first=0,last=-1,stride=100,waitfor=-1)
# print(Trajectory2.num_frames())
# com_Traj_bigdcd2=Trajectory2.center_mass("segname TOX")
# time_line=np.linspace(0, Trajectory2.num_frames(), num=Trajectory2.num_frames())
# #plt.plot(time_line, vector_distancia, label='sim 2')
# com_Traj_bigdcd2=np.transpose(com_Traj_bigdcd2)
# Trajectory2.close()


# plt.plot(time_line, com_Traj_bigdcd[2], label="x bigdcd")
# plt.plot(time_line, com_Traj_bigdcd2[2], label="x vmd pbc")
# #plt.plot(time_line, com_Traj_bigdcd[2], label="z")
# plt.legend()
# plt.show()
#vector_distancia=Trajectory1.distance_center_mass("protein","resname POPG")
#time_line=np.linspace(0, Trajectory1.num_frames()/10, num=Trajectory1.num_frames())


#rmds=Trajectory.rmsd_time("protein")
#rmsf=Trajectory.rmsf_time2("name CA")
#Trajectory1.close()
#plt.plot(time_line, vector_distancia, label='sim 2')


#plt.savefig('distanceMEMr2_200ns.jpg', format='jpg', dpi=900)
#Trajectory2=Trajectory('mem.tox.psf','sim3xpopc_1x_popg_prot_mem.dcd')
#vector_distancia=Trajectory2.distance_center_mass("resid 9 and protein","resname POPC or resname POPG")
#plt.plot(time_line, vector_distancia, label='sim 1')
#Trajectory2.close()


