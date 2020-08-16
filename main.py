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


dcd='/media/eniac/mdd/paper_membranas/POPG/sim1_popg/MD/sim1_POPG_500ns_origin.dcd'
psf='/media/eniac/mdd/paper_membranas/POPG/sim1_popg/MD/popg_snx1_ions.psf'
#llamo a mi clase tratectoria rep0
Trajectory1=Trajectory(psf,dcd,first=0,last=-1,stride=100,waitfor=-1)
print(Trajectory1.num_frames())


vector_distancia=Trajectory1.distance_center_mass("protein","resname POPG")
time_line=np.linspace(0, Trajectory1.num_frames()/10, num=Trajectory1.num_frames())


#rmds=Trajectory.rmsd_time("protein")
#rmsf=Trajectory.rmsf_time2("name CA")
Trajectory1.close()
plt.plot(time_line, vector_distancia, label='sim 2')


#plt.savefig('distanceMEMr2_200ns.jpg', format='jpg', dpi=900)
#Trajectory2=Trajectory('mem.tox.psf','sim3xpopc_1x_popg_prot_mem.dcd')
#vector_distancia=Trajectory2.distance_center_mass("resid 9 and protein","resname POPC or resname POPG")
#plt.plot(time_line, vector_distancia, label='sim 1')
#Trajectory2.close()
