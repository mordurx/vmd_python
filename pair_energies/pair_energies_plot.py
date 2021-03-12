#!/usr/bin/env python3
import re, sys, os
#from matplotlib import rc
#sys.path.insert(0, '/home/mordurx/vmd_python')
#from trajectory import Trajectory
#import numpy as np
import matplotlib.pyplot 
#import matplotlib.pyplot as plt
from pair_energies_method import *

energies_agua=pair_interations_energies('pair_energies/sim2_popc_pair_elec.out', 'ENERGY:')
#energies_membrane=pair_interations_energies('pair_energies/pair_int_membrane_snx.out', 'ENERGY:')

#index=np.linspace(0,20001, num=20001)
#plt.plot(index,energies_agua['VDW'],label='VDW0') 
#plt.savefig('pair_energies/membrane_popc_popg_vdw.png', dpi=900) 

# sim2_dcd= "/home/mordurx/resume_paper_membranas/sim_POPC_resume/sim2_unwrap_700nspopc.dcd"
# sim2_psf="/home/mordurx/resume_paper_membranas/sim_POPC_resume/sim2_unwrap_700nspopc.psf"

# traj2=Trajectory(sim2_dcd,sim2_psf,stride=1,last=7000)
# #time_line=np.linspace(0,traj1.num_frames() , num=traj1.num_frames())
# r2=traj2.distance_center_mass("protein","resname POPG POPC")
# traj2.close()
# label=['sim2_popc_VDW','sim2_popc_ELECT']
# plt.rcParams.update({'figure.figsize':(10,8), 'figure.dpi':100})
# #plt.scatter(x, y, c=energies_pair0['VDW'], cmap='Spectral')
# plt.scatter(r2,energies_pair0['ELECT'], label=label[1],c=energies_pair0['ELECT'], cmap='Spectral', s=4)
# plt.legend(prop={'size': 7})
#plt.colorbar()

# plt.ylabel('energy POPC-toxin (kcal/mol)')
# rc('text', usetex=True)
# plt.xlabel(r'COM distance')

# plt.colorbar(orientation='horizontal')
# plt.savefig('pair_energies/popc_vdw_distance.png', dpi=900) 
# %%
