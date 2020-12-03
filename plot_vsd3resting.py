
#%%%
#!/usr/bin/env python3
from trajectory import Trajectory
import matplotlib.pyplot as plt
#Latex Backend for Matplotlib
plt.rc('text', usetex=True)  
plt.rc('font', family='serif') 
import numpy as np

sim1_dcd= "/home/guido/cav23/vsd4_resting_mem_tox_align.dcd"
sim1_psf="/home/guido/cav23/vsd4_resting_mem_tox_align.psf"

sim2_dcd= "/mnt/e/cav2.3/vsd3_resting/eq18_22_850mv_vsd.dcd"
sim2_psf="/mnt/e/cav2.3/vsd3_resting/vsd3_last_frame.Wat.psf"

sim3_dcd= "/mnt/e/paper_membranas/3popc_popg/rep2_tox_mem500.dcd"
sim3_psf="/mnt/e/paper_membranas/3popc_popg/rep2.mem.tox.psf"

sim4_dcd= "/mnt/e/paper_membranas/3popc_popg/rep3_tox_mem500.dcd"
sim4_psf="/mnt/e/paper_membranas/3popc_popg/rep3.mem.tox.psf"

## vsdca1.1 embrionary
sim1_dcd= "/home/mordurx/cav1.1/vs4_cav1.1e_align2.5us.dcd"
sim1_psf="/home/mordurx/cav1.1/vs4_cav1.1e_align2.5us.psf"


traj1=Trajectory(sim1_dcd,sim1_psf,stride=10)
print (traj1.num_frames())

#contact =traj.porcentaje_contact_fit_list("protein","name P",3000,5000)
time_line=np.linspace(0, traj1.num_frames(), num=traj1.num_frames())
#time_line2=np.linspace(traj1.num_frames(), traj1.num_frames()+ traj2.num_frames(), num=traj2.num_frames())

r1=traj1.membrane_center_mass("resname GLU and resid 105 and name CA","resname ARG and resid 115 and name CA")
r2=traj1.membrane_center_mass("resname GLU and resid 105 and name CA","resname ARG and resid 118 and name CA")
r3=traj1.membrane_center_mass("resname GLU and resid 105 and name CA","resname ARG and resid 121 and name CA")
r4=traj1.membrane_center_mass("resname GLU and resid 105 and name CA","resname LYS and resid 124 and name CA")

#r1=traj1.distance_center_mass_Z("resname ARG and resid 115")
#r2=traj1.distance_center_mass_Z("resname ARG and resid 118")
#r3=traj1.distance_center_mass_Z("resname ARG and resid 121")
#r4=traj1.distance_center_mass_Z("resname LYS and resid 124")
traj1.close()
#r2=traj2.distance_center_mass_Z("resname ARG and resid 897")
#traj2.close()

plt.plot(time_line[:1001], r1[:1001], label='ARG 1236 -450mv')
plt.plot(time_line[1001:], r1[1001:], label='ARG 1236 -850mv')

plt.plot(time_line[:1001], r2[:1001], label='ARG 1239 -450mv')
plt.plot(time_line[1001:], r2[1001:], label='ARG 1239 -850mv')

plt.plot(time_line[:1001], r3[:1001], label='ARG 1242 -450mv')
plt.plot(time_line[1001:], r3[1001:], label='ARG 1242 -850mv')

plt.plot(time_line[:1001], r3[:1001], label='ARG 1242 -450mv')
plt.plot(time_line[1001:], r3[1001:], label='ARG 1242 -850mv')

plt.plot(time_line[:1001], r4[:1001], label='LYS 1245 -450mv')
plt.plot(time_line[1001:], r4[1001:], label='LYS 1245 -850mv')
plt.legend(prop={'size': 7})
plt.xlabel(r'\textbf{time (ns)}')
#plt.ylabel(r'\textbf{Z axis (\AA{})}')
plt.ylabel(r'\textbf{distance (\AA{})}')

plt.savefig("/home/mordurx/cav1.1/distance_vsd4_cav1.1.png", dpi=900) 
plt.show()







# %%
