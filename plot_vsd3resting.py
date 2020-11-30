
#%%%
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

traj1=Trajectory(sim1_dcd,sim1_psf,stride=10)
print (traj1.num_frames())

#contact =traj.porcentaje_contact_fit_list("protein","name P",3000,5000)
time_line=np.linspace(0, traj1.num_frames(), num=traj1.num_frames())
#time_line2=np.linspace(traj1.num_frames(), traj1.num_frames()+ traj2.num_frames(), num=traj2.num_frames())
r1=traj1.distance_center_mass_Z("resname LYS and resid 99")
r2=traj1.distance_center_mass_Z("resname ARG and resid 102")
r3=traj1.distance_center_mass_Z("resname ARG and resid 105")
traj1.close()
#r2=traj2.distance_center_mass_Z("resname ARG and resid 897")
#traj2.close()

plt.plot(time_line[:1001], r1[:1001], label='LYS 99 -450mv')
plt.plot(time_line[1001:], r1[1001:], label='LYS 99 -850mv')

plt.plot(time_line[:1001], r2[:1001], label='ARG 102 -450mv')
plt.plot(time_line[1001:], r2[1001:], label='ARG 102 -850mv')

plt.plot(time_line[:1001], r3[:1001], label='ARG 105 -450mv')
plt.plot(time_line[1001:], r3[1001:], label='ARG 105 -850mv')
plt.legend(prop={'size': 7})
plt.xlabel(r'\textbf{time (ns)}')
plt.ylabel(r'\textbf{Z axis (\AA{})}')

plt.savefig("/home/guido/cav23/distance_vsd4_arg.png", dpi=900) 
plt.show()







# %%
