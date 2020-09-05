
#%%%
from trajectory import Trajectory
import matplotlib.pyplot as plt
import numpy as np

sim1_dcd= "/media/eniac/mdd1/paper_membranas/POPC/sim1_popc/MD/unwrapfixpopc1_700ns.dcd"
sim1_psf="/media/eniac/mdd1/paper_membranas/POPC/sim1_popc/TOX_POPC.psf"

sim2_dcd= "/media/eniac/mdd1/paper_membranas/POPC/sim2_POPC/unwrapped2_700popc.dcd"
sim2_psf="/media/eniac/mdd1/paper_membranas/POPC/sim1_popc/TOX_POPC.psf"

sim3_dcd= "/media/eniac/mdd1/paper_membranas/POPC/sim3_popc/sim3_POPC_unwrap.dcd"
sim3_psf="/media/eniac/mdd1/paper_membranas/POPC/sim3_popc/protein_membrane.psf"

"""
sim4_dcd= "/mnt/e/3popc_popg/rep3_tox_mem500.dcd"
sim4_psf="/mnt/e/3popc_popg/rep3.mem.tox.psf"
"""
#traj1=Trajectory(sim1_dcd,sim1_psf,stride=1)
#traj2=Trajectory(sim2_dcd,sim2_psf,stride=1)
traj3=Trajectory(sim3_dcd,sim3_psf,stride=1)
#traj4=Trajectory(sim4_dcd,sim4_psf,stride=1)
print (traj3.get_molid())
#contact =traj.porcentaje_contact_fit_list("protein","name P",3000,5000)
time_line=np.linspace(0, traj3.num_frames(), num=traj3.num_frames())
#r1=traj1.membrane_center_mass("segname TOX","resname POPG POPC")
#traj1.close()
#r2=traj2.membrane_center_mass("segname TOX","resname POPG POPC")
#traj2.close()
r3=traj3.membrane_center_mass("segname TOX","resname POPG POPC")
traj3.close()
#r4=traj4.membrane_center_mass("segname TOX","resname POPG POPC")
#traj4.close()


#plt.plot(time_line, r1, label='sim1')
#plt.plot(time_line, r2, label='sim2')
plt.plot(time_line, r3, label='sim3')
#plt.plot(time_line, r4, label='sim4')
plt.legend()
#plt.savefig("/mnt/e/3popc_popg/distance_com.png", dpi=900) 
plt.show()







# %%
