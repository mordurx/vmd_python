
#%%%
from trajectory import Trajectory
import matplotlib.pyplot as plt
import numpy as np

sim1_dcd= "/mnt/e/3popc_popg/rep0_tox_mem500.dcd"
sim1_psf="/mnt/e/3popc_popg/rep0.mem.tox.psf"

sim2_dcd= "/mnt/e/3popc_popg/rep1_tox_mem500.dcd"
sim2_psf="/mnt/e/3popc_popg/rep1.mem.tox.psf"

sim3_dcd= "/mnt/e/3popc_popg/rep2_tox_mem500.dcd"
sim3_psf="/mnt/e/3popc_popg/rep2.mem.tox.psf"

sim4_dcd= "/mnt/e/3popc_popg/rep3_tox_mem500.dcd"
sim4_psf="/mnt/e/3popc_popg/rep3.mem.tox.psf"

traj1=Trajectory(sim1_dcd,sim1_psf,stride=100)
traj2=Trajectory(sim2_dcd,sim2_psf,stride=100)
traj3=Trajectory(sim3_dcd,sim3_psf,stride=100)
traj4=Trajectory(sim4_dcd,sim4_psf,stride=100)
print (traj1.get_molid())
#contact =traj.porcentaje_contact_fit_list("protein","name P",3000,5000)
time_line=np.linspace(0, traj1.num_frames(), num=traj1.num_frames())
r1=traj1.rmsd_time("segname TOX")
traj1.close()
r2=traj2.rmsd_time("segname TOX")
traj2.close()
r3=traj3.rmsd_time("segname TOX")

traj3.close()
r4=traj4.rmsd_time("segname TOX")

traj4.close()

plt.plot(time_line, r1, label='sim1')
plt.plot(time_line, r2, label='sim2')
plt.plot(time_line, r3, label='sim3')
plt.plot(time_line, r4, label='sim4')
plt.legend()
plt.savefig("/mnt/e/3popc_popg/rmds_3POPC_POPG.png", dpi=900) 
plt.show()







# %%
