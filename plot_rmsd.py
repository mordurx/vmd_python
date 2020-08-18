
#%% 
from trajectory import Trajectory
import matplotlib.pyplot as plt
import numpy as np
#%%
sim1_dcd= "/mnt/e/rep0_tox_mem500.dcd"
sim1_psf="/mnt/e/rep0.mem.tox.psf"

sim2_dcd= "/mnt/e/rep1_tox_mem500.dcd"
sim2_psf="/mnt/e/rep1.mem.tox.psf"

sim3_dcd= "/mnt/e/rep2_tox_mem500.dcd"
sim3_psf="/mnt/e/rep2.mem.tox.psf"

sim4_dcd= "/mnt/e/rep3_tox_mem500.dcd"
sim4_psf="/mnt/e/rep3.mem.tox.psf"

#%%
traj1=Trajectory(sim1_dcd,sim1_psf,stride=10)
traj2=Trajectory(sim2_dcd,sim2_psf,stride=100)
traj3=Trajectory(sim3_dcd,sim3_psf,stride=10)
traj4=Trajectory(sim4_dcd,sim4_psf,stride=10)

#%%
r1=traj1.rmsd_time("segname TOX")

r2=traj2.rmsd_time("segname TOX")

r3=traj3.rmsd_time("segname TOX")

r4=traj4.rmsd_time("segname TOX")

#contact =traj.porcentaje_contact_fit_list("protein","name P",3000,5000)
time_line=np.linspace(0, traj1.num_frames(), num=traj1.num_frames())

#%%
plt.plot(time_line, r1, label='sim')
plt.plot(time_line, r2, label='rep1')
plt.plot(time_line, r3, label='rep2')
plt.plot(time_line, r4, label='rep3')


plt.show()


# %%
traj1.close()
traj2.close()
traj3.close()
traj4.close()
# %%
