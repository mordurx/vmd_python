
#%%%
from trajectory import Trajectory
import matplotlib.pyplot as plt
import numpy as np
"""
3popg-popc
sim1_dcd= "/media/eniac/mdd1/paper_membranas/sim1_3popg_1popc/sim1_3popg_popc_fixed_unwrapped.dcd"
sim1_psf="/media/eniac/mdd1/paper_membranas/sim1_3popg_1popc/agua_tox_mem_Sim1_ions.Wat.psf"

sim2_dcd= "/media/eniac/mdd1/paper_membranas/sim2_3popg_popc/sim2_3popg_unwrap.dcd"
sim2_psf="/media/eniac/mdd1/paper_membranas/sim2_3popg_popc/3popg_popc_tox.psf"

sim3_dcd= "/media/eniac/mdd1/paper_membranas/sim3_3popg_popc/MD/sim3_3popg_popc.dcd"
sim3_psf="/media/eniac/mdd1/paper_membranas/sim3_3popg_popc/agua_tox_mem_Sim1_ions.Wat.psf"
sim4_dcd= "/media/eniac/mdd1/paper_membranas/sim4_3popg_popc/sim4_3popg_popc.dcd"
sim4_psf="/media/eniac/mdd1/paper_membranas/sim4_3popg_popc/agua_tox_mem_Sim1_ions.Wat.psf"
"""

sim1_dcd='/media/eniac/mdd1/paper_membranas/xi_jinghao_toxin/unwrapped_popg_XI.dcd'
sim1_psf='/media/eniac/mdd1/paper_membranas/xi_jinghao_toxin/sim1_IX_popg_ions.Wat.psf'

#sim2_dcd='/media/eniac/mdd1/paper_membranas/POPG/unwrap_sim2_popg.dcd'
#sim2_psf='/media/eniac/mdd1/paper_membranas/POPG/sim2_popg_ions.Wat.psf'
sim2_dcd='/media/eniac/mdd1/paper_membranas/xi_jinghao_toxin/unwrapped_popc_popg_XI.dcd'
sim2_psf='/media/eniac/mdd1/paper_membranas/xi_jinghao_toxin/sim1_ix_1_popc_popg.Wat.psf'


sim3_dcd= "/media/eniac/mdd1/paper_membranas/xi_jinghao_toxin/unwrapped_3popg_popg_XI.dcd"
sim3_psf="/media/eniac/mdd1/paper_membranas/xi_jinghao_toxin/ix_3popg.psf"

sim4_dcd='/media/eniac/mdd1/paper_membranas/xi_jinghao_toxin/simXI_3popc_popg.dcd'
sim4_psf='/media/eniac/mdd1/paper_membranas/xi_jinghao_toxin/IX_3popc.psf'


sim5_dcd='/media/eniac/mdd1/paper_membranas/xi_jinghao_toxin/unwrap_XI_sim3_popc.dcd'
sim5_psf='/media/eniac/mdd1/paper_membranas/xi_jinghao_toxin/sim3_ix_popc_ions.Wat.psf'




traj1=Trajectory(sim1_dcd,sim1_psf,stride=1)
time_line=np.linspace(0,traj1.num_frames() , num=traj1.num_frames())

traj2=Trajectory(sim2_dcd,sim2_psf,stride=1,first=0,last=5001)
time_line2=np.linspace(0,traj2.num_frames() , num=traj2.num_frames())


traj3=Trajectory(sim3_dcd,sim3_psf,stride=1)
time_line3=np.linspace(0,traj3.num_frames() , num=traj3.num_frames())

traj4=Trajectory(sim4_dcd,sim4_psf,stride=1)

traj5=Trajectory(sim5_dcd,sim5_psf,stride=1)
"""
print (traj3.get_molid())
#contact =traj.porcentaje_contact_fit_list("protein","name P",3000,5000)
"""
r1=traj1.membrane_center_mass("protein","resname POPG POPC")
traj1.close()
r2=traj2.membrane_center_mass("protein","resname POPG POPC")
traj2.close()
r3=traj3.membrane_center_mass("protein","resname POPG POPC")
#r33=Trajectory.delete_frame_wrap(r3,1)
traj3.close()

r4=traj4.membrane_center_mass("protein","resname POPG POPC")
time_line4=np.linspace(0,traj4.num_frames() , num=traj4.num_frames())

traj4.close()
r5=traj5.membrane_center_mass("protein","resname POPG POPC")
traj5.close()


plt.plot(time_line, r1, label='XI POPG')
plt.plot(time_line2, r2, label='XI POPG:POPC')
plt.plot(time_line3, r3, label='XI 3POPG:POPC')
#plt.plot(time_line3, r33, label='sim3f')
plt.plot(time_line4, r4, label='XI 3POPC:POPG')
plt.plot(time_line4, r5, label='XI POPC')
plt.legend()
plt.savefig("/media/eniac/mdd1/paper_membranas/xi_jinghao_toxin/distance_pOPg_XI_toxin.png", dpi=900) 
plt.show()







# %%
