
#%%%
from trajectory import Trajectory
import numpy as np
import matplotlib.pyplot as plt
#Latex Backend for Matplotlib
plt.rc('text', usetex=True)  
plt.rc('font', family='serif') 


# #3popg-popc

# sim1_dcd= "/home/guido/resume_paper_membranas/3POPG_popc_resume/sim1_3popg_popc_unwrap.dcd"
# sim1_psf="/home/guido/resume_paper_membranas/3POPG_popc_resume/sim1_3popg_popc_unwrap.psf"

# sim2_dcd= "/home/guido/resume_paper_membranas/3POPG_popc_resume/sim2_3popg_popc_unwrap.dcd"
# sim2_psf="/home/guido/resume_paper_membranas/3POPG_popc_resume/sim2_3popg_popc_unwrap.psf"

# sim3_dcd= "/home/guido/resume_paper_membranas/3POPG_popc_resume/sim3_3popg_popc.dcd"
# sim3_psf="/home/guido/resume_paper_membranas/3POPG_popc_resume/sim3_3popg_popc.psf"
# sim4_dcd= "/home/guido/resume_paper_membranas/3POPG_popc_resume/sim4_3popg_popc.dcd"
# sim4_psf="/home/guido/resume_paper_membranas/3POPG_popc_resume/sim4_3popg_popc.psf"
# img_distance="/home/guido/resume_paper_membranas/3POPG_popc_resume/distance_pOPg_XI_toxin.png"
# label=['sim1 3POPG:POPC','sim2 3POPG:POPC','sim3 3POPG:POPC','sim4 3POPG:POPC']

"""# Jingzhaotoxin-XI (JZTX-XI)

sim1_dcd='/home/guido/resume_paper_membranas/IX_TOXIN_resume/unwrap_XI_sim1_POPG.dcd'
sim1_psf='/home/guido/resume_paper_membranas/IX_TOXIN_resume/unwrap_XI_sim1_POPG.psf'

sim2_dcd='/home/guido/resume_paper_membranas/IX_TOXIN_resume/unwrap_POPC_POPG_sim1.dcd'
sim2_psf='/home/guido/resume_paper_membranas/IX_TOXIN_resume/unwrap_POPC_POPG_sim1.psf'

sim3_dcd= "/home/guido/resume_paper_membranas/IX_TOXIN_resume/unwrap_3popg_popg_XI.dcd"
sim3_psf="/home/guido/resume_paper_membranas/IX_TOXIN_resume/unwrap_3popg_popg_XI.psf"

sim4_dcd='/home/guido/resume_paper_membranas/IX_TOXIN_resume/unwrap_XI_3POPC_POPG.dcd'
sim4_psf='/home/guido/resume_paper_membranas/IX_TOXIN_resume/unwrap_XI_3POPC_POPG.psf'

sim5_dcd='/home/guido/resume_paper_membranas/IX_TOXIN_resume/unwrap_XI_sim3_POPC.dcd'
sim5_psf='/home/guido/resume_paper_membranas/IX_TOXIN_resume/unwrap_XI_sim3_POPC.psf'
img_distance="/home/guido/resume_paper_membranas/IX_TOXIN_resume/distance_pOPg_XI_toxin.png"

label=['JZTX-XI POPG','JZTX-XI POPC:POPG','JZTX-XI 3POPG:POPC','JZTX-XI 3POPC:POPG','JZTX-XI POPC']
"""
#popg
# sim1_dcd='/home/guido/resume_paper_membranas/POPG_resume/unwrap_sim1_popg.dcd'
# sim1_psf='/home/guido/resume_paper_membranas/POPG_resume/unwrap_sim1_popg.psf'

# sim2_dcd='/home/guido/resume_paper_membranas/POPG_resume/unwrap_sim2_popg.dcd'
# sim2_psf='/home/guido/resume_paper_membranas/POPG_resume/unwrap_sim2_popg.psf'

# sim3_dcd= "/home/guido/resume_paper_membranas/POPG_resume/unwrap_sim3_popg.dcd"
# sim3_psf="/home/guido/resume_paper_membranas/POPG_resume/unwrap_sim3_popg.psf"

# sim4_dcd='/home/guido/resume_paper_membranas/POPG_resume/sim4_POPG_filt_memtox.dcd'
# sim4_psf='/home/guido/resume_paper_membranas/POPG_resume/sim4_POPG_filt_memtox.psf'
# img_distance="/home/guido/resume_paper_membranas/POPG_resume/distance_pOPg__toxin.png"

# label=['sim1 POPG','sim2 POPG','sim3 POPG','sim4 POPG']


# """POPC"""
# sim1_dcd= "/home/guido/resume_paper_membranas/sim_POPC_resume/sim1_unwrap_700nspopc.dcd"
# sim1_psf="/home/guido/resume_paper_membranas/sim_POPC_resume/sim1_unwrap_700nspopc.psf"

# sim2_dcd= "/home/guido/resume_paper_membranas/sim_POPC_resume/sim2_unwrap_700nspopc.dcd"
# sim2_psf="/home/guido/resume_paper_membranas/sim_POPC_resume/sim2_unwrap_700nspopc.psf"

# sim3_dcd= "/home/guido/resume_paper_membranas/sim_POPC_resume/sim3_unwrap_500ns.dcd"
# sim3_psf="/home/guido/resume_paper_membranas/sim_POPC_resume/sim3_unwrap_500ns.psf"

# sim4_dcd='/home/guido/resume_paper_membranas/sim_POPC_resume/sim4_unwrap_500ns.dcd'
# sim4_psf='/home/guido/resume_paper_membranas/sim_POPC_resume/sim4_unwrap_500ns.psf'

# img_distance="/home/guido/resume_paper_membranas/sim_POPC_resume/sim_POPC_resume.png"
# label=['sim1 POPC','sim2 POPC','sim3 POPC','sim4 POPC']

# """3POPC:1POPG"""
# sim1_dcd= "/home/guido/resume_paper_membranas/3POPC_popg_resume/sim1_3popc_popg_unwrap.dcd"
# sim1_psf="/home/guido/resume_paper_membranas/3POPC_popg_resume/sim1_3popc_popg_unwrap.psf"

# sim2_dcd= "/home/guido/resume_paper_membranas/3POPC_popg_resume/sim2_3popc_popg_unwrap.dcd"
# sim2_psf="/home/guido/resume_paper_membranas/3POPC_popg_resume/sim2_3popc_popg_unwrap.psf"

# sim3_dcd= "/home/guido/resume_paper_membranas/3POPC_popg_resume/sim3_3popc_popg_unwrap.dcd"
# sim3_psf="/home/guido/resume_paper_membranas/3POPC_popg_resume/sim3_3popc_popg_unwrap.psf"

# sim4_dcd= "/home/guido/resume_paper_membranas/3POPC_popg_resume/sim4_3popc_popg_filt_unwrap.dcd.dcd"
# sim4_psf="/home/guido/resume_paper_membranas/3POPC_popg_resume/sim4_3popc_popg_filt_unwrap.dcd.psf"


# img_distance="/home/guido/resume_paper_membranas/3POPC_popg_resume/3POPC_popg_resume.png"
# label=['sim1 3POPC:POPG','sim2 3POPC:POPG','sim3 3POPC:POPG','sim4 3POPC:POPG']

"""1POPC:1POPG"""
sim1_dcd= "/home/guido/resume_paper_membranas/sim_POPC_POPG_resume/sim1_popc_popg_unwrap.dcd"
sim1_psf="/home/guido/resume_paper_membranas/sim_POPC_POPG_resume/sim1_popc_popg_unwrap.psf"

sim2_dcd= "/home/guido/resume_paper_membranas/sim_POPC_POPG_resume/sim2_popc_popg_unwrap.dcd"
sim2_psf="/home/guido/resume_paper_membranas/sim_POPC_POPG_resume/sim2_popc_popg_unwrap.psf"

sim3_dcd= "/home/guido/resume_paper_membranas/sim_POPC_POPG_resume/sim3_unwrap_popc_popg.dcd"
sim3_psf="/home/guido/resume_paper_membranas/sim_POPC_POPG_resume/sim3_unwrap_popc_popg.psf"

sim4_dcd= "/home/guido/resume_paper_membranas/sim_POPC_POPG_resume/sim4_unwrap_popc_popg.dcd"
sim4_psf="/home/guido/resume_paper_membranas/sim_POPC_POPG_resume/sim4_unwrap_popc_popg.psf"
img_distance="/home/guido/resume_paper_membranas/sim_POPC_POPG_resume/POPC_popg_resume.png"
label=['sim1 POPC:POPG','sim2 POPC:POPG','sim3 POPC:POPG','sim4 POPC:POPG']


traj1=Trajectory(sim1_dcd,sim1_psf,stride=10,last=5000)
time_line=np.linspace(0,traj1.num_frames() , num=traj1.num_frames())

traj2=Trajectory(sim2_dcd,sim2_psf,stride=10,last=5000)
time_line2=np.linspace(0,traj2.num_frames() , num=traj2.num_frames())


traj3=Trajectory(sim3_dcd,sim3_psf,stride=10)
time_line3=np.linspace(0,traj3.num_frames() , num=traj3.num_frames())

traj4=Trajectory(sim4_dcd,sim4_psf,stride=10)

#traj5=Trajectory(sim5_dcd,sim5_psf,stride=10)
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
#r5=traj5.membrane_center_mass("protein","resname POPG POPC")
#traj5.close()
#plt.figure(figsize=(3.6, 3.6))
#%%%
plt.plot(time_line, r1, label=label[0])
plt.plot(time_line2, r2, label=label[1])
plt.plot(time_line3, r3, label=label[2])
#plt.plot(time_line3, r33, label='sim3f')
plt.plot(time_line4, r4, label=label[3])
#plt.plot(time_line4, r5, label='JZTX-XI POPC')
plt.legend(prop={'size': 7})
plt.xlabel(r'\textbf{time (ns)}')
plt.ylabel(r'\textbf{COM distance (\AA{})}')
#plt.title('Legend inside')
#plt.title(r"\TeX\ is Number $\displaystyle\sum_{n=1}^\infty\frac{-e^{i\pi}}{2^n}$!", fontsize=16, color='r')
plt.savefig(img_distance, dpi=900) 
plt.show()







# %%
