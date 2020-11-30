from trajectory import Trajectory
import matplotlib.pyplot as plt
import numpy as np
from vmd import molecule,atomsel
from matplotlib  import rcParams
def color_snx():
    dict_color_snx_482={}
    dict_color_snx_482['1']='green'
    dict_color_snx_482['2']='white'
    dict_color_snx_482['3']='red'
    dict_color_snx_482['4']='blue'
    dict_color_snx_482['5']='white'
    dict_color_snx_482['6']='green'
    dict_color_snx_482['7']='white'
    dict_color_snx_482['8']='blue'
    dict_color_snx_482['9']='green'
    dict_color_snx_482['10']='white'
    dict_color_snx_482['11']='white'
    dict_color_snx_482['12']='green'
    dict_color_snx_482['13']='green'
    dict_color_snx_482['14']='white'
    dict_color_snx_482['15']='green'
    dict_color_snx_482['16']='white'
    dict_color_snx_482['17']='green'
    dict_color_snx_482['18']='red'
    dict_color_snx_482['19']='red'
    dict_color_snx_482['20']='white'
    dict_color_snx_482['21']='white'
    dict_color_snx_482['22']='white'
    dict_color_snx_482['23']='blue'
    dict_color_snx_482['24']='white'
    dict_color_snx_482['25']='green'
    dict_color_snx_482['26']='white'
    dict_color_snx_482['27']='green'
    dict_color_snx_482['28']='green'
    dict_color_snx_482['29']='white'
    dict_color_snx_482['30']='white'
    dict_color_snx_482['31']='green'
    dict_color_snx_482['32']='green'
    dict_color_snx_482['33']='white'
    dict_color_snx_482['34']='white'
    dict_color_snx_482['35']='white'
    dict_color_snx_482['36']='red'
    dict_color_snx_482['37']='white'
    dict_color_snx_482['38']='green'
    dict_color_snx_482['39']='white'
    dict_color_snx_482['40']='green'
    dict_color_snx_482['41']='red'
    return list(dict_color_snx_482.values())

    """ POPG"""
sim1_dcd='/mnt/e/paper_membranas/POPG/POPG_resume/unwrap_sim1_popg.dcd'
sim1_psf='/mnt/e/paper_membranas/POPG/POPG_resume/sim1_popg.psf'

sim2_dcd='/mnt/e/paper_membranas/POPG/POPG_resume/unwrap_sim2_popg.dcd'
sim2_psf='/mnt/e/paper_membranas/POPG/POPG_resume/sim2_popg_ions.Wat.psf'

sim3_dcd= "/mnt/e/paper_membranas/POPG/POPG_resume/unwrap_sim3_popg.dcd"
sim3_psf="/mnt/e/paper_membranas/POPG/POPG_resume/sim3_POPG_ions.Wat.psf"

sim4_dcd='/mnt/e/paper_membranas/POPG/POPG_resume/sim4_POPG_filt_memtox.dcd'
sim4_psf='/mnt/e/paper_membranas/POPG/POPG_resume/sim4_POPG_POPC.ions.Wat.psf'

traj1=Trajectory(sim1_dcd,sim1_psf,stride=1)
traj2=Trajectory(sim2_dcd,sim2_psf,stride=1)
traj3=Trajectory(sim3_dcd,sim3_psf,stride=1)
traj4=Trajectory(sim4_dcd,sim4_psf,stride=1)
time_line=np.linspace(0, 41, num=41)

r1_POPG =traj1.porcentaje_contact_fit("protein","name P",3000,5000)
r2_POPG =traj2.porcentaje_contact_fit("protein","name P",3000,5000)
r3_POPG =traj3.porcentaje_contact_fit("protein","name P",3000,5000)
r4_POPG =traj4.porcentaje_contact_fit("protein","name P",3000,5000)
traj1.close()
traj2.close()
traj3.close()
traj4.close()


""" 3POPG"""
sim1_3popg_dcd= "/mnt/e/paper_membranas/3POPG_popc_resume/sim1_3popg_popc_unwrap.dcd"
sim1_3popg_psf="/mnt/e/paper_membranas/3POPG_popc_resume/sim1_3popg_popc_unwrap.psf"

sim2_3popg_dcd= "/mnt/e/paper_membranas/3POPG_popc_resume/sim2_3popg_popc_unwrap.dcd"
sim2_3popg_psf="/mnt/e/paper_membranas/3POPG_popc_resume/sim2_3popg_popc_unwrap.psf"

sim3_3popg_dcd= "/mnt/e/paper_membranas/3POPG_popc_resume/sim3_3popg_popc.dcd"
sim3_3popg_psf="/mnt/e/paper_membranas/3POPG_popc_resume/sim3_3popg_popc.psf"

sim4_popg_dcd= "/mnt/e/paper_membranas/3POPG_popc_resume/sim4_3popg_popc.dcd"
sim4_popg_psf="/mnt/e/paper_membranas/3POPG_popc_resume/sim4_3popg_popc.psf"


traj1_3popg=Trajectory(sim1_3popg_dcd,sim1_3popg_psf,stride=1)
traj2_3popg=Trajectory(sim2_3popg_dcd,sim2_3popg_psf,stride=1)
traj3_3popg=Trajectory(sim3_3popg_dcd,sim3_3popg_psf,stride=1)
traj4_3popg=Trajectory(sim4_popg_dcd,sim4_popg_psf,stride=1)

r1_3POPG =traj1_3popg.porcentaje_contact_fit("protein","name P",3000,5000)
r2_3POPG =traj2_3popg.porcentaje_contact_fit("protein","name P",3000,5000)
r3_3POPG =traj3_3popg.porcentaje_contact_fit("protein","name P",3000,5000)
r4_3POPG =traj4_3popg.porcentaje_contact_fit("protein","name P",3000,5000)

traj1_3popg.close()
traj2_3popg.close()
traj3_3popg.close()
traj4_3popg.close()

"""POPC"""
sim1_popc_dcd= "/mnt/e/paper_membranas/POPC/sim_POPC_resume/sim1_unwrap_700nspopc.dcd"
sim1_popc_psf="/mnt/e/paper_membranas/POPC/sim_POPC_resume/sim1_unwrap_700nspopc.psf"

sim2_popc_dcd= "/mnt/e/paper_membranas/POPC/sim_POPC_resume/sim2_unwrap_700nspopc.dcd"
sim2_popc_psf="/mnt/e/paper_membranas/POPC/sim_POPC_resume/sim2_unwrap_700nspopc.psf"

traj1_popc=Trajectory(sim1_popc_dcd,sim1_popc_psf,stride=1)
traj2_popc=Trajectory(sim2_popc_dcd,sim2_popc_psf,stride=1)

r1_popc =traj1_popc.porcentaje_contact_fit("protein","name P",3000,5000)
r2_popc =traj2_popc.porcentaje_contact_fit("protein","name P",3000,5000)
traj1_popc.close()
traj2_popc.close()


"""1POPC:1POPG"""
sim1_popc_popg_dcd= "/mnt/e/paper_membranas/1POPC_1POPG/sim_POPC_POPG_resume/sim1_popc_popg_unwrap.dcd"
sim1_popc_popg_psf="/mnt/e/paper_membranas/1POPC_1POPG/sim_POPC_POPG_resume/sim1_popc_popg_unwrap.psf"

sim2_popc_popg_dcd= "/mnt/e/paper_membranas/1POPC_1POPG/sim_POPC_POPG_resume/sim2_popc_popg_unwrap.dcd"
sim2_popc_popg_psf="/mnt/e/paper_membranas/1POPC_1POPG/sim_POPC_POPG_resume/sim2_popc_popg_unwrap.psf"

sim3_popc_popg_dcd= "/mnt/e/paper_membranas/1POPC_1POPG/sim_POPC_POPG_resume/sim3_unwrap_popc_popg.dcd"
sim3_popc_popg_psf="/mnt/e/paper_membranas/1POPC_1POPG/sim_POPC_POPG_resume/sim3_unwrap_popc_popg.psf"

sim4_popc_popg_dcd= "/mnt/e/paper_membranas/1POPC_1POPG/sim_POPC_POPG_resume/sim4_unwrap_popc_popg.dcd"
sim4_popc_popg_psf="/mnt/e/paper_membranas/1POPC_1POPG/sim_POPC_POPG_resume/sim4_unwrap_popc_popg.psf"

traj1_popc_popg=Trajectory(sim1_popc_popg_dcd,sim1_popc_popg_psf,stride=1)
traj2_popc_popg=Trajectory(sim2_popc_popg_dcd,sim2_popc_popg_psf,stride=1)
traj3_popc_popg=Trajectory(sim3_popc_popg_dcd,sim3_popc_popg_psf,stride=1)
traj4_popc_popg=Trajectory(sim4_popc_popg_dcd,sim4_popc_popg_psf,stride=1)

r1_popc_popg =traj1_popc_popg.porcentaje_contact_fit("protein","name P",3000,5000)
r2_popc_popg =traj2_popc_popg.porcentaje_contact_fit("protein","name P",3000,5000)
r3_popc_popg =traj3_popc_popg.porcentaje_contact_fit("protein","name P",3000,5000)
r4_popc_popg =traj4_popc_popg.porcentaje_contact_fit("protein","name P",3000,5000)

traj1_popc_popg.close()
traj2_popc_popg.close()
traj3_popc_popg.close()
traj4_popc_popg.close()



"""ploteo"""

colors=color_snx()
index = np.arange(1,42)
#plt.rcParams['font.family'] = 'sans-serif'
#plt.rcParams['font.sans-serif'] = 'Comic Sans MS'

fig,a =  plt.subplots(2,2)
fig.tight_layout()
lines=a[0][0].bar(index, r1_POPG,color=colors,edgecolor='gray',width = (index[1]-index[0])*0.8)
a[0][0].set_title('sim 1 300-500 ns')
a[0][0].set_xticks(index)
a[0][0].set_xticklabels(index,fontsize=4,rotation=90)
#a[0][0].set_ylabel('occurrences')
#a[0][0].set_xlabel('residues')
a[0][0].legend(lines, ['polar', 'no polar', 'acid','basic'], loc='upper center',bbox_to_anchor=(0.22, 0.5, 0.5, 0.5),fontsize=5)

plt.savefig("/mnt/e/paper_membranas/POPG/POPG_resume/popg_residues_insert.png", dpi=900)