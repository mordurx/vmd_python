
#%%%
import plotly.graph_objects as go
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

    
sim1_dcd= "/mnt/e/3popc_popg/rep0.mem.tox.dcd"
sim1_psf="/mnt/e/3popc_popg/rep0.mem.tox.psf"

sim2_dcd= "/mnt/e/3popc_popg/rep1_tox_mem500.dcd"
sim2_psf="/mnt/e/3popc_popg/rep1.mem.tox.psf"

sim3_dcd= "/mnt/e/3popc_popg/rep2_tox_mem500.dcd"
sim3_psf="/mnt/e/3popc_popg/rep2.mem.tox.psf"

sim4_dcd= "/mnt/e/3popc_popg/rep3_tox_mem500.dcd"
sim4_psf="/mnt/e/3popc_popg/rep3.mem.tox.psf"

traj1=Trajectory(sim1_dcd,sim1_psf,stride=1)
traj2=Trajectory(sim2_dcd,sim2_psf,stride=1)
traj3=Trajectory(sim3_dcd,sim3_psf,stride=1)
traj4=Trajectory(sim4_dcd,sim4_psf,stride=1)
print (traj1.get_molid())
time_line=np.linspace(0, 41, num=41)

r1 =traj1.porcentaje_contact_fit("protein","name P",3000,5000)
r2 =traj2.porcentaje_contact_fit("protein","name P",3000,5000)
r3 =traj3.porcentaje_contact_fit("protein","name P",3000,5000)
r4 =traj4.porcentaje_contact_fit("protein","name P",3000,5000)

traj1.close()
traj2.close()
traj3.close()
traj4.close()

colors=color_snx()
index = np.arange(1,42)
#plt.rcParams['font.family'] = 'sans-serif'
#plt.rcParams['font.sans-serif'] = 'Comic Sans MS'

fig,a =  plt.subplots(2,2)
fig.tight_layout()
lines=a[0][0].bar(index, r1,color=colors,edgecolor='gray',width = (index[1]-index[0])*0.8)
a[0][0].set_title('sim 1 300-500 ns')
a[0][0].set_xticks(index)
a[0][0].set_xticklabels(index,fontsize=4,rotation=90)
a[0][0].set_ylabel('occurrences')
#a[0][0].set_xlabel('residues')
a[0][0].legend(lines, ['polar', 'no polar', 'acid','basic'], loc='upper center',bbox_to_anchor=(0.22, 0.5, 0.5, 0.5),fontsize=5)

#a.tick_params(axis="x",rotation=0,labelsize=20)
#lines=plt.bar(index, r1,color=colors,edgecolor='gray',width = (index[1]-index[0])*0.8)
lines=a[0][1].bar(index, r2,color=colors,edgecolor='gray',width = (index[1]-index[0])*0.8)
a[0][1].set_title('sim 2 300-500 ns')
a[0][1].set_xticks(index)
a[0][1].set_xticklabels(index,fontsize=4,rotation=90)
#a[0][1].set_ylabel('occurrences')
#a[0][0].set_xlabel('residues')
a[0][1].legend(lines, ['polar', 'no polar', 'acid','basic'], loc='upper center',bbox_to_anchor=(0.22, 0.5, 0.5, 0.5),fontsize=5)

lines=a[1][0].bar(index, r3,color=colors,edgecolor='gray',width = (index[1]-index[0])*0.8)
a[1][0].set_title('sim 3 300-500 ns')
a[1][0].set_xticks(index)
a[1][0].set_xticklabels(index,fontsize=4,rotation=90)
#a[0][1].set_ylabel('occurrences')
a[1][0].set_xlabel('residues')
a[1][0].legend(lines, ['polar', 'no polar', 'acid','basic'], loc='upper center',bbox_to_anchor=(0.22, 0.5, 0.5, 0.5),fontsize=5)

lines=a[1][1].bar(index, r4,color=colors,edgecolor='gray',width = (index[1]-index[0])*0.8)
a[1][1].set_title('sim 4 300-500 ns')
a[1][1].set_xticks(index)
a[1][1].set_xticklabels(index,fontsize=4,rotation=90)
a[1][1].set_ylabel('occurrences')
a[1][1].set_xlabel('residues')
a[1][1].legend(lines, ['polar', 'no polar', 'acid','basic'], loc='upper center',bbox_to_anchor=(0.22, 0.5, 0.5, 0.5),fontsize=5)


result=Trajectory.promedio_ocurrencias(r1,r2,r3,r4)
#plt.plot(time_line, r1, label='sim1')
#plt.plot(time_line, r2, label='sim2')
#plt.plot(time_line, r3, label='sim3')
#plt.plot(time_line, r4, label='sim4')
#plt.legend()
plt.savefig("/mnt/e/3popc_popg/residue.png", dpi=900)
plt.figure()
lower_error = np.divide(result[1], 5000000)
upper_error = result[1]
asymmetric_error = [lower_error, upper_error]
#plt.errorbar(index, result[0], result[1], linestyle='None', marker='^')
plt.bar(index, result[0],color=colors,edgecolor='gray',width = (index[1]-index[0])*0.8,yerr=asymmetric_error,capsize=2)
plt.xticks(index, index, fontsize=5, rotation=0)
plt.title('POPC promedio 300-500 ns')
plt.ylabel('relative ocurrences')
plt.xlabel('residue #')
plt.legend(lines, ['polar', 'no polar', 'acid','basic'], loc='upper center',bbox_to_anchor=(0.22, 0.5, 0.5, 0.5),fontsize=7)
plt.savefig("/mnt/e/3popc_popg/12.png", dpi=900) 
plt.show()

#### colorear por beta factor
pdb='/mnt/e/3popc_popg/snx482.pdb'
output='/mnt/e/3popc_popg/snx3popc1popg.pdb'
molid = molecule.load("pdb", pdb)
print(molid)
query = "segname TOX"
#vector=np.linspace(0, 1, num=len(protein.beta),)
Trajectory.set_beta_factor(result[0],query,molid,output)
#color scale
# small=red, middle=white, large=blue


# %%
