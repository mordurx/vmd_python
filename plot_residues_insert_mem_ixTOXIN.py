
#%%%

from trajectory import Trajectory
import matplotlib.pyplot as plt
import numpy as np
from vmd import molecule,atomsel
from matplotlib  import rcParams
def color_IX_toxin():
    dict_color_IX_toxin_482={}
    dict_color_IX_toxin_482['1']='red'
    dict_color_IX_toxin_482['2']='white'
    dict_color_IX_toxin_482['3']='blue'
    dict_color_IX_toxin_482['4']='blue'
    dict_color_IX_toxin_482['5']='white'
    dict_color_IX_toxin_482['6']='white'
    dict_color_IX_toxin_482['7']='green'
    dict_color_IX_toxin_482['8']='green'
    dict_color_IX_toxin_482['9']='white'
    dict_color_IX_toxin_482['10']='green'
    dict_color_IX_toxin_482['11']='white'
    dict_color_IX_toxin_482['12']='red'
    dict_color_IX_toxin_482['13']='green'
    dict_color_IX_toxin_482['14']='red'
    dict_color_IX_toxin_482['15']='white'
    dict_color_IX_toxin_482['16']='white'
    dict_color_IX_toxin_482['17']='white'
    dict_color_IX_toxin_482['18']='green'
    dict_color_IX_toxin_482['19']='white'
    dict_color_IX_toxin_482['20']='green'
    dict_color_IX_toxin_482['21']='white'
    dict_color_IX_toxin_482['22']='blue'
    dict_color_IX_toxin_482['23']='white'
    dict_color_IX_toxin_482['24']='green'
    dict_color_IX_toxin_482['25']='white'
    dict_color_IX_toxin_482['26']='blue'
    dict_color_IX_toxin_482['27']='green'
    dict_color_IX_toxin_482['28']='white'
    dict_color_IX_toxin_482['29']='white'
    dict_color_IX_toxin_482['30']='white'
    dict_color_IX_toxin_482['31']='red'
    dict_color_IX_toxin_482['32']='green'
    dict_color_IX_toxin_482['33']='green'
    dict_color_IX_toxin_482['34']='white'
 
    return list(dict_color_IX_toxin_482.values())

"""    
sim1_dcd= "/media/eniac/mdd1/paper_membranas/sim1_3popg_1popc/sim1_3popg_popc_fixed_unwrapped.dcd"
sim1_psf="/media/eniac/mdd1/paper_membranas/sim1_3popg_1popc/agua_tox_mem_Sim1_ions.Wat.psf"

sim2_dcd= "/media/eniac/mdd1/paper_membranas/sim2_3popg_popc/sim2_3popg_unwrap.dcd"
sim2_psf="/media/eniac/mdd1/paper_membranas/sim2_3popg_popc/3popg_popc_tox.psf"

sim3_dcd= "/media/eniac/mdd1/paper_membranas/sim3_3popg_popc/MD/sim3_3popg_popc.dcd"
sim3_psf="/media/eniac/mdd1/paper_membranas/sim3_3popg_popc/agua_tox_mem_Sim1_ions.Wat.psf"
sim4_dcd= "/media/eniac/mdd1/paper_membranas/sim4_3popg_popc/sim4_3popg_popc.dcd"
sim4_psf="/media/eniac/mdd1/paper_membranas/sim4_3popg_popc/agua_tox_mem_Sim1_ions.Wat.psf"
"""



sim1_dcd='/media/eniac/mdd1/paper_membranas/xi_jinghao_toxin/XI_POPG/unwrap_XI_sim1_POPG.dcd'
sim1_psf='/media/eniac/mdd1/paper_membranas/xi_jinghao_toxin/XI_POPG/unwrap_XI_sim1_POPG.psf'

sim2_dcd='/media/eniac/mdd1/paper_membranas/xi_jinghao_toxin/XI_POPC_POPG/unwrap_POPC_POPG_sim1.dcd'
sim2_psf='/media/eniac/mdd1/paper_membranas/xi_jinghao_toxin/XI_POPC_POPG/unwrap_POPC_POPG_sim1.psf'

sim3_dcd='/media/eniac/mdd1/paper_membranas/xi_jinghao_toxin/XI_3POPG_POPC/unwrapped_3popg_popg_XI.dcd'
sim3_psf='/media/eniac/mdd1/paper_membranas/xi_jinghao_toxin/XI_3POPG_POPC/unwrapped_3popg_popg_XI.psf'

sim4_dcd= "/media/eniac/mdd1/paper_membranas/xi_jinghao_toxin/XI_3POPC_POPG/unwrap_XI_3POPC_POPG.dcd"
sim4_psf="/media/eniac/mdd1/paper_membranas/xi_jinghao_toxin/XI_3POPC_POPG/unwrap_XI_3POPC_POPG.psf"


sim5_dcd='/media/eniac/mdd1/paper_membranas/xi_jinghao_toxin/XI_POPC/unwrap_XI_sim3_POPC.dcd'
sim5_psf='/media/eniac/mdd1/paper_membranas/xi_jinghao_toxin/XI_POPC/unwrap_XI_sim3_POPC.psf'






"""
sim1_dcd='/media/eniac/mdd1/paper_membranas/POPG/unwrap_sim1_popg.dcd'
sim1_psf='/media/eniac/mdd1/paper_membranas/POPG/sim1_popg.psf'

sim2_dcd='/media/eniac/mdd1/paper_membranas/POPG/unwrap_sim2_popg.dcd'
sim2_psf='/media/eniac/mdd1/paper_membranas/POPG/sim2_popg_ions.Wat.psf'

sim3_dcd= "/media/eniac/mdd1/paper_membranas/POPG/unwrap_sim3_popg.dcd"
sim3_psf="/media/eniac/mdd1/paper_membranas/POPG/sim3_POPG_ions.Wat.psf"

sim4_dcd='/media/eniac/mdd1/paper_membranas/POPG/sim4_POPG_filt_memtox.dcd'
sim4_psf='/media/eniac/mdd1/paper_membranas/POPG/sim4_POPG_POPC.ions.Wat.psf'
"""
traj1=Trajectory(sim1_dcd,sim1_psf,stride=1)
traj2=Trajectory(sim2_dcd,sim2_psf,stride=1)
traj3=Trajectory(sim3_dcd,sim3_psf,stride=1)
traj4=Trajectory(sim4_dcd,sim4_psf,stride=1)
traj5=Trajectory(sim5_dcd,sim5_psf,stride=1)

print (traj1.get_molid())
#time_line=np.linspace(0, 41, num=41)

r1 =traj1.porcentaje_contact_fit("protein","name P",3000,5000)
r2 =traj2.porcentaje_contact_fit("protein","name P",3000,5000)
r3 =traj3.porcentaje_contact_fit("protein","name P",3000,5000)
r4 =traj4.porcentaje_contact_fit("protein","name P",3000,5000)
r5 =traj5.porcentaje_contact_fit("protein","name P",3000,5000)

traj1.close()
traj2.close()
traj3.close()
traj4.close()
traj5.close()


colors=color_IX_toxin()
index = np.arange(1,35)
#plt.rcParams['font.family'] = 'sans-serif'
#plt.rcParams['font.sans-serif'] = 'Comic Sans MS'

fig,a =  plt.subplots(2,3)
#fig.tight_layout()
lines=a[0][0].bar(index, r1,color=colors,edgecolor='gray',width = (index[1]-index[0])*0.8)
a[0][0].set_title('sim POPG IX 300-500 ns')
a[0][0].set_xticks(index)
a[0][0].set_xticklabels(index,fontsize=8,rotation=90)
a[0][0].set_ylabel('occurrences')
a[0][0].set_ylim(0, 1)
#a[0][0].set_xlabel('residues')
a[0][0].legend(lines, ['polar', 'no polar', 'acid','basic'], loc='upper center',bbox_to_anchor=(0.22, 0.5, 0.5, 0.5),fontsize=10)
#a.tick_params(axis="x",rotation=0,labelsize=20)
#lines=plt.bar(index, r1,color=colors,edgecolor='gray',width = (index[1]-index[0])*0.8)
lines=a[0][1].bar(index, r2,color=colors,edgecolor='gray',width = (index[1]-index[0])*0.8)
a[0][1].set_title('sim 2 POPG_POPC IX 300-500 ns')
a[0][1].set_xticks(index)
a[0][1].set_ylim(0, 1) 
a[0][1].set_xticklabels(index,fontsize=8,rotation=90)
a[0][1].legend(lines, ['polar', 'no polar', 'acid','basic'], loc='upper center',bbox_to_anchor=(0.22, 0.5, 0.5, 0.5),fontsize=10)

#a[0][1].set_ylabel('occurrences')
#a[0][0].set_xlabel('residues')
#a[0][1].legend(lines, ['polar', 'no polar', 'acid','basic'], loc='upper center',bbox_to_anchor=(0.22, 0.5, 0.5, 0.5),fontsize=5)

lines=a[1][0].bar(index, r3,color=colors,edgecolor='gray',width = (index[1]-index[0])*0.8)
a[1][0].set_title('sim 3 3POPG_POPC IX 300-500 ns')
a[1][0].set_xticks(index)
a[1][0].legend(lines, ['polar', 'no polar', 'acid','basic'], loc='upper center',bbox_to_anchor=(0.22, 0.5, 0.5, 0.5),fontsize=10)
a[1][0].set_ylim(0, 1) 
a[1][0].set_xticklabels(index,fontsize=8,rotation=90)
a[1][0].set_ylabel('occurrences')
a[1][0].set_xlabel('residues')
#a[1][0].legend(lines, ['polar', 'no polar', 'acid','basic'], loc='upper center',bbox_to_anchor=(0.22, 0.5, 0.5, 0.5),fontsize=5)

lines=a[1][1].bar(index, r4,color=colors,edgecolor='gray',width = (index[1]-index[0])*0.8)
a[1][1].set_title('sim 4 3POPC_POPG IX  300-500 ns')
a[1][1].set_xticks(index)
a[1][1].set_xticklabels(index,fontsize=8,rotation=90)
#a[1][1].set_ylabel('occurrences')
a[1][1].set_xlabel('residues')
#a[1][1].legend(lines, ['polar', 'no polar', 'acid','basic'], loc='upper center',bbox_to_anchor=(0.22, 0.5, 0.5, 0.5),fontsize=5)
a[1][1].set_ylim(0, 1) 
a[1][1].legend(lines, ['polar', 'no polar', 'acid','basic'], loc='upper center',bbox_to_anchor=(0.22, 0.5, 0.5, 0.5),fontsize=10)

lines=a[0][2].bar(index, r5,color=colors,edgecolor='gray',width = (index[1]-index[0])*0.8)
a[0][2].set_title('sim 5 POPC IX  300-500 ns')
a[0][2].set_xticks(index)
a[0][2].set_xticklabels(index,fontsize=8,rotation=90)
#a[1][1].set_ylabel('occurrences')
#a[0][2].set_xlabel('residues')
a[0][2].legend(lines, ['polar', 'no polar', 'acid','basic'], loc='upper center',bbox_to_anchor=(0.22, 0.5, 0.5, 0.5),fontsize=10)

#a[1][1].legend(lines, ['polar', 'no polar', 'acid','basic'], loc='upper center',bbox_to_anchor=(0.22, 0.5, 0.5, 0.5),fontsize=5)

result=Trajectory.promedio_ocurrencias(r1,r2,r3,r4,r5)
#plt.plot(time_line, r1, label='sim1')
#plt.plot(time_line, r2, label='sim2')
#plt.plot(time_line, r3, label='sim3')
#plt.plot(time_line, r4, label='sim4')
#plt.legend()
plt.savefig("/media/eniac/mdd1/paper_membranas/xi_jinghao_toxin/residues_insert_IX.png", dpi=900)
lower_error = np.divide(result[1], 5000000)
upper_error = result[1]
asymmetric_error = [lower_error, upper_error]
#plt.errorbar(index, result[0], result[1], linestyle='None', marker='^')
lines=a[1][2].bar(index, result[0],color=colors,edgecolor='gray',width = (index[1]-index[0])*0.8,yerr=asymmetric_error,capsize=2)
a[1][2].set_title('sim Jingzhaotoxin-XI average 300-500 ns')
a[1][2].set_xticks(index)
a[1][2].set_xticklabels(index,fontsize=8,rotation=90)
#a[1][1].set_ylabel('occurrences')
a[1][2].set_xlabel('residues')
a[1][2].legend(lines, ['polar', 'no polar', 'acid','basic'], loc='upper center',bbox_to_anchor=(0.22, 0.5, 0.5, 0.5),fontsize=10)




#plt.bar(index, result[0],width = (index[1]-index[0])*0.8,yerr=asymmetric_error,capsize=2)
#plt.xticks(index, index, fontsize=5, rotation=0)
#plt.title('sim Jingzhaotoxin-XI average 300-500 ns')
#plt.ylabel('relative ocurrences')
#plt.xlabel('residue #')
#plt.legend(lines, ['polar', 'no polar', 'acid','basic'], loc='upper center',bbox_to_anchor=(0.22, 0.5, 0.5, 0.5),fontsize=7)
#plt.savefig("/media/eniac/mdd1/paper_membranas/POPG/popg_residues_all_contact.png", dpi=900) 
plt.show()

#### colorear por beta factor
#pdb='/media/eniac/mdd1/paper_membranas/xi_jinghao_toxin/XI_ingaotoxin.pdb'
#output='/media/eniac/mdd1/paper_membranas/XI_ingaotoxin_colour_beta.pdb'
#molid = molecule.load("pdb", pdb)
#print(molid)
#query = "protein"
#vector=np.linspace(0, 1, num=len(protein.beta),)
#Trajectory.set_beta_factor(result[0],query,molid,output)
#color scale
# small=red, middle=white, large=blue


# %%
