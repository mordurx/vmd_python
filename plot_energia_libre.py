from trajectory import Trajectory
import matplotlib.pyplot as plt
import numpy as np
from vmd import molecule,atomsel
from matplotlib  import rcParams
    
# sim1_dcd='/mnt/e/paper_membranas/xi_jinghao_toxin/unwrapped_popc_popg_XI.dcd'
# sim1_psf='/mnt/e/paper_membranas/xi_jinghao_toxin/sim1_ix_1_popc_popg.Wat.psf'

# sim2_dcd='/mnt/e/paper_membranas/xi_jinghao_toxin/unwrapped_3popg_popg_XI.dcd'
# sim2_psf='/mnt/e/paper_membranas/xi_jinghao_toxin/ix_fix_ions.Wat.psf'

# sim3_dcd= "/mnt/e/paper_membranas/xi_jinghao_toxin/simXI_3popc_popg.dcd"
# sim3_psf="/mnt/e/paper_membranas/xi_jinghao_toxin/IX_final1.Wat.psf"

# sim4_dcd='/mnt/e/paper_membranas/xi_jinghao_toxin/unwrapped_popg_XI.dcd'
# sim4_psf='/mnt/e/paper_membranas/xi_jinghao_toxin/sim1_IX_popg_ions.Wat.psf'



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

r1_popg =traj1.porcentaje_contact_free_energy("same residue as (protein within 5 of name P)",1000,5000,[[9,10,11],[27,29,30]])
r2_popg =traj2.porcentaje_contact_free_energy("same residue as (protein within 5 of name P)",1000,5000,[[9,10,11],[27,29,30]])
r3_popg =traj3.porcentaje_contact_free_energy("same residue as (protein within 5 of name P)",1000,5000,[[9,10,11],[27,29,30]])
r4_popg =traj4.porcentaje_contact_free_energy("same residue as (protein within 5 of name P)",1000,5000,[[9,10,11],[27,29,30]])

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

r1_3popg =traj1_3popg.porcentaje_contact_free_energy("same residue as (protein within 5 of name P)",1000,5000,[[9,10,11],[27,29,30]])
r2_3popg =traj2_3popg.porcentaje_contact_free_energy("same residue as (protein within 5 of name P)",1000,5000,[[9,10,11],[27,29,30]])
r3_3popg =traj3_3popg.porcentaje_contact_free_energy("same residue as (protein within 5 of name P)",1000,5000,[[9,10,11],[27,29,30]])
r4_3popg =traj4_3popg.porcentaje_contact_free_energy("same residue as (protein within 5 of name P)",1000,5000,[[9,10,11],[27,29,30]])

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
r1_popc =traj1_popc.porcentaje_contact_free_energy("same residue as (protein within 5 of name P)",3000,7000,[[9,10,11],[27,29,30]])
r2_popc =traj2_popc.porcentaje_contact_free_energy("same residue as (protein within 5 of name P)",3000,7000,[[9,10,11],[27,29,30]])

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

r1_popc_popg =traj1_popc_popg.porcentaje_contact_free_energy("same residue as (protein within 5 of name P)",1000,5000,[[9,10,11],[27,29,30]])
r2_popc_popg =traj2_popc_popg.porcentaje_contact_free_energy("same residue as (protein within 5 of name P)",1000,5000,[[9,10,11],[27,29,30]])
r3_popc_popg =traj3_popc_popg.porcentaje_contact_free_energy("same residue as (protein within 5 of name P)",1000,5000,[[9,10,11],[27,29,30]])
r4_popc_popg =traj4_popc_popg.porcentaje_contact_free_energy("same residue as (protein within 5 of name P)",1000,5000,[[9,10,11],[27,29,30]])

traj1_popc_popg.close()
traj2_popc_popg.close()
traj3_popc_popg.close()
traj4_popc_popg.close()

"""3POPC:1POPG"""
sim1_3popc_popg_dcd= "/mnt/e/paper_membranas/3POPC_popg_resume/sim1_3popc_popg_unwrap.dcd"
sim1_3popc_popg_psf="/mnt/e/paper_membranas/3POPC_popg_resume/sim1_3popc_popg_unwrap.psf"

sim2_3popc_popg_dcd= "/mnt/e/paper_membranas/3POPC_popg_resume/sim2_3popc_popg_unwrap.dcd"
sim2_3popc_popg_psf="/mnt/e/paper_membranas/3POPC_popg_resume/sim2_3popc_popg_unwrap.psf"

sim3_3popc_popg_dcd= "/mnt/e/paper_membranas/3POPC_popg_resume/sim3_3popc_popg_unwrap.dcd"
sim3_3popc_popg_psf="/mnt/e/paper_membranas/3POPC_popg_resume/sim3_3popc_popg_unwrap.psf"

sim4_3popc_popg_dcd= "/mnt/e/paper_membranas/3POPC_popg_resume/sim4_3popc_popg_filt_unwrap.dcd.dcd"
sim4_3popc_popg_psf="/mnt/e/paper_membranas/3POPC_popg_resume/sim4_3popc_popg_filt_unwrap.dcd.psf"

traj1_3popc_popg=Trajectory(sim1_3popc_popg_dcd,sim1_3popc_popg_psf,stride=1)
traj2_3popc_popg=Trajectory(sim2_3popc_popg_dcd,sim2_3popc_popg_psf,stride=1)
traj3_3popc_popg=Trajectory(sim3_3popc_popg_dcd,sim3_3popc_popg_psf,stride=1)
traj4_3popc_popg=Trajectory(sim4_3popc_popg_dcd,sim4_3popc_popg_psf,stride=1)

r1_3popc_popg =traj1_3popc_popg.porcentaje_contact_free_energy("same residue as (protein within 5 of name P)",1000,5000,[[9,10,11],[27,29,30]])
r2_3popc_popg =traj2_3popc_popg.porcentaje_contact_free_energy("same residue as (protein within 5 of name P)",1000,5000,[[9,10,11],[27,29,30]])
r3_3popc_popg =traj3_3popc_popg.porcentaje_contact_free_energy("same residue as (protein within 5 of name P)",1000,5000,[[9,10,11],[27,29,30]])
r4_3popc_popg =traj4_3popc_popg.porcentaje_contact_free_energy("same residue as (protein within 5 of name P)",1000,5000,[[9,10,11],[27,29,30]])

traj1_3popc_popg.close()
traj2_3popc_popg.close()
traj3_3popc_popg.close()
traj4_3popc_popg.close()

"""IX toxin"""
sim1_IX_popg_dcd= "/mnt/e/paper_membranas/xi_jinghao_toxin/IX_TOXIN_resume/unwrap_XI_sim1_POPG.dcd"
sim1_XI_popg_psf="/mnt/e/paper_membranas/xi_jinghao_toxin/IX_TOXIN_resume/unwrap_XI_sim1_POPG.psf"

sim2_3popg_popc_dcd= "/mnt/e/paper_membranas/xi_jinghao_toxin/IX_TOXIN_resume/unwrap_3popg_popg_XI.dcd"
sim2_3popg_popc_psf="/mnt/e/paper_membranas/xi_jinghao_toxin/IX_TOXIN_resume/unwrap_3popg_popg_XI.psf"

sim3_XI_popc_popg_dcd= "/mnt/e/paper_membranas/xi_jinghao_toxin/IX_TOXIN_resume/unwrap_POPC_POPG_sim1.dcd"
sim3_XI_popc_popg_psf="/mnt/e/paper_membranas/xi_jinghao_toxin/IX_TOXIN_resume/unwrap_POPC_POPG_sim1.psf"

sim4_3popc_popg_dcd= "/mnt/e/paper_membranas/xi_jinghao_toxin/IX_TOXIN_resume/unwrap_XI_3POPC_POPG.dcd"
sim4_3popc_popg_psf="/mnt/e/paper_membranas/xi_jinghao_toxin/IX_TOXIN_resume/unwrap_XI_3POPC_POPG.psf"

sim5_popc_dcd= "/mnt/e/paper_membranas/xi_jinghao_toxin/IX_TOXIN_resume/unwrap_XI_sim3_POPC.dcd"
sim5_popc_psf="/mnt/e/paper_membranas/xi_jinghao_toxin/IX_TOXIN_resume/unwrap_XI_sim3_POPC.psf"

traj1_XI_popg=Trajectory(sim1_IX_popg_dcd,sim1_XI_popg_psf,stride=1)
traj2_XI_3popg_popc=Trajectory(sim2_3popg_popc_dcd,sim2_3popg_popc_psf,stride=1)
traj3_XI_popc_popg=Trajectory(sim3_XI_popc_popg_dcd,sim3_XI_popc_popg_psf,stride=1)
traj4_XI_3popc_popg=Trajectory(sim4_3popc_popg_dcd,sim4_3popc_popg_psf,stride=1)
traj5_XI_popc=Trajectory(sim5_popc_dcd,sim5_popc_psf,stride=1)

r1_XI_popg =traj1_XI_popg.porcentaje_contact_free_energy("same residue as (protein within 5 of name P)",1000,5000,[[9,10,11],[27,29,30]])
r2_XI_3popg_popc =traj2_XI_3popg_popc.porcentaje_contact_free_energy("same residue as (protein within 5 of name P)",1000,5000,[[9,10,11],[27,29,30]])
r3_XI_popc_popg =traj3_XI_popc_popg.porcentaje_contact_free_energy("same residue as (protein within 5 of name P)",1000,5000,[[9,10,11],[27,29,30]])
r4_XI_3popc_popg =traj4_XI_3popc_popg.porcentaje_contact_free_energy("same residue as (protein within 5 of name P)",1000,5000,[[9,10,11],[27,29,30]])
r5_XI_popc =traj5_XI_popc.porcentaje_contact_free_energy("same residue as (protein within 5 of name P)",1000,5000,[[9,10,11],[27,29,30]])

traj1_XI_popg.close()
traj2_XI_3popg_popc.close()
traj3_XI_popc_popg.close()
traj4_XI_3popc_popg.close()
traj5_XI_popc.close()


# importing xlwt module 
import xlwt 

R=1.9872036e-3 #kcal/mol/k
T=298    #kelvin



workbook = xlwt.Workbook()  
sheet = workbook.add_sheet("sim pOPG") 
# Specifying style 
style = xlwt.easyxf('font: bold 1') 
# Specifying column 
sheet.write(0, 0, 'POPG', style) 

sheet.write(0, 1, 'unbound', style) 
sheet.write(0, 2, 'bound', style) 
sheet.write(0, 3, '-RT ln(bound/unbound)', style) 

sheet.write(1, 0, 'sim1', style) 
sheet.write(1, 1, r1_popg[0], style) 
sheet.write(1, 2, r1_popg[1], style) 
sheet.write(1, 3, -R*T * np.log(r1_popg[1]/r1_popg[0]), style) 


sheet.write(2, 0, 'sim2', style) 
sheet.write(2, 1, r2_popg[0], style) 
sheet.write(2, 2, r2_popg[1], style) 
sheet.write(2, 3, -R*T * np.log(r2_popg[1]/r2_popg[0]), style) 



sheet.write(3, 0, 'sim3', style) 
sheet.write(3, 1, r3_popg[0], style) 
sheet.write(3, 2, r3_popg[1], style) 
sheet.write(3, 3, -R*T * np.log(r3_popg[1]/r3_popg[0]), style) 


sheet.write(4, 0, 'sim4', style) 
sheet.write(4, 1, r4_popg[0], style) 
sheet.write(4, 2, r4_popg[1], style) 
"""3POPG"""
# Specifying column 
sheet.write(8, 0, '3POPG_POPC', style) 
sheet.write(8, 1, 'unbound', style) 
sheet.write(8, 2, 'bound', style)

sheet.write(9, 0, 'sim1', style) 
sheet.write(9, 1, r1_3popg[0], style) 
sheet.write(9, 2, r1_3popg[1], style)

sheet.write(10, 0, 'sim2', style) 
sheet.write(10, 1, r2_3popg[0], style) 
sheet.write(10, 2, r2_3popg[1], style)

sheet.write(11, 0, 'sim3', style) 
sheet.write(11, 1, r3_3popg[0], style) 
sheet.write(11, 2, r3_3popg[1], style)

sheet.write(12, 0, 'sim4', style) 
sheet.write(12, 1, r4_3popg[0], style) 
sheet.write(12, 2, r4_3popg[1], style)

sheet.write(16, 0, 'popc', style) 
sheet.write(16, 1, 'unbound', style) 
sheet.write(16, 2, 'bound', style) 

sheet.write(17, 0, 'sim1', style) 
sheet.write(17, 1, r1_popc[0], style) 
sheet.write(17, 2, r1_popc[1], style) 
sheet.write(18, 0, 'sim2', style) 
sheet.write(18, 1, r2_popc[0], style) 
sheet.write(18, 2, r2_popc[1], style)

sheet.write(20, 0, 'POPC_POPG', style) 
sheet.write(20, 1, 'unbound', style) 
sheet.write(20, 2, 'bound', style) 


sheet.write(21, 0, 'sim1', style) 
sheet.write(21, 1, r1_popc_popg[0], style) 
sheet.write(21, 2, r1_popc_popg[1], style) 
sheet.write(22, 0, 'sim2', style) 
sheet.write(22, 1, r2_popc_popg[0], style) 
sheet.write(22, 2, r2_popc_popg[1], style) 
sheet.write(23, 0, 'sim3', style) 
sheet.write(23, 1, r3_popc_popg[0], style) 
sheet.write(23, 2, r3_popc_popg[1], style) 
sheet.write(24, 0, 'sim4', style) 
sheet.write(24, 1, r4_popc_popg[0], style) 
sheet.write(24, 2, r4_popc_popg[1], style) 

sheet.write(27, 0, '3popc_popg', style) 
sheet.write(27, 1, 'unbound', style) 
sheet.write(27, 2, 'bound', style) 
sheet.write(28, 0, 'sim1', style) 
sheet.write(28, 1, r1_3popc_popg[0], style) 
sheet.write(28, 2, r1_3popc_popg[1], style) 
sheet.write(29, 0, 'sim2', style) 
sheet.write(29, 1, r2_3popc_popg[0], style) 
sheet.write(29, 2, r2_3popc_popg[1], style) 
sheet.write(30, 0, 'sim3', style) 
sheet.write(30, 1, r3_3popc_popg[0], style) 
sheet.write(30, 2, r3_3popc_popg[1], style) 
sheet.write(31, 0, 'sim4', style) 
sheet.write(31, 1, r4_3popc_popg[0], style) 
sheet.write(31, 2, r4_3popc_popg[1], style) 

sheet.write(33, 0, 'XI_toxin', style) 
sheet.write(33, 1, 'unbound', style) 
sheet.write(33, 2, 'bound', style) 

sheet.write(34, 0, 'POPG', style) 
sheet.write(34, 1, r1_XI_popg[0], style) 
sheet.write(34, 2, r1_XI_popg[1], style) 

sheet.write(35, 0, '3POPG_POPC', style) 
sheet.write(35, 1, r2_XI_3popg_popc[0], style) 
sheet.write(35, 2, r2_XI_3popg_popc[1], style) 

sheet.write(36, 0, 'POPG_POPC', style) 
sheet.write(36, 1, r3_XI_popc_popg[0], style) 
sheet.write(36, 2, r3_XI_popc_popg[1], style) 

sheet.write(37, 0, '3POPC_POPG', style) 
sheet.write(37, 1, r4_XI_3popc_popg[0], style) 
sheet.write(37, 2, r4_XI_3popc_popg[1], style) 

sheet.write(38, 0, 'POPC', style) 
sheet.write(38, 1, r5_XI_popc[0], style) 
sheet.write(38, 2, r5_XI_popc[1], style) 


workbook.save("sim_POPG_all.xls") 

