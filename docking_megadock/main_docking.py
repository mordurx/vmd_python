#!/usr/bin/env python3
#%%%
from blockeo import *
import numpy as np
import pandas as pd
import subprocess
from subprocess import Popen, PIPE
import matplotlib.pyplot as plt
sys.path.insert(0, '/home/mordurx/vmd_python')
import os

from trajectory import Trajectory
#resid 844 to 872  resid 902 to 912 and protein resid 794 to 808 and protein
ban='844-872,902-912,794-808,853A,853B,853C,853D,853E'
#pdb ='docking_megadock/vsd3_frame10000.pdb'
#receptor ='docking_megadock/vsd3_frame10000_ban.pdb'
ligand='docking_megadock/snx_pose4.pdb'
#ban='100-110'

## vsdca1.1 embrionary
sim1_dcd= "/home/mordurx/cav23/vsd3_resting_2.5us_align.dcd"
sim1_psf="/home/mordurx/cav23/vsd3_last_frame.Wat.psf"
#sim1_dcd= "/home/mordurx/cav23/vsd4_resting_mem_tox_align.dcd"
#sim1_psf="/home/mordurx/cav23/vsd4_resting_mem_tox_align.psf"


traj1=Trajectory(sim1_dcd,sim1_psf,stride=10000)
print (traj1.num_frames())
score_vec=[]
time=[]
for frame in range(traj1.num_frames()):
    new_receptor='docking_megadock/pdb_'+str(frame)+'.pdb'
    pdb=traj1.get_pdb_trajectory("output",frame)
    
    Bloqueo.block_resid(pdb,ban,new_receptor)
    megadock=subprocess.call("megadock-gpu -R "+new_receptor+" -L "+ ligand +" -o guido.out  -O ",shell=True)
    decoy=subprocess.call("decoygen docking_megadock/v.pdb "+ ligand+" guido.out 1",shell=True)
    cols=["total","ELEC "]
    df = pd.read_csv("guido.csv",usecols=cols)
    df.sort_values(by=['ELEC '], inplace=True)
    df.head()
    print (df['ELEC '][0])
    score_vec.append( df['ELEC '][0])
    time.append(frame)
    os.remove(pdb)
    os.remove(new_receptor)
    
df = pd.DataFrame([time,score_vec])
df=df.T
writer = pd.ExcelWriter('test.xlsx', engine='xlsxwriter')
df.to_excel(writer, sheet_name='welcome', index=False)
writer.save()

plt.plot(time,score_vec )
plt.show()
traj1.close()
plt.savefig("vsd3_score.png", dpi=900) 







#args2 = ['/home/eniac/Documentos/megadock-4.1.1/megadock-gpu','-R','VSD4_rep.cent_final_fix.pdb','-L','snx.notWater.pdb']
#p=subprocess.Popen(args2,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True) 
#p=subprocess.Popen(args2,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True) 

# %%
