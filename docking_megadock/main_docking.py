#!/usr/bin/env python3
from blockeo import *
import numpy as np
import pandas as pd
#resid 844 to 872  resid 902 to 912 and protein resid 794 to 808 and protein
ban='844-872,902-912,794-808,853A,853B,853C,853D,853E'
pdb ='docking_megadock/vsd3_frame10000.pdb'
receptor ='docking_megadock/vsd3_frame10000_ban.pdb'
ligand='docking_megadock/snx_pose4.pdb'
Bloqueo.block_resid(pdb,ban,receptor)

import subprocess
from subprocess import Popen, PIPE
args2 = ['/home/eniac/Documentos/megadock-4.1.1/megadock-gpu','-R','VSD4_rep.cent_final_fix.pdb','-L','snx.notWater.pdb']
#p=subprocess.Popen(args2,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True) 
#p=subprocess.Popen(args2,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True) 
megadock=subprocess.call("megadock-gpu -R "+receptor+" -L "+ ligand +" -o guido.out  -O ",shell=True)
decoy=subprocess.call("decoygen docking_megadock/v.pdb "+ ligand+" guido.out 1",shell=True)

cols=["total"]
df = pd.read_csv("guido.csv",usecols=cols)
df.head()
print (df['total'][0])