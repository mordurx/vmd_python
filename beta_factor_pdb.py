from vmd import molecule,atomsel
import numpy as np
def set_beta_factor(residue_occupancy:list):
    for x in protein.resid:
        atomsel("protein and resid "+str(x)).beta=residue_occupancy[x-1]
        print(x,residue_occupancy[x-1])
pdb='poses/snx_chanel.pdb'
molid = molecule.load("pdb", pdb)
print(molid)

protein = atomsel("segname TOX")
vector=np.linspace(0, 1, num=len(protein.beta))
set_beta_factor(vector)

#protein.beta=vector
print (len(protein.beta)) 
protein.write("pdb","caca.pdb")
molecule.delete(molid)
