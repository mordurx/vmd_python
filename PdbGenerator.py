from trajectory import Trajectory
from vmd import *
import random
class PdbGenerator(Trajectory):
    def __init__(self,pdb,seletions,output,psf=None):
        super().__init__(pdb,psf=None)
        self.__atom_seletion=seletions
        self.__output=output
        
    @property
    def box_size(self):
        return self.__box_size    
    @property
    def atom_seletion(self):
        return self.__atom_seletion
    @property
    def output(self):
        return self.__output
    
    @atom_seletion.setter
    def atom_seletion(self, value):
        self.__atom_seletion = value    
    @output.setter
    def output(self, value):
        self.__output = value 
    def self_assembly(self,box_size):
        x0=box_size[0][0]
        y0=box_size[0][1]
        z0=box_size[0][2]
        
        x1=box_size[1][0]
        y1=box_size[1][1]
        z1=box_size[1][2]
        
        #genera posicion random para un grupo de atomos
        dict_rotations_allowed={}
        dict_rotations_allowed["0"]=(1.0, 0.0, 0.0, 0.0, 0.0, 1.0, -0.0, 0.0,0.0, 0.0, 1.0, 0.0,0.0, 0.0, 0.0, 1.0)
        dict_rotations_allowed["1"]=(1.0, 0.0, 0.0, 0.0, 0.0, 0.8660254037844387, -0.49999999999999994, 0.0   ,0.0, 0.49999999999999994, 0.8660254037844387, 0.0,0.0, 0.0, 0.0, 1.0)
        dict_rotations_allowed["2"]=(1.0, 0.0, 0.0, 0.0, 0.0, 0.5000000000000001, -0.8660254037844386, 0.0   ,0.0, 0.8660254037844386 ,0.5000000000000001, 0.0,0.0, 0.0, 0.0, 1.0)
        dict_rotations_allowed["3"]=(1.0, 0.0, 0.0, 0.0, 0.0, 6.123233995736766e-17, -1.06, 0.0   ,0.0, 1.0, 6.123233995736766e-17, 0.0,0.0, 0.0, 0.0, 1.0)
        dict_rotations_allowed["4"]=(1.0, 0.0, 0.0, 0.0, 0.0, -0.4999999999999998, -0.8660254037844387, 0.0   ,0.0, 0.8660254037844387, -0.4999999999999998, 0.0,0.0, 0.0, 0.0, 1.0)
        dict_rotations_allowed["5"]=(1.0, 0.0, 0.0, 0.0, 0.0, -0.8660254037844387,   -0.49999999999999994, 0.0   ,0.0, 0.49999999999999994, -0.8660254037844387, 0.0,0.0, 0.0, 0.0, 1.0)
        dict_rotations_allowed["5"]=(1.0, 0.0, 0.0, 0.0, 0.0, -1.0, -1.2246467991473532e-16, 0.0   ,0.0, 1.2246467991473532e-16, -1.0, 0.0,0.0, 0.0, 0.0, 1.0)
        dict_rotations_allowed["6"]=(1.0, 0.0, 0.0, 0.0, 0.0, -0.8660254037844386, 0.5000000000000001, 0.0   ,0.0, -0.5000000000000001, -0.8660254037844386, 0.0,0.0, 0.0, 0.0, 1.0)
        dict_rotations_allowed["7"]=(1.0, 0.0, 0.0, 0.0, 0.0, -1.8369701987210297e-16,1.0, 0.0   ,0.0, -1.0, -1.8369701987210297e-16, 0.0,0.0, 0.0, 0.0, 1.0)
        dict_rotations_allowed["8"]=(1.0, 0.0, 0.0, 0.0, 0.0, 0.5000000000000001, 0.8660254037844386, 0.0   ,0.0, -0.8660254037844386, 0.5000000000000001, 0.0,0.0, 0.0, 0.0, 1.0)
        dict_rotations_allowed["9"]=(1.0, 0.0, 0.0, 0.0, 0.0, 0.8660254037844384, 0.5000000000000004, 0.0   ,0.0, -0.5000000000000004, 0.8660254037844384, 0.0,0.0, 0.0, 0.0, 1.0)
        dict_rotations_allowed["10"]=(1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.4492935982947064e-16, 0.0   ,0.0, -2.4492935982947064e-16, 1.0, 0.0,0.0, 0.0, 0.0, 1.0)
        
        memb = atomsel(self.__atom_seletion)
        num_lipids=memb.resid[-1]
        #sel1 = atomsel(selection="resid 1", molid=self.molID)
        list_exclude=[] 
        for resid in range(1,num_lipids+1):
            sel1 = atomsel(selection="resid "+str(resid), molid=self.molID)
            print ("resid ",resid)
            
            x=random.choice([i for i in range(x0,x1) if i not in list_exclude])
            #list_exclude.append(x)
            y=random.choice([i for i in range(y0,y1) if i not in list_exclude])
            #list_exclude.append(y)
            z=random.choice([i for i in range(z0,z1) if i not in list_exclude])
            #list_exclude.append(z)
            keys=dict_rotations_allowed.keys()
            key_random=random.choice(list(keys))
            #print(key_random)
            sel1.move(dict_rotations_allowed[key_random])
            sel1.moveby((x,y,z))
        molecule.write(self.molID,"pdb",self.__output)
        print (self.get_pdb_path())
    
        


