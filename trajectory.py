

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 15:33:39 2019

@author: eniac
"""
import collections
import random
from numpy import savetxt
import subprocess
from vmd import molecule,atomsel,vmdnumpy,animate,trans
import numpy as np
import pandas as pd
import statistics
from sklearn import preprocessing
import vmd
from sklearn.linear_model import LinearRegression
class Trajectory:
    def __init__(self,dcd=None,psf=None,first=0,last=-1,stride=1,waitfor=-1,pdb=None,type_format="NAMD"):
        
        if type_format=="GROMACS":
            traj_format="xtc"
            top_format="gro"
        if type_format=="NAMD":
            traj_format="dcd"
            top_format="psf"
        if type_format=="test":
            traj_format="dcd"
            top_format="gro"
        self.molID=0
        if psf is not None:
                self.psf = psf
                self.molID=molecule.load(top_format,psf)
        if  (pdb is None) and (dcd is not None):
            try:
                self.dcd = dcd
                self.molID=molecule.new('new')
                self.molID=molecule.load(top_format,self.psf) # load trajectory
                molecule.read(self.molID,traj_format,self.dcd,stride=stride,first=first,last=last,waitfor=waitfor) 
                print (traj_format+" file detected id ",self.molID)
            except IOError:
                print ("Could not read dcd file or psf:", dcd)
                raise Exception()
        else:
            try:
                self.molID=molecule.new('pdb')
                self.pdb=pdb
                molecule.read(molid=self.molID,filetype ='pdb',filename=self.pdb,first=0,last=0) 
                print ("pdb file detected id ",self.molID)
            except IOError:
                print ("Could not read dcd file or psf:",pdb)
                raise Exception()
    @classmethod
    def dcd(self,dcd,psf,first=0,last=-1,stride=1,waitfor=-1):
        self.psf = psf
        self.dcd = dcd
        self.molID=0
        try:
        
            self.molID=molecule.new('new')
            self.molID=molecule.load('psf',self.psf) # load trajectory
            molecule.read(self.molID,'dcd',self.dcd,stride=stride,first=first,last=last,waitfor=waitfor) 
            print (self.molID)
        except IOError:
            print ("Could not read dcd file or psf:", dcd)
            raise Exception()
    def pdb(self,pdb,psf=None):
        #sobrecarga para leer pdbs y no dcd
        try:
            self.molID=molecule.new('pdb')
            self.pdb=pdb
            if psf is not None:
                self.molID=molecule.load('psf',psf)
                self.psf=psf # add psf
            molecule.read(molid=self.molID,filetype ='pdb',filename=self.pdb,first=0,last=0) 
            print (self.molID)
        except IOError:
            print ("Could not read dcd file or psf:",pdb)
            raise Exception()
    
    def get_pdb_path(self):
        return self.pdb
    def tcl_dipoleZ(self,file_path,proc):
        #calcula el dipolo a partir de tcl file.
        Cosangle=[]
        magnitude_prot_vector=[]
        for frame in range(Trajectory.num_frames(self)):
            #seleciono la cabezas
            sel1 = atomsel(selection="protein", molid=self.molID, frame=frame)
            prot_z=sel1.center(sel1.mass)[2]
            vmd.evaltcl("source "+file_path)
            result=vmd.evaltcl(proc+" "+str(self.molID)+" "+str(frame))
            dip_prot=np.fromstring(result, dtype=float, sep=' ')
            if prot_z >0:
                axis_z=[0,0,1]
                #print ("arriba",frame)
                
                
            else:
                axis_z=[0,0,-1]
                #print ("abajo")
            #magnitude vector dipole
            magnitude_dip = np.linalg.norm(dip_prot)
            magnitude_z_axis = np.linalg.norm(axis_z)
            #print(magnitude_dip)
            #print(magnitude_z_axis)
            
            #prod punto entre dipolo protein y eje z
            prod_dot=np.dot(dip_prot,axis_z)
            cos_teta=prod_dot/(magnitude_dip*magnitude_z_axis)
            Cosangle.append(cos_teta)
            magnitude_prot_vector.append(magnitude_dip)
            #print(np.fromstring(result, dtype=float, sep=' '))
            
        return magnitude_prot_vector,Cosangle,np.degrees(np.arccos(Cosangle)).tolist()
        
        
    def mean_displacement(self,atomselect1):
       
        displacement_x=[]
        displacement_y=[]
        displacement_z=[]
        
        #sampleo=int(np.floor ((Trajectory.num_frames(self)-1)/2))
        promedio_x=[]
        promedio_y=[]
        promedio_z=[]
        sel0 = atomsel(selection=atomselect1, molid=self.molID, frame=0) 
        count=0   
        for frame in range(0, 1000,1):
            print (frame)
            for frame1 in range( Trajectory.num_frames(self)-5000, Trajectory.num_frames(self)-1):
                sel0 = atomsel(selection=atomselect1, molid=self.molID, frame=frame1-1) 
               
                sel1 = atomsel(selection=atomselect1, molid=self.molID, frame=frame1+frame) 
                sel1_x = np.array((sel1.center(sel1.mass)[0]))
                sel0_x = np.array((sel0.center(sel0.mass)[0]))
                
                sel1_y = np.array((sel1.center(sel1.mass)[1]))
                sel0_y = np.array((sel0.center(sel0.mass)[1]))
                
                
                sel1_z = np.array((sel1.center(sel1.mass)[2]))
                sel0_z = np.array((sel0.center(sel0.mass)[2]))
                
                
                
                t0_index=frame1-1
                t1=frame1+frame
                if t1>= Trajectory.num_frames(self)-1:
                    count=count+1
                    displacement_x.append((sel1_x-sel0_x)**2)
                    displacement_y.append((sel1_y-sel0_y)**2)
                    displacement_z.append((sel1_z-sel0_z)**2)
                    
                    #print (frame1-1,frame1+frame) 
                    break
                #print (-1+frame1 ,frame1+frame)
                count=count+1
                #print (count)
                #protein  = atomsel(selection="protein", molid=molid, frame=frame)
               
                #print (np.power(sel1_x-sel0_x, 2))
                #print ((sel1_x-sel0_x)**2)
                displacement_x.append((sel1_x-sel0_x)**2)
                displacement_y.append((sel1_y-sel0_y)**2)
                displacement_z.append((sel1_z-sel0_z)**2)
                
                #com_molecule.append(np.power(sel1_x, 2))
                #sel0_x = sel1_x
                #print (frame1-1,frame1+frame) 
            promedio_x.append(np.sum(displacement_x)/count)
            promedio_y.append(np.sum(displacement_y)/count)
            promedio_z.append(np.sum(displacement_z)/count)
            
            #print (count)
            count=0  
            displacement_x=[]
            displacement_y=[]
            displacement_z=[]
        return [promedio_x,promedio_y,promedio_z]
    @staticmethod
    def the_static_method(vector_pdb,output,outputDCD):
       files=' '.join(map(str,vector_pdb))
       #print (files)
       i=0
       with open('sort_ligand_order.txt', 'w') as filehandle:
           for listitem in vector_pdb:
               #print (listitem)
               nameFile=listitem.split("/")
               filehandle.write('%s' % str(i)+"\t"+nameFile[-1]+"\n")
               i+=1
       
       
       subprocess.call("cat "+files+"> "+ output,shell=True)
       molID=molecule.load2('pdb',output)
       print (molecule.numframes(molID))
       sel1 = atomsel(selection="segname TOX") 
       molecule.write(molID,'dcd',outputDCD,0,-1,1,sel1)
       molecule.cancel
       return output
        

        
    def close (self):
        print("cerrando molecule "+str(self.molID))
        molecule.cancel(self.molID)
        molecule.delete(self.molID)
         
    
    def num_frames(self):
        return molecule.numframes(self.molID)
    
       
    def get_molid(self):
        return self.molID
    def minmax(self,atomselect1):
        list_x=[]
        list_y=[]
        list_z=[]
        for frame in range(Trajectory.num_frames(self)):
            #seleciono la cabezas
            sel1 = atomsel(selection=atomselect1, molid=self.molID, frame=frame)
            sel_x= list_x.append(np.abs(sel1.minmax()[0][0]-sel1.minmax()[1][0]))
            sel_y= list_y.append(np.abs(sel1.minmax()[0][1]-sel1.minmax()[1][1]))
            sel_z= list_z.append(np.abs(sel1.minmax()[0][2]-sel1.minmax()[1][2]))
            sel_down= sel1.minmax()[1]
            
            
        #x = statistics.mean(data)
        #y = statistics.mean(data)
        #z = statistics.mean(data)
        return [list_x,list_y,list_z]   
    def get_membrane_outer_distance(self,sel_name_P,membrane):
        pos_z=[]
        #obtiene un promedio la las cabezas fosfolipidas en z por cada frame
        for frame in range(Trajectory.num_frames(self)):
            #seleciono la cabezas
            sel1 = atomsel(selection=sel_name_P, molid=self.molID, frame=frame)
            sel= sel1.centerperresidue(sel1.mass)
            name_p_z=(list(zip(*sel))[2])
            mean=np.mean(name_p_z)
            #membrane
            selM= atomsel(selection=membrane, molid=self.molID, frame=frame)
            sel= selM.center(selM.mass)[2]
            
            pos_z.append(np.abs(mean-sel))
        return pos_z    

    def delete_frames(self,atomselect1,cutoff,new_dcd_path,stride):
        #hay que limpiar y depurar bien 12-03-21
        distance_mass_weight=[]
        count=0
        for frame in range(Trajectory.num_frames(self)):
            #protein  = atomsel(selection="protein", molid=molid, frame=frame) 
            sel1 = atomsel(selection=atomselect1, molid=self.molID, frame=frame) 
            #sel2 =atomsel(selection=atomselect2, molid=self.molID, frame=frame)
            
            sel1_x = np.array((sel1.center(sel1.mass)[0]))
            sel2_y= np.array((sel1.center(sel1.mass)[1]))
            if sel1_x < 0:
                print ("borrando "+str(frame))
                count=count+1
                molecule.delframe(self.molID,frame,frame)
                continue
            distance_mass_weight.append(sel1_x)

            #sel1.write('dcd','vsd1000.dcd')
            #
            #
        print ("borrados ",count)    
        molecule.write(self.molID,"dcd",new_dcd_path,stride)
        
            
        return distance_mass_weight
    
    def center_simulations(self,atomselect1,new_dcd_path,stride):
        #centra el dcd. similiar a vecinvert pero para trayectoria.
        count=0
        frame=1
        for frame in range(Trajectory.num_frames(self)):
            print(frame)
            sel_all= atomsel(selection="all", molid=self.molID, frame=frame) 
            sel1 = atomsel(selection=atomselect1, molid=self.molID, frame=frame) 
            sel_all.moveby(-1*np.array((sel1.center(sel1.mass))))
        molecule.write(self.molID,"dcd",new_dcd_path,stride)
    @staticmethod
    def to_excel(output,sheet,mode,**kwargs):
        columns = ["frame","event","contact","occupancy %"]
        df = pd.DataFrame(columns = columns)
        df["frame"] = pd.Series(kwargs["numHbound"].keys())
        df["event"] = pd.Series(kwargs["numHbound"].values())
        df["contact"] = pd.Series(kwargs["resume"].keys())
        df["occupancy %"] = pd.Series(kwargs["resume"].values())
        
        if not df.empty:
            df.columns=columns
        else:
            pd.DataFrame(columns=columns)
        writer = pd.ExcelWriter(output, engine='openpyxl',mode=mode)
        df.to_excel(writer, sheet_name=sheet, index=False)
        writer.sheets[sheet].column_dimensions['A'].width = 10
        writer.sheets[sheet].column_dimensions['B'].width = 10
        writer.sheets[sheet].column_dimensions['C'].width = 30
        writer.sheets[sheet].column_dimensions['D'].width = 30
        writer.save()
    @staticmethod
    def delete_frame_wrap(vector_distancia,cutoff):
        "aun no completoooo no usarr 5 sep 2020"
        valor_anterior=vector_distancia[0]
        cont=0
        vector_new_distance=[]
        for i in vector_distancia:
            if np.abs(valor_anterior-j)>=cutoff:
                    del vector_distancia[cont]
                    print ('deleting.. frame',cont)        
            else:
                vector_new_distance.append(i)
                valor_anterior=j
            cont=cont+1
        return vector_new_distance


    def porcentaje_contact_radio_gyrations(self,atomselect1,atomselect2):
        for frame in range(Trajectory.num_frames(self)):
            #protein  = atomsel(selection="protein", molid=molid, frame=frame) 
            sel1 = atomsel(selection=atomselect1, molid=self.molID, frame=frame) 
            sel2 =atomsel(selection=atomselect2, molid=self.molID, frame=frame)
            sel1_z = np.array((sel1.center(sel1.mass)[2]))
            sel2_z= np.array((sel2.center(sel2.mass)[2]))
            distance_mass_weight.append(np.linalg.norm(sel2_z-sel1_z))
            #(x−a)2 + (y−b)2 = r2
    
    
    def average_radius_of_gyration(self,atomselect1):
        vector_radio=[]
        for frame in range(Trajectory.num_frames(self)):
            
            protein  = atomsel(selection=atomselect1, molid=self.molID, frame=frame) 
            radio_giro=protein.rgyr(protein.mass)
            vector_radio.append(radio_giro)
            #sel1 = atomsel(selection=atomselect1, molid=self.molID, frame=frame) 
            #sel2 =atomsel(selection=atomselect2, molid=self.molID, frame=frame)
            #sel1_z = np.array((sel1.center(sel1.mass)[2]))
            #sel2_z= np.array((sel2.center(sel2.mass)[2]))
            #distance_mass_weight.append(np.linalg.norm(sel2_z-sel1_z))
         
        return vector_radio  
    def center_mass(self,atomselect1,file=None):
        
        center_mass=[]
        for frame in range(Trajectory.num_frames(self)):
            #protein  = atomsel(selection="protein", molid=molid, frame=frame) 
        
            sel = atomsel(selection=atomselect1, molid=Trajectory.molID, frame=frame) 
            com_mass = np.array(sel.center(sel.mass))/10
            center_mass.append(com_mass)
        if file is not None:
            with open(file, 'wb') as file1:
                np.savetxt(file1, center_mass, delimiter=' ',fmt='%.6f')     
        return center_mass  
    
        return membrane_center_mass     
    def distance_center_mass_Z(self,atomselect1):
        membrane_center_mass=[]
        for frame in range(Trajectory.num_frames(self)):
            #protein  = atomsel(selection="protein", molid=molid, frame=frame) 
        
            sel1 = atomsel(selection=atomselect1, molid=self.molID, frame=frame)
           
            #sel3 = atomsel(selection=atomselect2, molid=self.molID, frame=frame)  
            sel1_z = np.array((sel1.center(sel1.mass)[2]))
            membrane_center_mass.append(sel1_z)
            
        return membrane_center_mass
    @staticmethod    
    def to_csv(df,mode="w",sep=",",output="",header=True,index=False):
        #las lista deben tener el mismo size
        df.to_csv(mode=mode,path_or_buf=output, header=header, index=index,sep=sep)    
                 
    def get_pdb_trajectory(self,output,first_frame,last_frame):
        
        path_output=output+str(first_frame)+".pdb"
        #protein  = atomsel(selection="protein", molid=molid, frame=frame) 
        molecule.write(self.molID, "pdb",path_output, first=first_frame,last=last_frame)
        return path_output
            
        
    def distance_center_mass(self,atomselect1,atomselect2):
        distance_mass_weight=[]
        for frame in range(Trajectory.num_frames(self)):
            #protein  = atomsel(selection="protein", molid=molid, frame=frame) 
        
            sel1 =atomsel(selection=atomselect1, molid=self.molID, frame=frame) 
            sel2 =atomsel(selection=atomselect2, molid=self.molID, frame=frame)
            sel1_z = np.array((sel1.center(sel1.mass)[2]))
            sel2_z= np.array((sel2.center(sel2.mass)[2]))
            
            sel1_x = np.array((sel1.center(sel1.mass)[0]))
            sel2_x= np.array((sel2.center(sel2.mass)[0]))
            
            sel1_y = np.array((sel1.center(sel1.mass)[1]))
            sel2_y= np.array((sel2.center(sel2.mass)[1]))
            
            
            dist = np.sqrt((sel2_x - sel1_x)**2 + (sel2_y - sel1_y)**2+ (sel2_z - sel1_z)**2)  
            distance_mass_weight.append(dist)
    
        return distance_mass_weight
    def distance_center_mass_Z(self,atomselect1,atomselect2):
        distance_mass_weight=[]
        for frame in range(Trajectory.num_frames(self)):
            #protein  = atomsel(selection="protein", molid=molid, frame=frame) 
        
            sel1 =atomsel(selection=atomselect1, molid=self.molID, frame=frame) 
            sel2 =atomsel(selection=atomselect2, molid=self.molID, frame=frame)
            sel1_z = np.array((sel1.center(sel1.mass)[2]))
            sel2_z= np.array((sel2.center(sel2.mass)[2]))
            
            #sel1_x = np.array((sel1.center(sel1.mass)[0]))
            #sel2_x= np.array((sel2.center(sel2.mass)[0]))
            
            #sel1_y = np.array((sel1.center(sel1.mass)[1]))
            #sel2_y= np.array((sel2.center(sel2.mass)[1]))
            
            
            dist = np.sqrt((sel2_z - sel1_z)**2)  
            distance_mass_weight.append(dist)
    
        return distance_mass_weight
        
    def porcentaje_contact(self,atomselect1,atomselect2):
        protein = atomsel(selection=atomselect1, molid=self.molID, frame=0) 
        #reference1 = atomsel(selection="protein", molid=self.molID, frame=0) 
             #num residuesatomselect1
        #num_residues=len(pd.factorize(protein.resid)[1])
        ocurrencias_vector=[]
        ocurrencias_temp=[]
              
        for frame in range(Trajectory.num_frames(self)):
             protein  = atomsel(selection="protein", molid=self.molID, frame=frame) 
             # use frame 0 for the reference
             #sel1 = atomsel(selection=atomselect1, molid=self.molID, frame=frame) 
             sel2 =atomsel(selection=atomselect2, molid=self.molID, frame=frame)
             resid_center=protein.centerperresidue()
             if frame==0:
                 ocurrencias_vector= np.zeros(len(resid_center))
                 ocurrencias_temp=np.zeros(len(resid_center))
        
        
            
             resid_z=list(map(lambda x: x[2], resid_center))
             nameP_z=list(map(lambda x: x[2], sel2.minmax()))
             for number_resid in range(len(resid_z)):
                 if resid_z[number_resid]>0.8*nameP_z[0] and resid_z[number_resid]<nameP_z[1]*0.8:
                     #print ("choca"+str(number_resid)+' '+str(frame))
                     ocurrencias_temp[number_resid]=1 
                 else:
                     ocurrencias_temp[number_resid]=0
                     #print ("no choca"+str(number_resid)+' '+str(frame))
             ocurrencias_vector+=ocurrencias_temp
             #print (ocurrencias_temp)
             ocurrencias_temp= np.zeros(len(resid_center))    
             #(>) y “menor que” (<). C
             #    sel2 =atomsel(selection=atomselect2, molid=self.molID, frame=frame)
             #sel1_z = np.array((sel1.center(sel1.mass)[2]))
             #sel2_z= np.array((sel2.center(sel2.mass)[2]))
             #distance_mass_weight.append(np.linalg.norm(sel2_z-sel1_z))
    
        return ocurrencias_vector/Trajectory.num_frames(self)
    def contact_time(self,atomselect1,atomselect2,cutoff):
        dict_resid_count={}
        # cuantos contactos tiene un ligando a cierto receptor by cutoff
        
        for frame in range(Trajectory.num_frames(self)):
             sel  = atomsel(selection=atomselect1, molid=self.molID, frame=frame)
             receptor  = atomsel(selection=atomselect2, molid=self.molID, frame=frame)
             contact=sel.contacts(receptor,cutoff)
             #la primera lista index del lig
             #segunda lista atom que hacen contacto
             aminoacid = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M','HSE': 'H'}
             #dejar solo valores unicos, no nos intereza el analisis de enlaces
             if len(contact[1])==0:
                 print("vacioo")
                 continue
             unique_contact=list(set(contact[1]))
             unique_index = ' '.join([str(elem) for elem in unique_contact])
             unique_resid = atomsel(selection="index "+unique_index, molid=self.molID, frame=frame)
             unique_contact=list(set(unique_resid.resid))
           
             for atom in unique_contact:
                 resid = atomsel(selection="resid "+str(atom), molid=self.molID, frame=frame)
                 resid=list(set(resid.resname))
                 get_one_letter_amino=aminoacid[resid[0]]
                 if str(atom)+get_one_letter_amino in dict_resid_count.keys():
                     #print("existe residuo")
                     dict_resid_count[str(atom)+get_one_letter_amino]=dict_resid_count[str(atom)+get_one_letter_amino]+1
                 else:
                     dict_resid_count[str(atom)+get_one_letter_amino]=1
        return dict_resid_count     
        
    def porcentaje_contact_fit(self,atomselect1,atomselect2,first,last):
        """
        insercion en membrane del residuo a los Posfato de la membrana, tambien devuelve el frame
        la primera insercion
        1) get prot and name P
        2) detectar cabezas fosfolipidas de la membrana y clasificar arriba y abajo
        2.1) the mass of center of membrane was obtained to get a measure than who is up and down, so 
        2.2) name P > centro de masa de menbrana , name P de arriba, caso contrario es abajo.
        2.3) obtener coord x,z y realizar un ajutes de minimos cuadrados sobre los puntos
        2.4) la coord pred(y) sera utilizada como limite para calcular la colision. 
        2.5) usar los coeficientes para calcular la distancia del punto a la recta.
        2.51) evaluo la distancia en x en la funcion polyfit 1 y los comparo con mi z del name P    
        2.6) todo lo que este dentro de las rectas es colision
        """
        
        
        protein = atomsel(selection=atomselect1, molid=self.molID, frame=0) 
        #reference1 = atomsel(selection="protein", molid=self.molID, frame=0) 
             #num residuesatomselect1
        #num_residues=len(pd.factorize(protein.resid)[1])
        ocurrencias_vector=[]
        ocurrencias_temp=[]
        frame_first_insertion=0      
        for frame in range(first,last):
             protein  = atomsel(selection=atomselect1, molid=self.molID, frame=frame) 
             # use frame 0 for the reference
             #sel1 = atomsel(selection=atomselect1, molid=self.molID, frame=frame) 
             sel2 =atomsel(selection=atomselect2, molid=self.molID, frame=frame)
             resid_center=protein.centerperresidue()
             if frame==first:
                 ocurrencias_vector= np.zeros(len(resid_center))
                 ocurrencias_temp=np.zeros(len(resid_center))
        
        
            
             resid_z=list(map(lambda x: x, resid_center))
             nameP_z=list(map(lambda x: x[2], sel2.minmax()))
             #print(nameP_z)
             membrane=atomsel(selection="resname POPC POPG", molid=self.molID, frame=frame)
             mem_mass_center_Z=membrane.center(membrane.mass)[2]
             #print (mem_mass_center_Z)
             
             #outerleaf
             name_P_down=atomsel(selection=atomselect2+" and z<"+str(mem_mass_center_Z) ,molid=self.molID, frame=frame)
             
             x_down=np.asarray(name_P_down.x) 
             y_down=np.asarray(name_P_down.z) 
             x_down=x_down.transpose() 
             y_down=y_down.transpose() 
             p_down = np.polyfit(x_down, y_down, 1)
             #print (p)
             #y_down_ajuste = p_down[0]*x_down + p_down[1]
             
             #inner leaf
             name_P_up=atomsel(selection=atomselect2+" and z>"+str(mem_mass_center_Z) ,molid=self.molID, frame=frame)
             x_up = name_P_up.x
             y_up = name_P_up.z
            
             x_up=np.asarray(x_up) 
             y_up=np.asarray(y_up) 
         
             p_up = np.polyfit(x_up, y_up, 1)
             #print (p)
             #y_up_ajuste = p_up[0]*x_up + p_up[1]
             
             
             #distancia=  (p[0]*x + p[1])/(np.sqrt(p[0]**2+p[1]**2))
             #print (distancia)
             #x.reshape((1, -1))
             #y.reshape((1,-1))
             #model = LinearRegression()
             #model.fit(x, y)
             
             #axis_z_ajust=model.predict(x)
             #print(statistics.mean(axis_z_ajust[0]))
             #print('Coefficients: \n', model.coef_)
             #recta=A*x+B*x
             #print (len(resid_z))
             for number_resid in range(len(resid_z)):
                 x_residue=resid_z[number_resid][0]
                 y_down_ajuste = p_down[0]*x_residue + p_down[1]
                 y_up_ajuste = p_up[0]*x_residue + p_up[1]
                 
                 if resid_z[number_resid][2]<y_up_ajuste and resid_z[number_resid][2]>y_down_ajuste:
                     if np.count_nonzero(ocurrencias_vector)==0:
                        frame_first_insertion=frame
                        #print ("choca"+str(number_resid)+' '+str(frame))
                     
                     ocurrencias_temp[number_resid]=1 
                 else:
                     ocurrencias_temp[number_resid]=0
                     #print ("no choca"+str(number_resid)+' '+str(frame))
             ocurrencias_vector+=ocurrencias_temp
             #print (ocurrencias_temp)
             ocurrencias_temp= np.zeros(len(resid_center))    
             #(>) y “menor que” (<). C
             #    sel2 =atomsel(selection=atomselect2, molid=self.molID, frame=frame)
             #sel1_z = np.array((sel1.center(sel1.mass)[2]))
             #sel2_z= np.array((sel2.center(sel2.mass)[2]))
             #distance_mass_weight.append(np.linalg.norm(sel2_z-sel1_z))
    
        return ocurrencias_vector/(last-first)
        
    def insertions_reach_membrane(self,atomselect1,atomselect2,first,last):
        """
        insercion en membrane del residuo a los Posfato de la membrana, tambien devuelve el frame
        la primera insercion
        1) get prot and name P
        2) detectar cabezas fosfolipidas de la membrana y clasificar arriba y abajo
        2.1) the mass of center of membrane was obtained to get a measure than who is up and down, so 
        2.2) name P > centro de masa de menbrana , name P de arriba, caso contrario es abajo.
        2.3) obtener coord x,z y realizar un ajutes de minimos cuadrados sobre los puntos
        2.4) la coord pred(y) sera utilizada como limite para calcular la colision. 
        2.5) usar los coeficientes para calcular la distancia del punto a la recta.
        2.51) evaluo la distancia en x en la funcion polyfit 1 y los comparo con mi z del name P    
        2.6) todo lo que este dentro de las rectas es colision
        """
        
        
        protein = atomsel(selection=atomselect1, molid=self.molID, frame=0) 
        #reference1 = atomsel(selection="protein", molid=self.molID, frame=0) 
             #num residuesatomselect1
        #num_residues=len(pd.factorize(protein.resid)[1])
        ocurrencias_vector=[]
        ocurrencias_temp=[]
        frame_first_insertion=0      
        for frame in range(first,last):
             protein  = atomsel(selection=atomselect1, molid=self.molID, frame=frame) 
             # use frame 0 for the reference
             #sel1 = atomsel(selection=atomselect1, molid=self.molID, frame=frame) 
             sel2 =atomsel(selection=atomselect2, molid=self.molID, frame=frame)
             resid_center=protein.centerperresidue()
             #center de masa del la proteina
             #protein_center=atomsel(selection="protein", molid=self.molID, frame=frame)
             protein_center_Z=protein.center(protein.mass)[2]
             
             if frame==first:
                 ocurrencias_vector= np.zeros(len(resid_center))
                 ocurrencias_temp=np.zeros(len(resid_center))
        
        
            
             resid_z=list(map(lambda x: x, resid_center))
             nameP_z=list(map(lambda x: x[2], sel2.minmax()))
             #print(nameP_z)
             membrane=atomsel(selection="resname POPC POPG", molid=self.molID, frame=frame)
             mem_mass_center_Z=membrane.center(membrane.mass)[2]
             #print (mem_mass_center_Z)
             
             #outerleaf
             name_P_down=atomsel(selection=atomselect2+" and z<"+str(mem_mass_center_Z) ,molid=self.molID, frame=frame)
             
             x_down=np.asarray(name_P_down.x) 
             y_down=np.asarray(name_P_down.z) 
             x_down=x_down.transpose() 
             y_down=y_down.transpose() 
             p_down = np.polyfit(x_down, y_down, 1)
             #print (p)
             #y_down_ajuste = p_down[0]*x_down + p_down[1]
             
             #inner leaf
             name_P_up=atomsel(selection=atomselect2+" and z>"+str(mem_mass_center_Z) ,molid=self.molID, frame=frame)
             x_up = name_P_up.x
             y_up = name_P_up.z
            
             x_up=np.asarray(x_up) 
             y_up=np.asarray(y_up) 
         
             p_up = np.polyfit(x_up, y_up, 1)
             #print (p)
             #y_up_ajuste = p_up[0]*x_up + p_up[1]
             
             
             #distancia=  (p[0]*x + p[1])/(np.sqrt(p[0]**2+p[1]**2))
             #print (distancia)
             #x.reshape((1, -1))
             #y.reshape((1,-1))
             #model = LinearRegression()
             #model.fit(x, y)
             
             #axis_z_ajust=model.predict(x)
             #print(statistics.mean(axis_z_ajust[0]))
             #print('Coefficients: \n', model.coef_)
             #recta=A*x+B*x
             #print (len(resid_z))
             for number_resid in range(len(resid_z)):
                 x_residue=resid_z[number_resid][0]
                 y_down_ajuste = p_down[0]*x_residue + p_down[1]
                 y_up_ajuste = p_up[0]*x_residue + p_up[1]
                 
                 if resid_z[number_resid][2]<y_up_ajuste and resid_z[number_resid][2]>y_down_ajuste:
                     #if np.count_nonzero(ocurrencias_vector)==0:
                    if protein_center_Z<y_up_ajuste+5 and protein_center_Z> y_down_ajuste-5: 
                        frame_first_insertion=frame
                        print ("choca",str(number_resid),str(frame),protein_center_Z,y_up_ajuste,y_down_ajuste)
                        return frame_first_insertion
                     
        return frame_first_insertion
    
    def get_residue_insert_mem_per_frame(self,frame,membrane,fosfolipid_head_select,protein_select):
        inserted_residues=[]
        protein  = atomsel(selection=protein_select, molid=self.molID, frame=frame) 
        resid_center=protein.centerperresidue()
        resid_z=list(map(lambda x: x, resid_center))     
        membrane  = atomsel(selection=membrane, molid=self.molID, frame=frame)
        #outerleaf
        name_P_down=atomsel(selection=fosfolipid_head_select+" and z<"+str(membrane.center(membrane.mass)[2]) ,molid=self.molID, frame=frame) 
        #inner leaf
        name_P_up=atomsel(selection=fosfolipid_head_select+" and z>"+str(membrane.center(membrane.mass)[2]) ,molid=self.molID, frame=frame)
             
        x_down=np.asarray(name_P_down.x).transpose() 
        y_down=np.asarray(name_P_down.z).transpose()
        p_down = np.polyfit(x_down, y_down, 1)
             
        x_up=np.asarray(name_P_up.x).transpose()  
        y_up=np.asarray(name_P_up.z).transpose()  
        p_up = np.polyfit(x_up, y_up, 1)
        
        for number_resid in range(len(resid_z)):
            x_residue=resid_z[number_resid][0]
            z_down_ajuste = p_down[0]*x_residue + p_down[1]
            z_up_ajuste = p_up[0]*x_residue + p_up[1]
            if resid_z[number_resid][2]<z_up_ajuste and resid_z[number_resid][2]>z_down_ajuste:
                    inserted_residues.append(number_resid)
                 
        return inserted_residues
    def porcentaje_contact_free_energy(self,atomselect1,first,last,residue_to_find):
        #calcular la energia libre a partir de eventos de pegado no pegado.
        #criterio:puede tener varios criterios para decir que es pegado o no
        #1) si el parche esta hacia abajo o mas bien cerca de la membrana es 1
        #es decir si resid 10 or 11 esta cerca de men and 29 or  30 esta cerca de nem 
        
        #get selecion at frame 0
        protein = atomsel(selection=atomselect1, molid=self.molID, frame=0) 
        
        #divido el residuos del parche en 2 grupos
        parche1=residue_to_find[0]
        parche2=residue_to_find[1]

        unbound=0
        bound=0      
        for frame in range(first,last):
             protein  = atomsel(selection=atomselect1, molid=self.molID, frame=frame)
             inserted_residues=Trajectory.get_residue_insert_mem_per_frame(self,frame,"resname POPC POPG","name P","protein")
             #print (inserted_residues)
        
             

             if not protein.resid:
                 #print ("vacio frame "+str(frame))
                 unbound=unbound+1
             else: 
                
                residues_in_contact_membrane=sorted(protein.resid+inserted_residues)
                residues_in_contact_membrane=np.unique(residues_in_contact_membrane).tolist()
                
                print  (residues_in_contact_membrane)
                parche1_ON=0
                parche2_ON=0
                for ind in residues_in_contact_membrane:
                    if (ind in parche1):
                        parche1_ON=1
                        continue
                    elif (ind in parche2):
                        parche2_ON=1
                        continue
                if (parche1_ON==1 and parche2_ON==1):
                    bound=bound+1
                    #print (residues_in_contact_membrane)
                    print (frame)
                        
                else:
                    unbound=unbound+1
                residues_in_contact_membrane=[]                
        return [unbound,bound]


    def porcentaje_contact_fit_list(self,atomselect1,atomselect2,first,last):
        """
        este metodo la colision  del residuo a la cabezas de Posfato de la membrana.
        1) get prot and name P
        2) detectar cabezas fosfolipidas de la membrana y clasificar arriba y abajo
        2.1) the mass of center of membrane was obtained to get a measure than who is up and down, so 
        2.2) name P > centro de masa de menbrana , name P de arriba, caso contrario es abajo.
        2.3) obtener coord x,z y realizar un ajutes de minimos cuadrados sobre los puntos
        2.4) la coord pred(y) sera utilizada como limite para calcular la colision. 
        2.5) usar los coeficientes para calcular la distancia del punto a la recta.
        2.51) evaluo la distancia en x en la funcion ax+b= z y los comparo con mi z del name P    
        2.6) todo lo que este dentro de las rectas es colision
        """
        
        
        protein = atomsel(selection=atomselect1, molid=self.molID, frame=0) 
        #reference1 = atomsel(selection="protein", molid=self.molID, frame=0) 
             #num residuesatomselect1
        #num_residues=len(pd.factorize(protein.resid)[1])
        ocurrencias_vector=[]
        ocurrencias_temp=[]
              
        for frame in range(first,last):
             protein  = atomsel(selection=atomselect1, molid=self.molID, frame=frame) 
             # use frame 0 for the reference
             #sel1 = atomsel(selection=atomselect1, molid=self.molID, frame=frame) 
             sel2 =atomsel(selection=atomselect2, molid=self.molID, frame=frame)
             resid_center=protein.centerperresidue()
             if frame==first:
                 ocurrencias_vector= np.zeros(len(resid_center))
                 ocurrencias_temp=np.zeros(len(resid_center))
        
        
            
             resid_z=list(map(lambda x: x, resid_center))
             nameP_z=list(map(lambda x: x[2], sel2.minmax()))
             #print(nameP_z)
             membrane=atomsel(selection="resname POPC POPG", molid=self.molID, frame=frame)
             mem_mass_center_Z=membrane.center(membrane.mass)[2]
             #print (mem_mass_center_Z)
             
             #outerleaf
             name_P_down=atomsel(selection="name P and z<"+str(mem_mass_center_Z) ,molid=self.molID, frame=frame)
             
             x_down=np.asarray(name_P_down.x) 
             y_down=np.asarray(name_P_down.z) 
             x_down=x_down.transpose() 
             y_down=y_down.transpose() 
             p_down = np.polyfit(x_down, y_down, 1)
             #print (p)
             #y_down_ajuste = p_down[0]*x_down + p_down[1]
             
             #inner leaf
             name_P_up=atomsel(selection="name P and z>"+str(mem_mass_center_Z) ,molid=self.molID, frame=frame)
             x_up = name_P_up.x
             y_up = name_P_up.z
            
             x_up=np.asarray(x_up) 
             y_up=np.asarray(y_up) 
             x_up=x_up.transpose() 
             y_up=y_up.transpose() 
             p_up = np.polyfit(x_up, y_up, 1)
             #print (p)
             #y_up_ajuste = p_up[0]*x_up + p_up[1]
             
             
             #distancia=  (p[0]*x + p[1])/(np.sqrt(p[0]**2+p[1]**2))
             #print (distancia)
             #x.reshape((1, -1))
             #y.reshape((1,-1))
             #model = LinearRegression()
             #model.fit(x, y)
             
             #axis_z_ajust=model.predict(x)
             #print(statistics.mean(axis_z_ajust[0]))
             #print('Coefficients: \n', model.coef_)
             #recta=A*x+B*x
             #print (len(resid_z))
             for number_resid in range(len(resid_z)):
                 x_residue=resid_z[number_resid][0]
                 y_down_ajuste = p_down[0]*x_residue + p_down[1]
                 y_up_ajuste = p_up[0]*x_residue + p_up[1]
                 
                 if resid_z[number_resid][2]<y_up_ajuste and resid_z[number_resid][2]>y_down_ajuste:
                     #print ("choca"+str(number_resid)+' '+str(frame))
                     
                     ocurrencias_temp[number_resid]=1 
                 else:
                     ocurrencias_temp[number_resid]=0
                     #print ("no choca"+str(number_resid)+' '+str(frame))
             ocurrencias_vector+=ocurrencias_temp
             #print (ocurrencias_temp)
             ocurrencias_temp= np.zeros(len(resid_center))    
             #(>) y “menor que” (<). C
             #    sel2 =atomsel(selection=atomselect2, molid=self.molID, frame=frame)
             #sel1_z = np.array((sel1.center(sel1.mass)[2]))
             #sel2_z= np.array((sel2.center(sel2.mass)[2]))
             #distance_mass_weight.append(np.linalg.norm(sel2_z-sel1_z))
    
        return ocurrencias_vector
    
    def time_contact(self,atomselect1,atomselect2):
        protein = atomsel(selection=atomselect1, molid=self.molID, frame=0) 
        #reference1 = atomsel(selection="protein", molid=self.molID, frame=0) 
             #num residuesatomselect1
        #num_residues=len(pd.factorize(protein.resid)[1])
        ocurrencias_vector=[]
        ocurrencias_temp=[]
              
        for frame in range(Trajectory.num_frames(self)):
             protein  = atomsel(selection="protein", molid=self.molID, frame=frame) 
             # use frame 0 for the reference
             #sel1 = atomsel(selection=atomselect1, molid=self.molID, frame=frame) 
             sel2 =atomsel(selection=atomselect2, molid=self.molID, frame=frame)
             resid_center=protein.centerperresidue()
             if frame==0:
                 ocurrencias_vector= np.zeros(len(resid_center))
                 ocurrencias_temp=np.zeros(len(resid_center))
        
        
            
             resid_z=list(map(lambda x: x[2], resid_center))
             nameP_z=list(map(lambda x: x[2], sel2.minmax()))
             for number_resid in range(len(resid_z)):
                 if resid_z[number_resid]>0.8*nameP_z[0] and resid_z[number_resid]<nameP_z[1]*0.8:
                     #print ("choca"+str(number_resid)+' '+str(frame))
                     ocurrencias_temp[number_resid]=1 
                 else:
                     ocurrencias_temp[number_resid]=0
                     #print ("no choca"+str(number_resid)+' '+str(frame))
             
             #np.append(ocurrencias_vector,ocurrencias_temp)
             #print (ocurrencias_temp)
             #np.append(ocurrencias_vector, ocurrencias_temp, axis = 1)
             ocurrencias_vector=np.vstack((ocurrencias_vector, ocurrencias_temp))
             ocurrencias_temp= np.zeros(len(resid_center))    
             #(>) y “menor que” (<). C
             #    sel2 =atomsel(selection=atomselect2, molid=self.molID, frame=frame)
             #sel1_z = np.array((sel1.center(sel1.mass)[2]))
             #sel2_z= np.array((sel2.center(sel2.mass)[2]))
             #distance_mass_weight.append(np.linalg.norm(sel2_z-sel1_z))
             #ocurrencias_vector=ocurrencias_vector.tolist()
             #type(ocurrencias_vector)
             #ocurrencias_vector.append(ocurrencias_temp)
                 
        return ocurrencias_vector
    
    def rmsd_time(self,atomselect):
         rmsd_array=[]
         # use frame 0 for the reference
         reference = atomsel(selection=atomselect, molid=self.molID, frame=0) 
         #compare = atomsel(atomselect)
         #set reference [atomselect $mol "protein" frame 0]
         # the frame being compared 
         #set compare [atomselect $mol "protein"]
         for frame in range(Trajectory.num_frames(self)):
             #protein  = atomsel(selection="protein", molid=molid, frame=frame) 
              # the frame being compared 
             
             compare = atomsel(selection=atomselect, molid=self.molID, frame=frame) 
             #set trans_mat [measure fit $compare $reference]
             trans_mat=atomsel.fit(compare,reference)
             # do the alignment
             compare.move(trans_mat)
             #$compare move $trans_mat
             # compute the RMSD
             #set rmsd [measure rmsd $compare $reference]
             rmsd_array.append(atomsel.rmsd(compare,reference))
    
         return rmsd_array
    def rmsd_time_references(self,atomselect,references1):
         rmsd_array=[]
         # use frame 0 for the reference
         reference = atomsel(selection=atomselect, molid=self.molID, frame=0)
         reference2 = atomsel(selection=references1, molid=self.molID, frame=0) 
         #compare = atomsel(atomselect)
         #set reference [atomselect $mol "protein" frame 0]
         # the frame being compared 
         #set compare [atomselect $mol "protein"]
         for frame in range(Trajectory.num_frames(self)):
             #protein  = atomsel(selection="protein", molid=molid, frame=frame) 
              # the frame being compared 
             compare = atomsel(selection=atomselect, molid=self.molID, frame=frame) 
             
             compare2 = atomsel(selection=references1, molid=self.molID, frame=frame) 
             #set trans_mat [measure fit $compare $reference]
             trans_mat=atomsel.fit(compare,reference)
             # do the alignment
             #compare.move(trans_mat)
             compare2.move(trans_mat)
             
             #$compare move $trans_mat
             # compute the RMSD
             #set rmsd [measure rmsd $compare $reference]
             rmsd_array.append(atomsel.rmsd(compare2,reference2))
    
         return rmsd_array
    @staticmethod     
    def promedio_ocurrencias(*args):
        #ocurrence= np.array(args)
        ocurrence=np.transpose(args)
        value = 0
        mean_residues=[]
        std_error_residues=[]
        for x in ocurrence:
            
            std_error_residues.append(np.std(x, ddof=1) / np.sqrt(np.size(x)))
            mean_residues.append(np.mean(x))
            #0.42:0.1805
            #print (np.sum(x[0]))
        return [mean_residues,std_error_residues]
    @staticmethod     
    def set_beta_factor(residue_occupancy:list, query,molid:int,output:str):
        """
        metodo que colorea por beta factor
        residue_occupancy: vector con valores a colocar en beta columm
        query: atomselecion "protein"
        output "pbd salida"
        molid: id de la proteina  
        """
        protein = atomsel(query)
        for x in protein.resid:
            atomsel(query+" and resid "+str(x)).beta=residue_occupancy[x-1]
            print(x,residue_occupancy[x-1])
        protein.write("pdb",output)
        print("cerrando molecula ",molid) 
        molecule.delete(molid)        
    def rmsf_residue(self,atomselect,traj1):
        rmsd_array=[]
        

     
        reference = atomsel(selection=atomselect, molid=traj1.molID)
        
        rmsf=reference.rmsfperresidue(first=1, last=traj1.num_frames()-1)
         #rmsf+=compare2.rmsf(frame) 
        #reference1 = atomsel(selection="protein", molid=self.molID, frame=0) 
        #num residues
        #um_residues=len(pd.factorize(reference.resid)[1])
        #rmsf = np.zeros(num_residues)
        
        #compare2 = atomsel(atomselect)
        #set reference [atomselect $mol "protein" frame 0]
        # the frame being compared 
        #set compare [atomselect $mol "protein"]
       
        #newList = [x /Trajectory.num_frames(self)  for x in rmsf]
        #newList = map(lambda rmsf: rmsf/int(max(rmsf)), rmsf)
        #rmsf[:]=[rmsf / int(max(rmsf)) for x in rmsf]
        return rmsf
    def rmsd_vmd(self,atomselect):
         rmsd_array=[]
         # use frame 0 for the reference
         reference = atomsel(selection=atomselect, molid=self.molID, frame=0)
         reference1 = atomsel(selection=atomselect, molid=self.molID, frame=1)  
         x=reference.rmsd(selection=reference1)
         print (x) 
         #compare = atomsel(atomselect)
         #set reference [atomselect $mol "protein" frame 0]

    def rmsf_time2(self,atomselect):
            
            
            #Leftraru2018
            
            # use frame 0 for the reference
            reference = atomsel(selection=atomselect, molid=self.molID, frame=0) 
            #reference1 = atomsel(selection="protein", molid=self.molID, frame=0) 
            #num residues
            num_residues=len(pd.factorize(reference.resid)[1])
            print (num_residues)
            rmsf = np.zeros(num_residues)
            
            #compare2 = atomsel(atomselect)
            #set reference [atomselect $mol "protein" frame 0]
            # the frame being compared 
            #set compare [atomselect $mol "protein"]
            mask = vmdnumpy.atomselect(molid=self.molID, frame=0,selection=atomselect)
            ref = np.compress(mask, vmdnumpy.timestep(self.molID, 0), axis=0)
            
            for frame in range(Trajectory.num_frames(self)):
                 #protein  = atomsel(selection="protein", molid=molid, frame=frame) 
                  # the frame being compared 
                 frame = np.compress(mask, vmdnumpy.timestep(self.molID, frame), axis=0)
                 rmsf += np.sqrt(np.sum((frame-ref)**2, axis=1))
                 #compare = atomsel(selection=atomselect, molid=self.molID, frame=frame) 
                 #set trans_mat [measure fit $compare $reference]
                 #trans_mat=atomsel.fit(compare,reference)
                 # do the alignment
                 #compare.move(trans_mat)
                 #$compare move $trans_mat
                 #compute the RMSD
                 #set rmsd [measure rmsd $compare $reference]
                 #rmsf+=compare2.rmsf(frame)
            rmsf /= float(Trajectory.num_frames(self)-0)
            rmsf = np.sqrt(rmsf)     
            #newList = [x /Trajectory.num_frames(self)  for x in rmsf]
            #newList = map(lambda rmsf: rmsf/int(max(rmsf)), rmsf)
            #rmsf[:]=[rmsf / int(max(rmsf)) for x in rmsf]
            return rmsf
    def rmsf_ligand(self,atomselect):
            # rmsf para ligando no residues
            # use frame 0 for the reference
            reference = atomsel(selection=atomselect, molid=self.molID, frame=0) 
            #reference1 = atomsel(selection="protein", molid=self.molID, frame=0) 
            #num residues
           
            rmsf = np.zeros(len(reference.index))
            
            #compare2 = atomsel(atomselect)
            #set reference [atomselect $mol "protein" frame 0]
            # the frame being compared 
            #set compare [atomselect $mol "protein"]
            mask = vmdnumpy.atomselect(molid=self.molID, frame=0,selection=atomselect)
            ref = np.compress(mask, vmdnumpy.timestep(self.molID, 0), axis=0)
            
            for frame in range(Trajectory.num_frames(self)):
                 #protein  = atomsel(selection="protein", molid=molid, frame=frame) 
                  # the frame being compared 
                 frame = np.compress(mask, vmdnumpy.timestep(self.molID, frame), axis=0)
                 rmsf += np.sqrt(np.sum((frame-ref)**2, axis=1))
                 #compare = atomsel(selection=atomselect, molid=self.molID, frame=frame) 
                 #set trans_mat [measure fit $compare $reference]
                 #trans_mat=atomsel.fit(compare,reference)
                 # do the alignment
                 #compare.move(trans_mat)
                 #$compare move $trans_mat
                 #compute the RMSD
                 #set rmsd [measure rmsd $compare $reference]
                 #rmsf+=compare2.rmsf(frame)
            rmsf /= float(Trajectory.num_frames(self)-0)
            rmsf = np.sqrt(rmsf)     
            #newList = [x /Trajectory.num_frames(self)  for x in rmsf]
            #newList = map(lambda rmsf: rmsf/int(max(rmsf)), rmsf)
            #rmsf[:]=[rmsf / int(max(rmsf)) for x in rmsf]
            return rmsf    
    def SASA_percent(self,complejo,ligand_sel):
        #get percent solvent accssibility of a ligand
        sasa_vector=[]
        for frame in range(Trajectory.num_frames(self)):
                 lig_sel  = atomsel(selection=ligand_sel, molid=self.molID, frame=frame) 
                 complejo_lig_prot  = atomsel(selection=complejo, molid=self.molID, frame=frame) 
                 
                 lig_sel_sasa = lig_sel.sasa(srad=1.4, restrict=lig_sel)
                 complejo_lig_prot_sasa = complejo_lig_prot.sasa(srad=1.4, restrict=complejo_lig_prot)
                 if frame%100:
                     print("frame "+str(frame)+"  " +str(complejo_lig_prot_sasa )+"-" +str(lig_sel_sasa ))
                 sasa_vector.append(complejo_lig_prot_sasa - lig_sel_sasa)
        return sasa_vector
    def Water_Hbond(self,donor,aceptor,output,cutoff,sheet,angle,mode,water,only_polar=True):
        #get percent solvent accssibility of a ligand
        H_bond_vector=[]
        dict_hbond={}
        summary_hbound={}
        number_hbond={}
        for frame in range(Trajectory.num_frames(self)):
                 sel_1  = atomsel(selection=donor, molid=self.molID, frame=frame) 
                 sel_2  = atomsel(selection=aceptor, molid=self.molID, frame=frame) 
                 lig_sel_sasa = sel_1.hbonds(cutoff=cutoff,maxangle=angle,acceptor=sel_2)
                 #H_bond_vector.append(lig_sel_sasa)
                 number_hbond[frame]=0
                 dict_aux={}
                 for data in range(len(lig_sel_sasa[1])):
                     sel_A  = atomsel(selection="index "+str(lig_sel_sasa[0][data]), molid=self.molID, frame=frame)
                     sel_B  = atomsel(selection="index "+str(lig_sel_sasa[1][data]), molid=self.molID, frame=frame)
                     sel_H  = atomsel(selection="index "+str(lig_sel_sasa[2][data]), molid=self.molID, frame=frame)
                     
                     if water=="acceptor":
                         ligand_atom=str(sel_B.resname[0])+"-"+str(sel_B.name[0])+"---"+str(sel_A.resname[0])+"-"+str(sel_A.name[0])+"-"+str(sel_A.index[0])
                     elif water=="donor":
                         ligand_atom=str(sel_A.resname[0])+"-"+str(sel_A.name[0])+"---"+str(sel_B.resname[0])+"-"+str(sel_B.name[0])+"-"+str(sel_B.index[0])
                     else:
                         return print("water, only canbe donor or acceptor")
                     if only_polar==True:
                         polar_allowing=['N','O','S','F']
                         #print ("only polar active",polar_allowing)
                         if any(x in sel_A.type[0][0] for x in polar_allowing) and any(x in sel_B.type[0][0] for x in polar_allowing):
                             number_hbond[frame]=number_hbond[frame]+1
                             if ligand_atom not in dict_aux.keys():
                                dict_aux[ligand_atom]=1
                           
                      
                     else:
                         if ligand_atom not in dict_aux.keys():
                            dict_aux[ligand_atom]=1
                         number_hbond[frame]=number_hbond[frame]+1
                 dict_aux = collections.Counter(dict_aux)
                 summary_hbound = collections.Counter(summary_hbound)        
                 summary_hbound=summary_hbound+dict_aux    
                            
        for key in summary_hbound:
            summary_hbound[key]= (summary_hbound[key]/ Trajectory.num_frames(self)*100)               
        self.to_excel(output=output,resume=summary_hbound,numHbound=number_hbond,sheet=sheet, mode=mode)               
    def Dewetting(self,sel1,sel2, first, last):
        dew_vector=[]
        radio_vect=Trajectory.average_radius_of_gyration(self, sel1)
        print (radio_vect)

        for frame in range(first,last):
                 protein  = atomsel(selection=sel1, molid=self.molID, frame=frame) 
                 query_name_p=sel2+" and ( x>= "+str(protein.center()[0]-radio_vect)+" and x<= "+str(protein.center()[0]+radio_vect)+")"
                 radio_giro=protein.rgyr(protein.mass)
                 print (radio_giro)
                 print (query_name_p)
                 #print (protein.mass)
                 print (protein.minmax())
                 zmax=protein.minmax()[0][2]
                 
                 name_p  = atomsel(selection=query_name_p, molid=self.molID, frame=frame)
                 
                 
                 print (name_p.minmax())
                 zmin=name_p.minmax()[1][2]
                 print (zmin)
                 x=protein.center(protein.mass)[0]
                 y=protein.center(protein.mass)[1]
                 z=protein.center(protein.mass)[2]
                 cilinder_water= "water and  same residue as (((x-"+str(x)+")^2" +"+(y-"+str(y)+")^2 <="+str(radio_giro)+"^2)and z >="+str(zmin)+" and z <="+str(zmax)+")"
                 print (cilinder_water)
                 agua  = atomsel(selection=cilinder_water, molid=self.molID, frame=frame) 
                 print (len(agua.index))
                 num_agua=len(agua.index)
                 h=np.abs(zmin-zmax)
                 volumen=np.pi*radio_giro**2*h
                 densidad_agua=num_agua/volumen
                 #>>> big_sel = atomsel('protein or resname LIG')
                 #>>> lig_sel = atomsel('resname LIG')
                 #ligand_in_protein_sasa = big_sel.sasa(srad=1.4, restrict=lig_sel)
                 #ligand_in_protein_sasa = big_sel.sasa(srad=1.4)
                 
                 #ligand_alone_sasa= lig_sel.sasa(srad=1.4, points=True)
                 if (frame % 1000)==0:
                     print("frame "+str(frame)+"  " +str(protein.center()))
                 dew_vector.append(densidad_agua)
        return dew_vector              
    def pdbs_from_dcd(self,atomselect1,output):
        vector_radio=[]
        for frame in range(Trajectory.num_frames(self)):
            
            #protein  = atomsel(selection=atomselect1, molid=self.molID, frame=frame) 
            #last_frame = molecule.numframes - 1
            #atomsel.write()
            molecule.write(self.molID, "pdb", "last_frame.pdb", first=frame)
            #vector_radio.append(radio_giro)