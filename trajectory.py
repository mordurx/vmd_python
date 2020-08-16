#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 15:33:39 2019

@author: eniac
"""

import subprocess
from vmd import molecule,atomsel,vmdnumpy,animate
import numpy as np
import pandas as pd
import statistics
from sklearn import preprocessing
from sklearn.linear_model import LinearRegression
class Trajectory:
    
    molID=0
    def __init__(self, psf,dcd,first=0,last=-1,stride=1,waitfor=-1):
        self.psf = psf
        self.dcd = dcd
        try:
        
            Trajectory.molID=molecule.new('new2')
            Trajectory.molID=molecule.load('psf',psf) # load trajectory
            molecule.read(Trajectory.molID,'dcd',dcd,stride=stride,first=first,last=last,waitfor=waitfor) 
            print (Trajectory.molID)
        except IOError:
            print ("Could not read dcd file or psf:", dcd)
            raise Exception()
    def mean_displacement(self,atomselect1):
       
        displacement_x=[]
        displacement_y=[]
        displacement_z=[]
        
        #sampleo=int(np.floor ((Trajectory.num_frames(self)-1)/2))
        promedio_x=[]
        promedio_y=[]
        promedio_z=[]
        sel0 = atomsel(selection=atomselect1, molid=Trajectory.molID, frame=0) 
        count=0   
        for frame in range(0, 1000,1):
            print (frame)
            for frame1 in range( Trajectory.num_frames(self)-5000, Trajectory.num_frames(self)-1):
                sel0 = atomsel(selection=atomselect1, molid=Trajectory.molID, frame=frame1-1) 
               
                sel1 = atomsel(selection=atomselect1, molid=Trajectory.molID, frame=frame1+frame) 
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
        
#        for frame in vector_pdb:
#            
#            molID=molecule.load('pdb',frame)
#            print (animate.is_active(molID)) 
#            #molecule.read(molID,'pdb',frame)
#            sel1 = atomsel(selection="segname TOX",molid=molID) 
#        molecule.listall()
#        molecule.cancel
#        sel1.write('dcd','vsd1000.dcd')
#       
        #return molecule.listall()

        
    def close (self):
        print("cerrando molecule "+str(self.molID))
        molecule.cancel(self.molID)
        molecule.delete(self.molID)
         
    
    def num_frames(self):
        return molecule.numframes(Trajectory.molID)
    
       
    def get_molid(self):
        return Trajectory.molID
    
    def porcentaje_contact_radio_gyrations(self,atomselect1,atomselect2):
        for frame in range(Trajectory.num_frames(self)):
            #protein  = atomsel(selection="protein", molid=molid, frame=frame) 
            sel1 = atomsel(selection=atomselect1, molid=Trajectory.molID, frame=frame) 
            sel2 =atomsel(selection=atomselect2, molid=Trajectory.molID, frame=frame)
            sel1_z = np.array((sel1.center(sel1.mass)[2]))
            sel2_z= np.array((sel2.center(sel2.mass)[2]))
            distance_mass_weight.append(np.linalg.norm(sel2_z-sel1_z))
            #(x−a)2 + (y−b)2 = r2
    
    
    def average_radius_of_gyration(self,atomselect1):
        vector_radio=[]
        for frame in range(Trajectory.num_frames(self)):
            
            protein  = atomsel(selection=atomselect1, molid=Trajectory.molID, frame=frame) 
            radio_giro=protein.rgyr(protein.mass)
            vector_radio.append(radio_giro)
            #sel1 = atomsel(selection=atomselect1, molid=Trajectory.molID, frame=frame) 
            #sel2 =atomsel(selection=atomselect2, molid=Trajectory.molID, frame=frame)
            #sel1_z = np.array((sel1.center(sel1.mass)[2]))
            #sel2_z= np.array((sel2.center(sel2.mass)[2]))
            #distance_mass_weight.append(np.linalg.norm(sel2_z-sel1_z))
         
        return vector_radio  
    
    def membrane_center_mass(self,atomselect1):
        membrane_center_mass=[]
        for frame in range(Trajectory.num_frames(self)):
            #protein  = atomsel(selection="protein", molid=molid, frame=frame) 
        
            sel1 = atomsel(selection=atomselect1, molid=Trajectory.molID, frame=frame) 
            sel1_z = np.array((sel1.center(sel1.mass)[2]))
            membrane_center_mass.append(np.linalg.norm(sel1_z))
        return membrane_center_mass     
    def distance_center_mass(self,atomselect1,atomselect2):
        distance_mass_weight=[]
        for frame in range(Trajectory.num_frames(self)):
            #protein  = atomsel(selection="protein", molid=molid, frame=frame) 
        
            sel1 =atomsel(selection=atomselect1, molid=Trajectory.molID, frame=frame) 
            sel2 =atomsel(selection=atomselect2, molid=Trajectory.molID, frame=frame)
            sel1_z = np.array((sel1.center(sel1.mass)[2]))
            sel2_z= np.array((sel2.center(sel2.mass)[2]))
            
            sel1_x = np.array((sel1.center(sel1.mass)[0]))
            sel2_x= np.array((sel2.center(sel2.mass)[0]))
            
            sel1_y = np.array((sel1.center(sel1.mass)[1]))
            sel2_y= np.array((sel2.center(sel2.mass)[1]))
            
            
            dist = np.sqrt((sel2_x - sel1_x)**2 + (sel2_y - sel1_y)**2+ (sel2_z - sel1_z)**2)  
            distance_mass_weight.append(dist)
    
        return distance_mass_weight
    def porcentaje_contact(self,atomselect1,atomselect2):
        protein = atomsel(selection=atomselect1, molid=Trajectory.molID, frame=0) 
        #reference1 = atomsel(selection="protein", molid=Trajectory.molID, frame=0) 
             #num residuesatomselect1
        #num_residues=len(pd.factorize(protein.resid)[1])
        ocurrencias_vector=[]
        ocurrencias_temp=[]
              
        for frame in range(Trajectory.num_frames(self)):
             protein  = atomsel(selection="protein", molid=Trajectory.molID, frame=frame) 
             # use frame 0 for the reference
             #sel1 = atomsel(selection=atomselect1, molid=Trajectory.molID, frame=frame) 
             sel2 =atomsel(selection=atomselect2, molid=Trajectory.molID, frame=frame)
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
             #    sel2 =atomsel(selection=atomselect2, molid=Trajectory.molID, frame=frame)
             #sel1_z = np.array((sel1.center(sel1.mass)[2]))
             #sel2_z= np.array((sel2.center(sel2.mass)[2]))
             #distance_mass_weight.append(np.linalg.norm(sel2_z-sel1_z))
    
        return ocurrencias_vector/Trajectory.num_frames(self)
    def porcentaje_contact_fit(self,atomselect1,atomselect2,first,last):
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
        
        
        protein = atomsel(selection=atomselect1, molid=Trajectory.molID, frame=0) 
        #reference1 = atomsel(selection="protein", molid=Trajectory.molID, frame=0) 
             #num residuesatomselect1
        #num_residues=len(pd.factorize(protein.resid)[1])
        ocurrencias_vector=[]
        ocurrencias_temp=[]
              
        for frame in range(first,last):
             protein  = atomsel(selection=atomselect1, molid=Trajectory.molID, frame=frame) 
             # use frame 0 for the reference
             #sel1 = atomsel(selection=atomselect1, molid=Trajectory.molID, frame=frame) 
             sel2 =atomsel(selection=atomselect2, molid=Trajectory.molID, frame=frame)
             resid_center=protein.centerperresidue()
             if frame==first:
                 ocurrencias_vector= np.zeros(len(resid_center))
                 ocurrencias_temp=np.zeros(len(resid_center))
        
        
            
             resid_z=list(map(lambda x: x, resid_center))
             nameP_z=list(map(lambda x: x[2], sel2.minmax()))
             #print(nameP_z)
             membrane=atomsel(selection="resname POPC POPG", molid=Trajectory.molID, frame=frame)
             mem_mass_center_Z=membrane.center(membrane.mass)[2]
             #print (mem_mass_center_Z)
             
             #outerleaf
             name_P_down=atomsel(selection="name P and z<"+str(mem_mass_center_Z) ,molid=Trajectory.molID, frame=frame)
             
             x_down=np.asarray(name_P_down.x) 
             y_down=np.asarray(name_P_down.z) 
             x_down=x_down.transpose() 
             y_down=y_down.transpose() 
             p_down = np.polyfit(x_down, y_down, 1)
             #print (p)
             #y_down_ajuste = p_down[0]*x_down + p_down[1]
             
             #inner leaf
             name_P_up=atomsel(selection="name P and z>"+str(mem_mass_center_Z) ,molid=Trajectory.molID, frame=frame)
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
             #    sel2 =atomsel(selection=atomselect2, molid=Trajectory.molID, frame=frame)
             #sel1_z = np.array((sel1.center(sel1.mass)[2]))
             #sel2_z= np.array((sel2.center(sel2.mass)[2]))
             #distance_mass_weight.append(np.linalg.norm(sel2_z-sel1_z))
    
        return ocurrencias_vector/(last-first)
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
        
        
        protein = atomsel(selection=atomselect1, molid=Trajectory.molID, frame=0) 
        #reference1 = atomsel(selection="protein", molid=Trajectory.molID, frame=0) 
             #num residuesatomselect1
        #num_residues=len(pd.factorize(protein.resid)[1])
        ocurrencias_vector=[]
        ocurrencias_temp=[]
              
        for frame in range(first,last):
             protein  = atomsel(selection=atomselect1, molid=Trajectory.molID, frame=frame) 
             # use frame 0 for the reference
             #sel1 = atomsel(selection=atomselect1, molid=Trajectory.molID, frame=frame) 
             sel2 =atomsel(selection=atomselect2, molid=Trajectory.molID, frame=frame)
             resid_center=protein.centerperresidue()
             if frame==first:
                 ocurrencias_vector= np.zeros(len(resid_center))
                 ocurrencias_temp=np.zeros(len(resid_center))
        
        
            
             resid_z=list(map(lambda x: x, resid_center))
             nameP_z=list(map(lambda x: x[2], sel2.minmax()))
             #print(nameP_z)
             membrane=atomsel(selection="resname POPC POPG", molid=Trajectory.molID, frame=frame)
             mem_mass_center_Z=membrane.center(membrane.mass)[2]
             #print (mem_mass_center_Z)
             
             #outerleaf
             name_P_down=atomsel(selection="name P and z<"+str(mem_mass_center_Z) ,molid=Trajectory.molID, frame=frame)
             
             x_down=np.asarray(name_P_down.x) 
             y_down=np.asarray(name_P_down.z) 
             x_down=x_down.transpose() 
             y_down=y_down.transpose() 
             p_down = np.polyfit(x_down, y_down, 1)
             #print (p)
             #y_down_ajuste = p_down[0]*x_down + p_down[1]
             
             #inner leaf
             name_P_up=atomsel(selection="name P and z>"+str(mem_mass_center_Z) ,molid=Trajectory.molID, frame=frame)
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
             #    sel2 =atomsel(selection=atomselect2, molid=Trajectory.molID, frame=frame)
             #sel1_z = np.array((sel1.center(sel1.mass)[2]))
             #sel2_z= np.array((sel2.center(sel2.mass)[2]))
             #distance_mass_weight.append(np.linalg.norm(sel2_z-sel1_z))
    
        return ocurrencias_vector
    
    def time_contact(self,atomselect1,atomselect2):
        protein = atomsel(selection=atomselect1, molid=Trajectory.molID, frame=0) 
        #reference1 = atomsel(selection="protein", molid=Trajectory.molID, frame=0) 
             #num residuesatomselect1
        #num_residues=len(pd.factorize(protein.resid)[1])
        ocurrencias_vector=[]
        ocurrencias_temp=[]
              
        for frame in range(Trajectory.num_frames(self)):
             protein  = atomsel(selection="protein", molid=Trajectory.molID, frame=frame) 
             # use frame 0 for the reference
             #sel1 = atomsel(selection=atomselect1, molid=Trajectory.molID, frame=frame) 
             sel2 =atomsel(selection=atomselect2, molid=Trajectory.molID, frame=frame)
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
             #    sel2 =atomsel(selection=atomselect2, molid=Trajectory.molID, frame=frame)
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
         reference = atomsel(selection=atomselect, molid=Trajectory.molID, frame=0) 
         #compare = atomsel(atomselect)
         #set reference [atomselect $mol "protein" frame 0]
         # the frame being compared 
         #set compare [atomselect $mol "protein"]
         for frame in range(Trajectory.num_frames(self)):
             #protein  = atomsel(selection="protein", molid=molid, frame=frame) 
              # the frame being compared 
             
             compare = atomsel(selection=atomselect, molid=Trajectory.molID, frame=frame) 
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
         reference = atomsel(selection=atomselect, molid=Trajectory.molID, frame=0)
         reference2 = atomsel(selection=references1, molid=Trajectory.molID, frame=0) 
         #compare = atomsel(atomselect)
         #set reference [atomselect $mol "protein" frame 0]
         # the frame being compared 
         #set compare [atomselect $mol "protein"]
         for frame in range(Trajectory.num_frames(self)):
             #protein  = atomsel(selection="protein", molid=molid, frame=frame) 
              # the frame being compared 
             compare = atomsel(selection=atomselect, molid=Trajectory.molID, frame=frame) 
             
             compare2 = atomsel(selection=references1, molid=Trajectory.molID, frame=frame) 
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
    def rmsf_time(self,atomselect):
        rmsd_array=[]
        
        
        # use frame 0 for the reference
        reference = atomsel(selection=atomselect, molid=Trajectory.molID, frame=0) 
        #reference1 = atomsel(selection="protein", molid=Trajectory.molID, frame=0) 
        #num residues
        num_residues=len(pd.factorize(reference.resid)[1])
        rmsf = np.zeros(num_residues)
        
        compare2 = atomsel(atomselect)
        #set reference [atomselect $mol "protein" frame 0]
        # the frame being compared 
        #set compare [atomselect $mol "protein"]
        for frame in range(Trajectory.num_frames(self)):
             #protein  = atomsel(selection="protein", molid=molid, frame=frame) 
              # the frame being compared 
        
             compare = atomsel(selection=atomselect, molid=Trajectory.molID, frame=frame) 
             #set trans_mat [measure fit $compare $reference]
             trans_mat=atomsel.fit(compare,reference)
             # do the alignment
             compare.move(trans_mat)
             #$compare move $trans_mat
             #compute the RMSD
             #set rmsd [measure rmsd $compare $reference]
             rmsf+=compare2.rmsf(frame)
        newList = [x /Trajectory.num_frames(self)  for x in rmsf]
        #newList = map(lambda rmsf: rmsf/int(max(rmsf)), rmsf)
        #rmsf[:]=[rmsf / int(max(rmsf)) for x in rmsf]
        return newList
   
    def rmsf_time2(self,atomselect):
            
            
            #Leftraru2018
            
            # use frame 0 for the reference
            reference = atomsel(selection=atomselect, molid=Trajectory.molID, frame=0) 
            #reference1 = atomsel(selection="protein", molid=Trajectory.molID, frame=0) 
            #num residues
            num_residues=len(pd.factorize(reference.resid)[1])
            print (num_residues)
            rmsf = np.zeros(num_residues)
            
            #compare2 = atomsel(atomselect)
            #set reference [atomselect $mol "protein" frame 0]
            # the frame being compared 
            #set compare [atomselect $mol "protein"]
            mask = vmdnumpy.atomselect(molid=Trajectory.molID, frame=0,selection=atomselect)
            ref = np.compress(mask, vmdnumpy.timestep(Trajectory.molID, 0), axis=0)
            
            for frame in range(Trajectory.num_frames(self)):
                 #protein  = atomsel(selection="protein", molid=molid, frame=frame) 
                  # the frame being compared 
                 frame = np.compress(mask, vmdnumpy.timestep(Trajectory.molID, frame), axis=0)
                 rmsf += np.sqrt(np.sum((frame-ref)**2, axis=1))
                 #compare = atomsel(selection=atomselect, molid=Trajectory.molID, frame=frame) 
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
        
    def SASA(self,atomselect):
        sasa_vector=[]
        for frame in range(Trajectory.num_frames(self)):
                 lig_sel  = atomsel(selection="resname LEU ILE VAL ALA PHE TRP MET", molid=Trajectory.molID, frame=frame) 
                 big_sel  = atomsel(selection=atomselect, molid=Trajectory.molID, frame=frame) 
                 #>>> big_sel = atomsel('protein or resname LIG')
                 #>>> lig_sel = atomsel('resname LIG')
                 ligand_in_protein_sasa = big_sel.sasa(srad=1.4, restrict=lig_sel)
                 #ligand_in_protein_sasa = big_sel.sasa(srad=1.4)
                 
                 #ligand_alone_sasa= lig_sel.sasa(srad=1.4, points=True)
                 if (frame % 100)==0:
                     print("frame "+str(frame)+"  " +str(ligand_in_protein_sasa ))
                 sasa_vector.append(ligand_in_protein_sasa)
        return sasa_vector               
    def Dewetting(self,sel1,sel2, first, last):
        dew_vector=[]
        radio_vect=Trajectory.average_radius_of_gyration(self, sel1)
        print (radio_vect)

        for frame in range(first,last):
                 protein  = atomsel(selection=sel1, molid=Trajectory.molID, frame=frame) 
                 query_name_p=sel2+" and ( x>= "+str(protein.center()[0]-radio_vect)+" and x<= "+str(protein.center()[0]+radio_vect)+")"
                 radio_giro=protein.rgyr(protein.mass)
                 print (radio_giro)
                 print (query_name_p)
                 #print (protein.mass)
                 print (protein.minmax())
                 zmax=protein.minmax()[0][2]
                 
                 name_p  = atomsel(selection=query_name_p, molid=Trajectory.molID, frame=frame)
                 
                 
                 print (name_p.minmax())
                 zmin=name_p.minmax()[1][2]
                 print (zmin)
                 x=protein.center(protein.mass)[0]
                 y=protein.center(protein.mass)[1]
                 z=protein.center(protein.mass)[2]
                 cilinder_water= "water and  same residue as (((x-"+str(x)+")^2" +"+(y-"+str(y)+")^2 <="+str(radio_giro)+"^2)and z >="+str(zmin)+" and z <="+str(zmax)+")"
                 print (cilinder_water)
                 agua  = atomsel(selection=cilinder_water, molid=Trajectory.molID, frame=frame) 
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
