

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 15:33:39 2019

@author: eniac
"""
from numpy import savetxt
import subprocess
from vmd import molecule,atomsel,vmdnumpy,animate
import numpy as np
import pandas as pd
import statistics
from sklearn import preprocessing
from sklearn.linear_model import LinearRegression
class Trajectory:
    
    
    def __init__(self,dcd,psf,first=0,last=-1,stride=1,waitfor=-1):
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
        return molecule.numframes(self.molID)
    
       
    def get_molid(self):
        return self.molID
    
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
        for frame in range(Trajectory.num_frames(self)):
            sel1 = atomsel(selection=atomselect1, molid=self.molID, frame=frame) 
            sel1.moveby(-1*np.array((sel1.center(sel1.mass))))
        molecule.write(self.molID,"dcd",new_dcd_path,stride)
        
    
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
    
        return ocurrencias_vector/(last-first)
        
    
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
        ocurrence=np.transpose(args)
        value = 0
        mean_residues=[]
        std_residues=[]
        for x in ocurrence:
            
            std_residues.append(np.std(x))
            mean_residues.append(np.mean(x))
            #0.42:0.1805
            print (np.sum(x[0]))
        return [mean_residues,std_residues]
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
            atomsel("protein and resid "+str(x)).beta=residue_occupancy[x-1]
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
        
    def SASA(self,atomselect):
        sasa_vector=[]
        for frame in range(Trajectory.num_frames(self)):
                 lig_sel  = atomsel(selection="resname LEU ILE VAL ALA PHE TRP MET", molid=self.molID, frame=frame) 
                 big_sel  = atomsel(selection=atomselect, molid=self.molID, frame=frame) 
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