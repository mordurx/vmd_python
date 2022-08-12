#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 15:33:39 2019

@author: eniac
"""
import os
import glob
import natsort
import collections
import random
from numpy import savetxt
import subprocess
from vmd import molecule,atomsel,vmdnumpy,animate,trans
import numpy as np
import pandas as pd
import statistics
import re
from sklearn import preprocessing
import vmd
from sklearn.linear_model import LinearRegression
from pathlib2 import Path
class Trajectory:
    def __init__(self,dcd=None,psf=None,first=0,last=-1,stride=1,waitfor=-1,pdb=None,traj_format="dcd",top_format="psf"):
        
        # if type_format=="GROMACS":
        #     traj_format="xtc"
        #     top_format="gro"
        # if type_format=="NAMD":
        #     traj_format="dcd"
        #     top_format="psf"
        # if type_format=="test":
        #     traj_format="dcd"
        #     top_format="gro"
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
    @classmethod    
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
        return self.molID
    
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
    def tcl_volmap(self,file_path, maptypes, sel, res, output):
        #calcula el volmap density promedio y entrega un dx file.
        #sel1 = atomsel(selection=sel, molid=self.molID)
        vmd.evaltcl("source "+file_path)
        print("source "+file_path)
        print("MapVolVmd"+" "+str(self.molID)+" "+str(maptypes)+" "+f'"{sel}"'+" "+str(res)+" "+str(output))
        result=vmd.evaltcl("MapVolVmd"+" "+str(self.molID)+" "+str(maptypes)+" "+f'"{sel}"'+" "+str(res)+" "+str(output))
        #return result
    @staticmethod
    def pair_interaction_energies(energies_file,outputfolder, cutoff=0,request=['ELECT','VDW']):
        TS=[]
        ELECT=[]
        VDW=[]
        BOND=[]
        ANGLE=[]
        DIHED=[]
        IMPRP=[]
        KINEC=[]
        TOTAL=[]
        phrase='ENERGY:'
        for line in open(energies_file).readlines():
            if re.match('ETITLE:', line):
                head=re.split('\W+',line)
                #print(head)
            if re.match(phrase, line):            
                line_array=re.split(r'/^[+-]?\d+(\.\d+)?$/',line)
                line_array=line_array[0].split()
                TS.append(float(line_array[1]))
                BOND.append(float(line_array[2]))
                ANGLE.append(float(line_array[3]))
                DIHED.append(float(line_array[4]))
                IMPRP.append(float(line_array[5]))
                ELECT.append(float(line_array[6]))
                VDW.append(float(line_array[7]))
                KINEC.append(float(line_array[10]))
                TOTAL.append(float(line_array[11]))
        #zipbObj = zip([head[6],head[7],head[1]], [ELECT,VDW,TS])
        zipbObj = zip([head[1],head[2],head[3],head[4],head[5],head[6],head[7],head[10],head[11]], [TS,BOND,ANGLE,DIHED,IMPRP,ELECT,VDW,KINEC,TOTAL])
        dictOfWords = dict(zipbObj)
        
        size=len(dictOfWords["TS"][cutoff:])
        index=np.linspace(0,size, num=size)
        x=np.asarray(dictOfWords['ELECT'][cutoff:])
        print ('ELECT',np.mean(x))
        x=np.asarray(dictOfWords['VDW'][cutoff:])
        print ('VDW',np.mean(x))
        x=np.asarray(dictOfWords['BOND'][cutoff:])
        print ('BOND',np.mean(x))
        x=np.asarray(dictOfWords['ANGLE'][cutoff:])
        print ('ANGLE',np.mean(x))
        x=np.asarray(dictOfWords['DIHED'][cutoff:])
        print ('DIHED',np.mean(x))
        x=np.asarray(dictOfWords['IMPRP'][cutoff:])
        print ('IMPRP',np.mean(x))
        x=np.asarray(dictOfWords['KINETIC'][cutoff:])
        print ('KINETIC',np.mean(x))
        x=np.asarray(dictOfWords['TOTAL'][cutoff:])
        print ('TOTAL',np.mean(x))
        #return file txt
        if 'ELECT' in request:
            df=pd.DataFrame.from_dict(dictOfWords['ELECT'][cutoff:])
            df.to_csv(outputfolder+"/ELECT.txt",index=False,sep='\t',header=False)
        if 'VDW' in request:
            df=pd.DataFrame.from_dict(dictOfWords['VDW'][cutoff:])
            df.to_csv(outputfolder+"/VDW.txt",index=False,sep='\t',header=False)
        if 'BOND' in request:
            df=pd.DataFrame.from_dict(dictOfWords['BOND'][cutoff:])
            df.to_csv(outputfolder+"/BOND.txt",index=False,sep='\t',header=False)
        if 'ANGLE' in request:
            df=pd.DataFrame.from_dict(dictOfWords['ANGLE'][cutoff:])
            df.to_csv(outputfolder+"/ANGLE.txt",index=False,sep='\t',header=False)
        if 'DIHED' in request:
            df=pd.DataFrame.from_dict(dictOfWords['DIHED'][cutoff:])
            df.to_csv(outputfolder+"/DIHED.txt",index=False,sep='\t',header=False)    
        if 'IMPRP' in request:
            df=pd.DataFrame.from_dict(dictOfWords['IMPRP'][cutoff:])
            df.to_csv(outputfolder+"/IMPRP.txt",index=False,sep='\t',header=False)
        if 'KINETIC' in request:
            df=pd.DataFrame.from_dict(dictOfWords['KINETIC'][cutoff:])
            df.to_csv(outputfolder+"/KINETIC.txt",index=False,sep='\t',header=False)
        if 'TOTAL' in request:
            df=pd.DataFrame.from_dict(dictOfWords['TOTAL'][cutoff:])
            df.to_csv(outputfolder+"/TOTAL.txt",index=False,sep='\t',header=False)
    
        return dictOfWords   
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
            
    def distance_center_massXY(self,atomselect1,atomselect2):
        distance_mass_weight=[]
        for frame in range(Trajectory.num_frames(self)):
            #protein  = atomsel(selection="protein", molid=molid, frame=frame) 
        
            sel1 =atomsel(selection=atomselect1, molid=self.molID, frame=frame) 
            sel2 =atomsel(selection=atomselect2, molid=self.molID, frame=frame)
            
            
            sel1_x = np.array((sel1.center(sel1.mass)[0]))
            sel2_x= np.array((sel2.center(sel2.mass)[0]))
            
            sel1_y = np.array((sel1.center(sel1.mass)[1]))
            sel2_y= np.array((sel2.center(sel2.mass)[1]))
            
            
            dist = np.sqrt((sel2_x - sel1_x)**2 + (sel2_y - sel1_y)**2)  
            distance_mass_weight.append(dist)
        return distance_mass_weight        
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
    def contact_tamra(self,atomselect1,atomselect2,cutoff):
        for frame in range(Trajectory.num_frames(self)):
             sel  = atomsel(selection=atomselect1, molid=self.molID, frame=frame)
             receptor  = atomsel(selection=atomselect2, molid=self.molID, frame=frame)
             contact=sel.contacts(receptor,cutoff)
             if len(contact[1])==0:
                 print("vacioo")
                 continue
             unique_contact=list(set(contact[1]))
             unique_index = ' '.join([str(elem) for elem in unique_contact])
             unique_resid = atomsel(selection="index "+unique_index, molid=self.molID, frame=frame)
             unique_contact=list(set(unique_resid.resid))
             aminoacid = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M','HSE': 'H','TMR':'TMR'}
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

    def del_frame_not_contact(self,ligand,receptor,cutoff,output,output_gro,sel):
        """elimina los frame donde no hubo contacto
        """
        sel  = atomsel(selection=sel, molid=self.molID)
        delete_index=False
        range_frame=range(molecule.get_frame(self.molID))
        frame=1
        while frame <= molecule.get_frame(self.molID):
                 print(molecule.get_frame(self.molID))
                 print(frame)
                 #if delete_index==True:
                     #frame=frame-1
                 sel_1  = atomsel(selection=ligand, molid=self.molID, frame=frame) 
                 sel_2  = atomsel(selection=receptor, molid=self.molID, frame=frame) 
                 lig_contact = sel_1.contacts(sel_2,cutoff=cutoff)
                 if len(lig_contact[1])==0 or len(lig_contact[0])==0:
                     print ("delete",frame)
                     molecule.delframe(self.molID,first=frame,last=frame,stride=0)
                     
                 else:
                    frame=frame+1
        molecule.write(self.molID, "gro", output_gro,first=molecule.get_frame(self.molID),last=-1,stride=1,selection=sel)    
        molecule.write(self.molID,'trr',output,first=1,last=-1,stride=1,selection=sel)
        molecule.cancel             

                 
    def contact_map_protein(self,ligand,receptor,cutoff):
        dict_resid_count={}
        summary_hbound={}
        aminoacid = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M','HSE': 'H','TMR':'TMR'}
        # cuantos contactos tiene un ligando a cierto receptor by cutoff
        for frame in range(Trajectory.num_frames(self)):
                 sel_1  = atomsel(selection=ligand, molid=self.molID, frame=frame) 
                 sel_2  = atomsel(selection=receptor, molid=self.molID, frame=frame) 
                 lig_contact = sel_1.contacts(sel_2,cutoff=cutoff)
                 #H_bond_vector.append(lig_sel_sasa)
                 #number_hbond[frame]=0
                 dict_aux={}
                 for data in range(len(lig_contact[1])):
                     sel_A  = atomsel(selection="index "+str(lig_contact[0][data]), molid=self.molID, frame=frame)
                     sel_B  = atomsel(selection="index "+str(lig_contact[1][data]), molid=self.molID, frame=frame)
                     #sel_H  = atomsel(selection="index "+str(lig_contact[2][data]), molid=self.molID, frame=frame)
                     
                     ligand_atom=str(sel_A.resid[0])+""+str(aminoacid[sel_A.resname[0]])+"--"+str(sel_B.resid[0])+""+str(aminoacid[sel_B.resname[0]])
                     if ligand_atom not in dict_aux.keys():
                            dict_aux[ligand_atom]=1
                     else:
                         if ligand_atom not in dict_aux.keys():
                            dict_aux[ligand_atom]=1
                         
                 dict_aux = collections.Counter(dict_aux)
                 summary_hbound = collections.Counter(summary_hbound)        
                 summary_hbound=summary_hbound+dict_aux    
        
        for key in summary_hbound:
            summary_hbound[key]= (summary_hbound[key]/ Trajectory.num_frames(self)*100)            
        
        sorted_dict = dict(collections.OrderedDict(natsort.natsorted(summary_hbound.items())))             
        return sorted_dict             
    def contact_time(self,ligand,protein,cutoff):
        dict_resid_count={}
        # cuantos contactos tiene un ligando a cierto receptor by cutoff
        
        for frame in range(Trajectory.num_frames(self)):
             sel  = atomsel(selection=ligand, molid=self.molID, frame=frame)
             receptor  = atomsel(selection=protein, molid=self.molID, frame=frame)
             contact=sel.contacts(receptor,cutoff)
             #la primera lista index del lig
             #segunda lista atom que hacen contacto
             aminoacid = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M','HSE': 'H','TMR':'TMR'}
             #dejar solo valores unicos, no nos intereza el analisis de enlaces
             if len(contact[1])==0:
                 #print("vacioo")
                 continue
             unique_contact=list(set(contact[1]))
             unique_index = ' '.join([str(elem) for elem in unique_contact])
             unique_resid = atomsel(selection="index "+unique_index, molid=self.molID, frame=frame)
             unique_contact=list(set(unique_resid.resid))
           
             for atom in unique_contact:
                 resid = atomsel(selection="resid "+str(atom) +" and "+protein, molid=self.molID, frame=frame)
                 resid=list(set(resid.resname))
                 get_one_letter_amino=aminoacid[resid[0]]
                 if str(atom)+get_one_letter_amino in dict_resid_count.keys():
                     #print("existe residuo")
                     dict_resid_count[str(atom)+get_one_letter_amino]=dict_resid_count[str(atom)+get_one_letter_amino]+1
                 else:
                     dict_resid_count[str(atom)+get_one_letter_amino]=1
        for i in dict_resid_count.keys():
            dict_resid_count[i]=dict_resid_count[i]/Trajectory.num_frames(self)
        
        sorted_dict = dict(collections.OrderedDict(natsort.natsorted(dict_resid_count.items())))             
        return sorted_dict     
        
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
    def max_insert_residue(self,atomselect1,atomselect2,first,last):
        """
        entrega el maximo de residuos insertados
        """
        sel1 = atomsel(selection=atomselect1, molid=self.molID, frame=0) 
        frame_first_insertion=0
        insert_resid_by_frame={}     
        for frame in range(first,last):
             sel1  = atomsel(selection=atomselect1, molid=self.molID, frame=frame) 
             sel2 =atomsel(selection=atomselect2, molid=self.molID, frame=frame)
             resid_center=sel1.centerperresidue()
            
                 
             resid_z=list(map(lambda x: x, resid_center))
             nameP_z=list(map(lambda x: x[2], sel2.minmax()))
             membrane=atomsel(selection="resname POPC POPG", molid=self.molID, frame=frame)
             mem_mass_center_Z=membrane.center(membrane.mass)[2]
             
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
             cont=0
             for number_resid in range(len(resid_z)):
                 x_residue=resid_z[number_resid][0]
                 y_down_ajuste = p_down[0]*x_residue + p_down[1]
                 y_up_ajuste = p_up[0]*x_residue + p_up[1]
                 
                 if resid_z[number_resid][2]<y_up_ajuste and resid_z[number_resid][2]>y_down_ajuste:
                     cont=cont+1
             insert_resid_by_frame[str(frame)]=cont
             #print (ocurrencias_temp)
             #ocurrencias_temp= np.zeros(len(resid_center))    
          
    
        return max(insert_resid_by_frame.values())
         
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
    
    def rmsd_time(self,atomselect,first=0,last=-1):
         rmsd_array=[]
         if last==-1:
            last=Trajectory.num_frames(self)-1
         if first==-1:
            first=Trajectory.num_frames(self)
         else:
             last=last+1   
         # use frame 0 for the reference
         reference = atomsel(selection=atomselect, molid=self.molID, frame=0) 
         #compare = atomsel(atomselect)
         #set reference [atomselect $mol "protein" frame 0]
         # the frame being compared 
         #set compare [atomselect $mol "protein"]
         for frame in range(first,last):
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
    def colour_pdb_by_column(sel_res_values:dict, query,molid:int,output:str,colunm='beta'):
        """
        metodo que colorea por beta factor
        residue_occupancy: vector con valores a colocar en beta columm
        query: atomselecion "protein"
        output "pbd salida"
        molid: id de la proteina  
        """
        aminoacid = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M','HSE': 'H','TMR':'TMR'}
        protein = atomsel(query)
        if colunm=='beta':
            protein.beta=0
        elif colunm=='occupancy':
            protein.occupancy=0                
        keys=list(sel_res_values.keys())
        for key, value in sel_res_values.items():
            resid_letter=key[-1]
            resid_index=key[:-1]
            for key_A, value_A in aminoacid.items():
                if resid_letter== value_A:
                    if colunm=='beta':
                        atomsel(query+" and resid "+str(resid_index)+" and resname "+key_A).beta=value
                    elif colunm=='occupancy':
                        atomsel(query+" and resid "+str(resid_index)+" and resname "+key_A).occupancy=value
                    
                        
        
                 
        #for x in protein.resid:
        #    atomsel(query+" and resid "+str(x)).beta=residue_occupancy[x-1]
        #    print(x,residue_occupancy[x-1])
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
                       
        self.to_excel(output=output,resume=summary_hbound,numHbound=number_hbond,sheet=sheet, mode=mode)                             
    def wrap_CG(self,wrapselect,output_traj_trr,output_gro,it):
        #con evaltcl puedo acceder a la consola tk
        vmd.evaltcl('package require pbctools')
        for frame in range(Trajectory.num_frames(self)):
                i=1
                for i in range(it):
                    vmd.evaltcl('pbc wrap -centersel "'+wrapselect+'" -center com  -first '+str(frame)+' -last '+str(frame)+' -molid '+str(self.molID))
                #vmd.evaltcl('set all2 [atomselect '+str(molID)+ " all " 'frame '+str(frame)+']')
                sel_all= atomsel(selection="all", molid=self.molID, frame=frame) 
                sel_all.moveby(-1*np.array((sel_all.center())))
                #vmd.evaltcl(sel_all+'moveby [vecinvert [measure center'$all']]')
                #$all moveby [vecinvert [measure center $all]]
                print('pbc wrap -centersel "not resid 1 to 41 and name BB SC1 SC2 SC3 SC4 SC5" -center com  -first '+str(frame)+' -last '+str(frame)+' -molid '+str(self.molID))
                if frame==Trajectory.num_frames(self)-1:
                        molecule.write(self.molID, "gro", output_gro,first=frame,last=frame,stride=1)

                    #sel1 = atomsel(selection=atomselect1, molid=self.molID, frame=frame)
        molecule.write(self.molID,'trr',output_traj_trr,first=1,last=-1,stride=1)
        molecule.cancel
    @staticmethod    
    def pdb_to_traj(folder_pdb_input,name_output,folder_pdb_output,format_traj_output='pdb',format_pdb_input='pdb',psf_gro_file=None):
        #convierte una carpeta de pdb o gro a una trajectoria, la trayectoria puede ser de pdb o gro o dcd
        path_output=os.path.join(folder_pdb_output, name_output)
        path_input=os.path.join(folder_pdb_input)
        pattern_to_found=os.path.join(path_input,'*.'+format_pdb_input)
        array_pdb=glob.glob(pattern_to_found)
        
        molID=molecule.new('pdb')
        for pdb in array_pdb:
            molecule.read(molID,filetype=format_pdb_input,filename=pdb)
        molecule.write(molid =molID,filetype=format_traj_output,filename=path_output) 
        molecule.delete(molID)
        
    
    @staticmethod    
    def gro_to_pdb(gro,output,atomselect_rec="protein",atomselect_lig="not protein",fixchain=False):
        molID=molecule.load('gro',gro) 
        receptor  = atomsel(selection=atomselect_rec, molid=molID) 
        lig  = atomsel(selection=atomselect_lig, molid=molID) 
        all_complex  = atomsel(selection="all", molid=molID) 
        if fixchain==True:
                receptor.chain='R'
                lig.chain='L'
        all_complex.write('pdb',output)  
        molecule.delete(molID)        
    def traj_to_pdb(self,outputfolder,output_format,name_output,atomselect_lig,atomselect_rec,fixchain=False):
        #a partir de la trajectoria genera pdb
        pdb_vector=[]
        for frame in range(Trajectory.num_frames(self)):
            path_output=os.path.join(outputfolder, name_output+str(frame)+"."+output_format)
        
            all_complex  = atomsel(selection="all", molid=self.molID,frame=frame)
            if fixchain==True:
                receptor  = atomsel(selection=atomselect_rec, molid=self.molID ,frame=frame) 
                lig  = atomsel(selection=atomselect_lig, molid=self.molID,frame=frame) 
                receptor.chain='R'
                lig.chain='L'
                
            all_complex.write(output_format,path_output)
            pdb_vector.append(path_output)
        return pdb_vector
    @staticmethod 
    def reformat_pdb(pdb):
        path = Path(pdb)
        text = path.read_text()
        text = text.replace('HSE', 'HIS')
        text = text.replace('HSP', 'HIS')
        text = text.replace('CD  ILE', 'CD1 ILE')
        text = text.replace('OT1 GLN', 'O   GLN')
        text = text.replace('OT2 GLN', 'OXT GLN')
        text = text.replace('OT1 ASP', 'O   ASP')
        text = text.replace('OT2 ASP', 'OXT ASP')
        path.write_text(text)
    def rmsd_matrix(self,atomselect,file_name,ref=None):
        """genera una matriz de rmds"""
        rmsd_array=[]
        matrix=[]
        # use frame 0 for the reference
        
        if ref==None:
            reference = atomsel(selection=atomselect, molid=self.molID, frame=0)
        else:
            reference = atomsel(selection=atomselect, molid=ref)    
            
        for frame_row in range(Trajectory.num_frames(self)):
            rmsd_array=[]
            
            for frame_col in range(frame_row,Trajectory.num_frames(self)):
                
                reference2 = atomsel(selection=atomselect, molid=self.molID, frame=frame_col) 
                #set trans_mat [measure fit $compare $reference]
                trans_mat=atomsel.fit(reference2,reference)
                # do the alignment
                #compare.move(trans_mat)
                reference2.move(trans_mat)
                #$compare move $trans_mat
                # compute the RMSD
                #set rmsd [measure rmsd $compare $reference]
                matrix.append(int(atomsel.rmsd(reference2,reference)*1000))
            #matrix.append(rmsd_array)
            reference=atomsel(selection=atomselect, molid=self.molID, frame=frame_row)     
        df = pd.DataFrame(matrix)
        df.to_csv(file_name, sep='\t', index = None, header=False)
        return df        