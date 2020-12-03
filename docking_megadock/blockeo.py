#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 12:28:23 2019

@author: mordurx
"""

import os
import sys
import csv

class Bloqueo(object):
    reslist = []
   
    def list_ban():
        return Bloqueo.reslist
    def __lista_resid_block(ban_residue):
        resarg = ban_residue.split(",")
        print(resarg)
        for r in resarg:
        #restrcion de residuos
            if "-" in r:
                rseq = r.split("-")
                for i in range(int(rseq[0]), int(rseq[1])+1):
                    Bloqueo.reslist.append(str(i)) #entrega el rango completo de block residues 100,101,102.....10N
            else:
                Bloqueo.reslist.append(r)    
        
    @staticmethod
    def block_resid(pdb,banlist,salida):
        with open(salida, 'w') as f:
            Bloqueo.__lista_resid_block(banlist)
            fp = open(pdb, "r")
            try:
                for l in fp.readlines():
                    #l = l.strip()
    
                    if l[0:4] != "ATOM" and l[0:6] != "HETATM":
                        print (l, file=f)
                        continue
        
               #    if l[23] != ch:
               #        print (l)
               #        continue
        
                    ll = l
                    if ll[23:27].strip() not in Bloqueo.reslist:
                        print (l, file=f)
                        continue
                
                    print (l[0:17] + "BLK " + l[21:], file=f)
            
            finally:
                fp.close()

        
            
