#!/usr/bin/env python3
import re, sys, os
import numpy as np

def pair_interations_energies(f,phrase):
    TS=[]
    ELECT=[]
    VDW=[]
    for line in open(f).readlines():
        if re.match('ETITLE:', line):
            head=re.split('\W+',line)
            #print(head)
        if re.match(phrase, line):            
            line_array=re.split(r'/^[+-]?\d+(\.\d+)?$/',line)
            line_array=line_array[0].split()
            ELECT.append(float(line_array[6]))
            VDW.append(float(line_array[7]))
            TS.append(float(line_array[1]))
    zipbObj = zip([head[6],head[7],head[1]], [ELECT,VDW,TS])
    dictOfWords = dict(zipbObj)
    return dictOfWords