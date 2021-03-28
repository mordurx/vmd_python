#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from trajectory import Trajectory
import numpy as np
from PdbGenerator import PdbGenerator
pdb='test/step4_lipid.pdb'
tcl='tcl_script/01-coarse-grain.tcl'
psf='test/step4_lipid.psf'
#llamo a mi clase tratectoria rep0
box_size = [(-120,-120,-120),(120,120,120)]
PdbGenerator1=PdbGenerator(pdb,"resname POPC POPG",'test/popg_popc_random.pdb',psf)
PdbGenerator1.self_assembly(box_size)
#x=PdbGenerator1.call_tcl("play "+tcl,"frame")
#print(x)
PdbGenerator1.close()
