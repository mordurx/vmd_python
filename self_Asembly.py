#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from trajectory import Trajectory
import numpy as np
from PdbGenerator import PdbGenerator
pdb='test/step4_lipid.pdb'
lipid='test/lipid_one.pdb'
tcl='tcl_script/01-coarse-grain.tcl'
psf='test/step4_lipid.psf'
test="test/generate_ramdom_pose.pdb"
snxpdb="test/AA-snx_482.pdb"
snxpsf="test/AA-snx_482.psf"
#llamo a mi clase tratectoria rep0
box_size = [(-520,-520,-520),(520,520,520)]
PdbGenerator1=PdbGenerator(snxpdb,"protein",'test/generate_ramdom_pose2.pdb',snxpsf)
ratio=15
PdbGenerator1.generate_random_poses_from_element(box_size,ratio)
PdbGenerator1.close()

PdbGenerator2=PdbGenerator('test/generate_ramdom_pose2.pdb',"segname TOX",'test/generate_ramdom_pose3.pdb')
PdbGenerator2.sort_random_sel(box_size)

#print(outList)
#x=PdbGenerator1.call_tcl("play "+tcl,"frame")
#print(x)

