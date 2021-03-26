#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from trajectory import Trajectory
import numpy as np

pdb='docking_megadock/step1_pdbreader.pdb'
psf='docking_megadock/step1_pdbreader.psf'
#llamo a mi clase tratectoria rep0
Trajectory1=Trajectory(pdb,psf)
Trajectory1.close()
