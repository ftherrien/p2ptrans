#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
JVC algorithm using scipy optimize. This part of the code solves the linear assingment problem and saves the distance and the minimum atom-to-atom assignment into two files, for further processing by the fortran routine. 
"""
from scipy.io import FortranFile
import numpy as np
from scipy.optimize import linear_sum_assignment
import os
ff = FortranFile('file.dat', 'r') # Specifying read-only access.
#os.system('rm *.dat')
data = ff.read_reals(dtype=np.float)
n=int(np.sqrt(len(data)))
dataReshaped = data.reshape(n, n)
Assignment=linear_sum_assignment(dataReshaped)
cost=0
for i in range(n):
   cost=cost+dataReshaped[Assignment[0][i],Assignment[1][i]]
with open('cost.csv', 'w') as f:
    f.write(str(cost))

Map=[int(Assignment[1][i])+1 for i in range(n)]
np.savetxt("map.csv", Map, delimiter=",", fmt="%d")
   
