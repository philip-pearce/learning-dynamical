#! /usr/bin/python2.7
from __future__ import division
import numpy as np

weight = np.genfromtxt('data/fortran/weight.out')
matrix = np.genfromtxt('data/fortran/matrix1.out')

if( weight.size == 1 ): nump_comp = 1
else: num_comp = len(weight)

if( matrix.size == 1 ): num_dim = 1
else: num_dim = len(matrix[0,:])

print "Loading GMM with", num_comp, "components and", num_dim, "dimensions..."

file = open('source/fitness.f90', 'r') 
lines = file.readlines()
temp = '    integer, parameter :: NDIMENSION = ' + str(num_dim) + ', NCOMPONENT = ' + str(num_comp) + '\n'
lines[6] = temp
file.close()

file = open('source/fitness.f90', 'w')
for i in range(0, len(lines)):	
	file.write(lines[i])
file.close()

