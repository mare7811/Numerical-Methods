# -*- coding: utf-8 -*-
"""
Miguel Mares
Numerical Methods
Final Exam problem 2
"""
#This program uses the power method to obtain the smallest 
#eigenvalue of the given matrix and using the given initial
#guess for x.

import numpy as np
import math


a = 0
b = 3
h = 1
M = 5
nx = math.floor((b-a)/h) - 1


a = -1/h**2
def b(x):
    b = (2+(h**2)*(x**2))/h**2
    return b
c = -1/h**2

#matrix, A, from second order ODE, finite difference equation
A = np.zeros((nx,nx))
A[0][0] = b(1)
A[0][1] = c
A[nx-1][nx-2] = a
A[nx-1][nx-1] = b(1)
for i in range(1,nx-1):
    A[i][i-1] = a
    A[i][i] = b(0)
    A[i][i+1] = c

print('A = \n',A)


constant = 1/np.sqrt(2)
x = [1*constant,1*constant]
for i in range(0,M):
    x = np.dot(A,x)
    x_norm = np.linalg.norm(x)
    x = x/x_norm

Ax = np.dot(A,x)
num = np.dot(Ax,x)
denom = np.dot(x,x)
l = num/denom

print('iterations = ',M)
print('lambda = ',l)
print('\nThis is the smaller eignenvalue')