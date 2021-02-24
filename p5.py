# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 18:35:21 2020

@author: admin
"""
import numpy as np
import math


h = .1 #step size
mu0 = .000954   #average mass in kilograms
T = 1000    #Tension in Newtons
L = 1   #Length of string
delta = .0005   #

#Create points along string based on step size
nx = math.floor(L/h) - 1
x = np.zeros(nx)
for i in range(1,nx+1):
    x[i-1] = i*h

#'A' matrix elements
mu = mu0+(x-L/2)*delta
a = -T/(h*h*mu)
b = -2*a
c = a


y = np.zeros((nx,1))
#choose initial guess to be y = (1,1)T
for i in range(0,nx):
    y[i][0] = 1


#matrix, A, from second order ODE, finite difference equation
A = np.zeros((nx,nx))
A[0][0] = b[0]
A[0][1] = c[0]
A[nx-1][nx-2] = a[nx-1]
A[nx-1][nx-1] = b[nx-1]
for i in range(1,nx-1):
    A[i][i-1] = a[i]
    A[i][i] = b[i]
    A[i][i+1] = c[i]

print(A)

for i in range(0,10):
    z = np.dot(A,y);
    #print(z);
    eigenval = z[1][0]

    for j in range(1,nx):
        if z[j][0] > eigenval:
            eigenval = z[j][0];

    for j in range(0,nx):
        z[j][0] = z[j][0]/eigenval;
    
    eigenval = 1/eigenval;
    y = z;

print(eigenval)