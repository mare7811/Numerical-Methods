# -*- coding: utf-8 -*-
"""
Miguel Mares
Engr 428
Homework 
Due 2/7/2020
"""
#This program solves the given general matrix using LU decomposition.

import numpy as np

N = 3       #For a nxn matrix, N=n

#Creates the given A matrix
A = np.zeros((N,N))
A[0][0] = .1
A[0][1] = .2
A[0][2] = .4
A[1][0] = 4
A[1][1] = 1
A[1][2] = -1
A[2][0] = 2
A[2][1] = 5
A[2][2] = 2
print('A = \n',A,'\n')

#Creates the given b matrix
b = np.zeros((N,1))
b[0][0] = 1.1
b[1][0] = 6
b[2][0] = 3
print('b = \n',b,'\n')

#Creates the L matrix and sets known entries.
L = np.zeros((N,N))
L[0][0] = A[0][0]
L[1][0] = A[1][0]
L[2][0] = A[2][0]

#Creates the U matrix and sets known entries.
U = np.zeros((N,N))
U[0][0] = A[0][0]/L[0][0]
U[0][1] = A[0][1]/L[0][0]
U[0][2] = A[0][2]/L[0][0]


#The following functions calculate the sums for the specified entries. 
#For example, l_sum() sums up the appropriate entries for L[i][k].
def l_sum(i,k):
    l_sum = 0
    for j in range(0,k):
        l_sum += L[i][j]*U[j][k]
    return l_sum

def u_sum(k,j):
    u_sum = 0
    for i in range(0,k):
        u_sum += L[k][i]*U[i][j]
    return u_sum

def z_sum(i):
    z_sum = 0
    for k in range(0,i):
        z_sum += L[i][k]*z[k][0]
    return z_sum

def x_sum(i):
    x_sum = 0;
    for k in range(N-i+1,N):
        x_sum += U[N-i][k]*x[k][0]
    return x_sum



#This section calculates the unkown L and U entries using the known values and
#sum functions, then prints the L and U matrices. 
for k in range(1,N):
    for z in range(0,N):
        L[z][k] = A[z][k] - l_sum(z,k)   
        U[k][z] = (A[k][z] - u_sum(k,z))/L[k][k]
print('L = \n',L,'\n')
print('U = \n',U,'\n')

#This section calculates the z matrix using the sum function and b matrix.
z = np.zeros((N,1))
for i in range(0,N):
    z[i][0] = (b[i][0] - z_sum(i))/L[i][i]   
print('z = \n',z,'\n')

#This section calculates the x values (solution) using the sum function and
#z matrix and prints the result.
x = np.zeros((N,1))
x[N-1][0]=z[N-1][0]
for i in range(2,N+1):
    x[N-i][0] = z[N-i][0] - x_sum(i)
print('x = \n',x,'\n')

#This section checks the solution by multiplying the A matrix with the x matrix
r = np.zeros((N,1))
for i in range(0,N):
    r[i][0] = A[i][0]*x[0][0] + A[i][1]*x[1][0] + A[i][2]*x[2][0]
print('solution check: \n',r)