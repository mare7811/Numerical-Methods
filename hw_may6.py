# -*- coding: utf-8 -*-
"""
Miguel Mares
Numerical Methods
Homework
Due 5/8/2020
"""
#This program attempts to implement the ADI method to solve for the heat 
#distribution of a given plate and boundary conditions. 
#I could not figure out the error in my code.

import numpy as np
import math


#The following functions calculate the sums for the specified entries. 
#For example, l_sum() sums up the appropriate entries for L[i][k].
def l_sum(i,k,L,U):
    l_sum = 0
    for j in range(0,k):
        l_sum += L[i][j]*U[j][k]
    return l_sum

def u_sum(k,j,L,U):
    u_sum = 0
    for i in range(0,k):
        u_sum += L[k][i]*U[i][j]
    return u_sum

def z_sum(i,L,z):
    z_sum = 0
    for k in range(0,i):
        z_sum += L[i][k]*z[k][0]
    return z_sum

def x_sum(i,x,U):
    x_sum = 0;
    for k in range(N-i+1,N):
        x_sum += U[N-i][k]*x[k][0]
    return x_sum


#This function solves the system using the global A matrix and global b vector
#The b vector is the right hand side vector (either S or R) 
def LU():
    #print(A,'\n',b)
    L = np.zeros((N,N))
    for i in range(0,N):
        L[i][0] = A[i][0]
        
#Creates the U matrix and sets known entries.
    U = np.zeros((N,N))
    for i in range(0,N):
        U[0][i] = A[0][i]/L[0][0]

#This section calculates the unkown L and U entries using the known values and
#sum functions, then prints the L and U matrices. 
    for k in range(1,N):
        for z in range(0,N):
            L[z][k] = A[z][k] - l_sum(z,k,L,U)
            if(L[k][k] == 0):
                U[k][z] = 0
            else:
                U[k][z] = (A[k][z] - u_sum(k,z,L,U))/L[k][k]
    #print('L = \n',L,'\n')
    #print('U = \n',U,'\n')
                
#This section calculates the z matrix using the sum function and b matrix.
    z = np.zeros((N,1))
    for i in range(0,N):
        z[i][0] = (b[i][0] - z_sum(i,L,z))/L[i][i]   
#print('z = \n',z,'\n')

#This section calculates the x values (solution) using the sum function and
#z matrix and prints the result.
    x = np.zeros((N,1))
    x[N-1][0]=z[N-1][0]
    for i in range(2,N+1):
        x[N-i][0] = z[N-i][0] - x_sum(i,x,U)
    #print('x = \n',x,'\n')
    return x

#This function updates the u matrix "plate" with the solved tridiagonal system
def update():
    for i in range(1,Nx-1):
        for j in range(1,Nx-1):
            u[i][j] = s[i][j]
##############################################################################
#plate conditions
a = 0 
b = 1
h = 1/3 #step size
d = 1/3 #time step size
p = 1   #equivalent to rho in algorithm of A matrix

Nx = math.floor((b-a)/h)+1  #Number of points on plate
N = Nx     #for LU decomp function purposes


u = np.zeros((Nx,Nx))   #Array of plate vaues
x = np.zeros(Nx)        #array to hold discretized space values
y = np.zeros(Nx)        #This is not necessary but representative of discretized space

#Fill arrays with discretized space values
for i in range(0,Nx):
    x[i] = i*h    
for i in range(0,Nx):
    y[i] = i*h

#Set boundary conditions
for i in range(1,Nx):
    u[0][i] = 10
    u[Nx-1][i] = 10
    u[i][Nx-1] = 10
print('initial u \n',u)
        
#Create A matrix
A = np.zeros((Nx,Nx))
for i in range(0,Nx):
    A[i][i] = 2+(1/p)
    if(i == Nx-1):
        break
    A[i+1][i] = -1
    A[i][i+1] = -1
print('A = \n',A)

maxit = 3   #max number of iterations
S = np.zeros((Nx,1))    #create array for S and R vectors
R = np.zeros((Nx,1))
#Boundary conditions
S[0][0] = 10
S[Nx-1][0] = 10
R[0][0]= 0
R[Nx-1][0] = 10

iteration=0 #iteration counter
for n in range(0,maxit):
    print('n:',iteration)
    s = np.zeros((Nx,1))
    for i in range(0,Nx-1): #This loop sweeps through for each S
        for j in range(0,Nx-1): #This loop fills R.H.S vector for S
            S[j][0] = u[i+1][j] - (2-h**2/d)*u[i][j] + u[i-1][j]
        b = S
        solution = LU()
        s = np.concatenate((s,solution),axis=1)
    s = np.delete(s,0,1)
    update()
    print(u,'\n')
    
        
    iteration += 1
    print('n:',iteration)
    s = np.zeros((1,Nx))
    for j in range(0,Nx-1):
        for i in range(0,Nx-1):
            R[i][0] = u[i][j+1] - (2-h**2/d)*u[i][j] + u[i][j-1]
        b = R
        solution = LU()
        solution = np.resize(solution,(1,Nx))
        s = np.concatenate((s,solution),axis=0)
    s = np.delete(s,0,0)
    update()
    print(u,'\n')
    iteration += 1


print("I am not sure where the error is in my code but I am assuming that it",
      " has something to do with how I sweep through the columns and rows",
      " and subsequently update the u matrix.")
    
