# -*- coding: utf-8 -*-
"""
Miguel Mares
Numerical Methods
Homework
Due 5/1/2020
"""
import numpy as np
import math 
import pandas as pd
pd.options.display.float_format='{:,.2f}'.format    #format table for two decimal places
pd.set_option('display.max_columns',500)    

#Define functions for partial derivatives
def A(j):
    return -700

def B(j):
    return -200
    
def C(i):
    return 400

def D(i):
    return -100

#Set step size, boundary limits, number of points
a = 0
b = 1
h = 1/11
N = math.floor((b-a)/h) + 1

#Set tolerance and alpha 
tolerance = 1e-6
alpha = .5

#Create matrix to hold old values and a matrix to represent plate
Y = np.zeros((N,N))
U = np.zeros((N,N))

#Corner plate values
U[0][0] = 750
U[0][N-1] = (U[0][N-2]+h*D(0)+U[1][N-1]-h*A(N-1))/2
U[N-1][0] = (U[N-2][0]+h*B(0)+U[N-1][1]-h*C(N-1))/2
U[N-1][N-1] = (U[N-2][N-1]+h*B(N-1)+U[N-1][N-2]+h*D(N-1))/2


#Boundary values
for j in range(1,N-1):
    U[0][j] = .25*(2*U[1][j] - 2*h*A(j) + U[0][j+1] + U[0][j-1])
    U[N-1][j] = .25*(2*U[N-2][j]+2*h*B(j)+U[N-1][j+1]+U[N-1][j-1])

for i in range(1,N-1):
    U[i][0] = .25*(U[i+1][0]+U[i-1][0]+2*U[i][1]-2*h*C(i))
    U[i][N-1] = .25*(U[i+1][N-1]+U[i-1][N-1]+2*U[i][N-2]+2*h*D(i))


#This section fills in the rest of the plate values
iteration = 0   #iteration counter
done  = False
while(done == False):
    for i in range(1,N-1):
        for j in range(1,N-1):
            U[i][j] = .25*(U[i+1][j]+U[i-1][j]+U[i][j+1]+U[i][j-1])
            U[i][j] = U[i][j] + alpha*(U[i][j] - U[i][j])
            if(np.abs((U[i][j]-Y[i][j])/U[i][j]) < tolerance):    #check tolerance
                done = True 
    for i in range(1,N-1):
        for j in range(1,N-1):
            Y[i][j] = U[i][j]
    if(done == False):
        iteration += 1
    
    
#Output array of values that represent plate temps (U array)    
print('\n Number of iterations:',iteration)
U = np.flip(U,0)
U = pd.DataFrame(U)
print(U.to_string(index=False,header=False))

print('The tolerance was set to .000001. The h value was set small enough that'
      ' a large set of values could be seen.')
