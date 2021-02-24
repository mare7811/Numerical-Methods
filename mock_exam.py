# -*- coding: utf-8 -*-
"""
Miguel Mares
Numerical Methods
Mock Exam
4/27/2020
"""
import numpy as np
import math 
import pandas as pd


def A(j):
    #A = (U[0][j-1]-U[-2][j-1])/(2*h)
    return -700

def B(j):
    #B = (U[N-1][j-1]-U[N-3][j-1])/(2*h)
    return -200
    
def C(i):
    #C = (U[i-1][0]-U[i-1][N-2])/(2*h)
    return 400

def D(i):
    #D = (U[i-1][N-1]-U[i-1][N-3])/(2*h)
    return -100

#Set step size, boundary limits, number of intervals
a = 0
b = 1
h = 1/3
N = math.floor((b-a)/h) + 1

#Set tolerance and alpha 
tolerance = 1e-4
alpha = .5


Y = np.zeros((N,N))
U = np.zeros((N,N))

U[0][0] = 750
#U[0][0] = (U[0][1]-h*C(0)+U[1][0]-h*A(0))/2
U[0][N-1] = (U[0][N-2]+h*D(0)+U[1][N-1]-h*A(N-1))/2
U[N-1][0] = (U[N-2][0]+h*B(0)+U[N-1][1]-h*C(N-1))/2
U[N-1][N-1] = (U[N-2][N-1]+h*B(N-1)+U[N-1][N-2]+h*D(N-1))/2


for j in range(1,N-1):
    U[0][j] = .25*(2*U[1][j] - 2*h*A(j) + U[0][j+1] + U[0][j-1])
    U[N-1][j] = .25*(2*U[N-2][j]+2*h*B(j)+U[N-1][j+1]+U[N-1][j-1])

for i in range(1,N-1):
    U[i][0] = .25*(U[i+1][0]+U[i-1][0]+2*U[i][1]-2*h*C(i))
    U[i][N-1] = .25*(U[i+1][N-1]+U[i-1][N-1]+2*U[i][N-2]+2*h*D(i))


iteration = 0
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
    
print('Number of iterations:',iteration)
U = np.flip(U,0)
print(U)