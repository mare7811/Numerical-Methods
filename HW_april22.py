# -*- coding: utf-8 -*-
"""
Miguel Mares
Numerical Methods
Homework
Due 4/22/2020
"""
import numpy as np
import math 
import pandas as pd

h = 1
d = .1

#initial conditions
u_0t = 0
u_3t = 0
c = 1

#Creating grids 
Nx = math.floor((3-0)/h)
Nt = math.floor((.5-0)/d)
x = np.zeros(Nx+1)
t = np.zeros(Nt+1)

#x grid
for i in range(0,Nx+1):
    x[i] = i*h

#boundary conditions
u = np.zeros((Nx+1,Nt+1))
for i in range(0,Nt):
    u[0][i] = u_0t
    u[1][i] = u_3t
    
#forward difference algorithm
for i in range(0,Nx):
    u[i][0] = x[i]*(3-x[i])
    for j in range(0,Nt):
        u[i][j+1] = (c**2)*((d**2)/(h**2))*(u[i+1][j]-2*u[i][j]+u[i-1][j])+2*u[i][j]-u[i][j-1]
    
    
T = pd.DataFrame(u)
T.columns = ['0','.1','.2','.3','.4','.5']
print(T)
print('\nThe solution is stable')
