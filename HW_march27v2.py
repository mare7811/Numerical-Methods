# -*- coding: utf-8 -*-
"""
Miguel Mares
Numerical Methods
Homework 
Due 3/30/2020
"""

import numpy as np
import matplotlib.pyplot as plt
import math
import pandas as pd

#Defining functions
def F1(x,S1,S2):
    f1 = S2
    return f1
    
def F2(x,S1,S2):
    f2 = -(S2)**2 - S1 + np.log(x)
    return f2

#Defining step size and boundary values
h = 0.25
x_0 = 1
S1_0 = 0

#Interval limits and number of intervals
a = 1
b = 2
N = math.floor((b-a)/h)

#Defining arrays to hold values
x = np.zeros((N+1,1))
S1 = np.zeros((N+1,1))
S2 = np.zeros((N+1,1))

#Defining range of initial conditions to test
R = np.arange(-10,10,1)
g = np.zeros((R.size,1))    #Holds difference between final value of solution and final boundary value for each initial value tested


#Function that finds root using secant method
def secant(x1,x2,y1,y2):
    #Create lists and set first entries
    r = []
    y = []
    r.append(x1)
    r.append(x2)
    y.append(y1)
    y.append(y2)

    #loops through x values until the difference between two y values is too small to use
    i = 2
    while(r[i-1] != 0):
        if((y[i-1]-y[i-2]) != 0 ):
            a = (r[i-2]*y[i-1] - r[i-1]*y[i-2])/(y[i-1]-y[i-2])
            r.append(a)
            y.append(modified_euler(r[i]) - np.log(2))
        else:
            r.append(0)
        i += 1
    return(r[i-2])

#Functions that finds final value of solution using Modified Euler method
def modified_euler(S2_0):
    #Filling time array with step values and other arrays with initial values
    for i in range(0,N+1):
        x[i] = x_0 + i*h
    S1[0] = S1_0
    S2[0] = S2_0

    #Filling/calculating rest of array values using Modified Euler 
    for i in range(1,N+1):
        y_mid = S1[i-1] + (h/2)*F1(x[i-1],S1[i-1],S2[i-1])
        v_mid = S2[i-1] + (h/2)*F2(x[i-1],S1[i-1],S2[i-1])
    
        S1[i] = S1[i-1] + h*F1(x[i-1],y_mid,v_mid)
        S2[i] = S2[i-1] + h*F2(x[i-1],y_mid,v_mid)

    return S1[N]


#Sweeps through range of initial values and finds closest points to root
for j in range(0,R.size):
    g[j] = modified_euler(R[j]) - np.log(2)
    
    if(((g[j-1] < 0) and (g[j] > 0)) or ((g[j-1] > 0) and (g[j] < 0))):
        root = secant(R[j-1],R[j],g[j-1],g[j])
        break
  
    
print('root = ', root)


#Compares solution that is found with actual solution with table and plot
modified_euler(root)
y = np.log(x)

T = np.concatenate((x,S1,y),axis=1)
T = pd.DataFrame(T)
T.columns = ['x','S1','ln(t)']
print(T.to_string(index=False))

myFigSize = (15,15)
plt.figure(figsize=myFigSize)    
plt.subplot(1,1,1)
plt.plot(x,y)
plt.plot(x,S1)
plt.grid(True)
plt.show()



