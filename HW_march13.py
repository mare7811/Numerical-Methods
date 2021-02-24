# -*- coding: utf-8 -*-
"""
Miguel Mares
Numerical Methods
Homework 
Due 3/13/2020
"""
#This program applies the RK4 method with h = 0.1 to to the "consider" exercise
#from (3/11) on the interval 0 <= t <= 1.


import numpy as np
import math
import pandas as pd


#Defining functions 
def F1(x,S1,S2):
    f1 = S2
    return f1

def F2(x,S1,S2):
    f2 = -4*S1
    return f2



#Interval limits
a = 0
b = 2*np.pi

#Specifying step size and initial conditions
h = 0.1
x_0 = 0
S_01 = 1
S_02 = 0

#Creating arrays to hold step sizes and results
N = math.floor((b-a)/h)
X = np.zeros((N+1,1))
S1 = np.zeros((N+1,1))
S2 = np.zeros((N+1,1))

#Filling arrays with initial values
for i in range(0,N+1):
    X[i] = x_0 + i*h
S1[0] = S_01
S2[0] = S_02


#Calculating/filling in rest of table values
for i in range(1,N+1):
    f01 = F1(X[i-1],S1[i-1],S2[i-1])
    f02 = F2(X[i-1],S1[i-1],S2[i-1])
    
    f11 = F1((X[i-1]+(h/2)),(S1[i-1]+(h/2)*f01),(S2[i-1]+(h/2)*f02))
    f12 = F2((X[i-1]+(h/2)),(S1[i-1]+(h/2)*f01),(S2[i-1]+(h/2)*f02))

    f21 = F1((X[i-1]+(h/2)),(S1[i-1]+(h/2)*f11),(S2[i-1]+(h/2)*f12))
    f22 = F2((X[i-1]+(h/2)),(S1[i-1]+(h/2)*f11),(S2[i-1]+(h/2)*f12))
    
    f31 = F1((X[i-1]+h),(S1[i-1]+h*f21),(S2[i-1]+h*f22))
    f32 = F2((X[i-1]+h),(S1[i-1]+h*f21),(S2[i-1]+h*f22))
    
    S1[i] = S1[i-1] + (h/6)*(f01 + (2*f11) + (2*f21) + f31) 
    S2[i] = S2[i-1] + (h/6)*(f02 + (2*f12) + (2*f22) + f32) 

    


#Outputting table
T = np.concatenate((X,S1,S2),axis=1)
T = pd.DataFrame(T)
T.columns = ['t','y1','y2']
print(T.to_string(index=False))
