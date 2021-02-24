# -*- coding: utf-8 -*-
"""
Created on Wed Mar 25 17:02:17 2020

@author: admin
"""
import numpy as np
import matplotlib.pyplot as plt
import math
import pandas as pd

def F1(x,S1,S2):
    f1 = S2
    return f1
    
def F2(x,S1,S2):
    f2 = -(S2)**2 - S1 + np.log(x)
    return f2


h = 0.25
x_0 = 1
S1_0 = 0

a = 1
b = 2
N = math.floor((b-a)/h)

#Defining arrays to hold values
x = np.zeros((N+1,1))
S1 = np.zeros((N+1,1))
S2 = np.zeros((N+1,1))

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



t = np.arange(0,2.1,.1)
y = np.zeros((t.size,1))
for i in range(0,t.size):
    y[i] = modified_euler(t[i]) - np.log(2)


t = np.array(t)
t = np.resize(t,(t.size,1))
T = np.concatenate((t,y),axis=1)
T = pd.DataFrame(T)
T.columns = ['x','y']
print(T.to_string(index=False))
    

r = np.zeros((10,1))
z = np.zeros((10,1))
r[0] = 1
r[1] = 1.1
z[0] = -.029054
z[1] = .012179
for i in range(2,10):
    if((z[i-1]-z[i-2]) != 0 ):
        r[i] = (r[i-2]*z[i-1] - r[i-1]*z[i-2])/(z[i-1]-z[i-2])
    else:
        break

T = np.concatenate((r,z),axis=1)
T = pd.DataFrame(T)
T.columns = ['roots','y']
print(T.to_string(index=False))

y = np.log(x)
modified_euler(1.07046298)

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
