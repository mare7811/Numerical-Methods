# -*- coding: utf-8 -*-
"""

"""

import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt

#Define constants
g = 9.8
k = 10e-4
m = 10e-2

#Defining functions 
def F1(t,y,v):
    f1 = v
    return f1

def F2(t,y,v):
    f2 = (g - (k/m)*(v**2))
    return f2



#Interval limits
a = 0
b = 10

#Specifying step size and initial conditions
h = 0.1
t_0 = 0
y_0 = 0
v_0 = 0

#Creating arrays to hold step sizes and results
N = math.floor((b-a)/h)
t = np.zeros((N+1,1))
y = np.zeros((N+1,1))
v = np.zeros((N+1,1))

#Filling arrays with initial values
for i in range(0,N+1):
    t[i] = t_0 + i*h
y[0] = y_0
v[0] = v_0


#Calculating/filling in rest of table values
for i in range(1,N+1):
    f01 = F1(t[i-1],y[i-1],v[i-1])
    f02 = F2(t[i-1],y[i-1],v[i-1])
    
    f11 = F1((t[i-1]+(h/2)),(y[i-1]+(h/2)*f01),(v[i-1]+(h/2)*f02))
    f12 = F2((t[i-1]+(h/2)),(y[i-1]+(h/2)*f01),(v[i-1]+(h/2)*f02))

    f21 = F1((t[i-1]+(h/2)),(y[i-1]+(h/2)*f11),(v[i-1]+(h/2)*f12))
    f22 = F2((t[i-1]+(h/2)),(y[i-1]+(h/2)*f11),(v[i-1]+(h/2)*f12))
    
    f31 = F1((t[i-1]+h),(y[i-1]+h*f21),(v[i-1]+h*f22))
    f32 = F2((t[i-1]+h),(y[i-1]+h*f21),(v[i-1]+h*f22))
    
    y[i] = (y[i-1] + (h/6)*(f01 + (2*f11) + (2*f21) + f31)) 
    v[i] = (v[i-1] + (h/6)*(f02 + (2*f12) + (2*f22) + f32))

    

#Outputting table
T = np.concatenate((t,y,v),axis=1)
T = pd.DataFrame(T)
T.columns = ['t','y(t)','v(t)']
print('h = ',h)
#print(T.to_string(index=False))

#Plotting results
myFigSize = (12,12)
plt.figure(figsize=myFigSize)

plt.subplot(3,1,1)
plt.scatter(t,-y)
plt.grid(True)
plt.ylabel('y(t) [m]')
plt.xlabel('t [sec]')
plt.title('Height vs Time')

plt.subplot(3,1,2)
plt.scatter(t,-v)
plt.grid(True)
plt.ylabel('v(t) [m/s]')
plt.xlabel('t [sec]')
plt.title('Velocity vs Time')

plt.subplot(3,1,3)
plt.plot(y,v)
plt.grid(True)
plt.ylabel('v(t) [m/s]')
plt.xlabel('x(t) [m]')
plt.title('Phase Space')




N = math.floor((b-a)/h)
t_mid = np.zeros((N+1,1))
y1 = np.zeros((N+1,1))
v1 = np.zeros((N+1,1))


for i in range(0,N+1):
    t_mid[i] = t_0 + .5*h
y[0] = y_0
v[0] = v_0


for i in range(1,N+1):
    y_mid = y1[i-1] + (h/2)*F1(t[i-1],y1[i-1],v1[i-1])
    v_mid = v1[i-1] + (h/2)*F2(t[i-1],y1[i-1],v1[i-1])
    
    y1[i] = y1[i-1] + h*F1(t_mid,y_mid,v_mid)
    v1[i] = v1[i-1] + h*F2(t_mid,y_mid,v_mid)


myFigSize = (12,12)
plt.figure(figsize=myFigSize)

plt.subplot(3,1,1)
plt.scatter(t,-y1)
plt.grid(True)
plt.ylabel('y(t) [m]')
plt.xlabel('t [sec]')
plt.title('Height vs Time')


plt.subplot(3,1,2)
plt.scatter(t,-v1)
plt.grid(True)
plt.ylabel('y(t) [m]')
plt.xlabel('t [sec]')
plt.title('Height vs Time')



