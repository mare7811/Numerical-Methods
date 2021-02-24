# -*- coding: utf-8 -*-
"""
Miguel Mares
Numerical Methods
Homework 
Due 4/1/2020
"""
#This program solves the projectile motion problem presented in the textbook 
#using the RK4 method the problem given is a second order differential equation
#in two dimensions. - Exercise 5.24

import numpy as np
import matplotlib.pyplot as plt
import math
import pandas as pd

#Defining functions
def F1(t,S1,S2,S3,S4):  #x(t)
    f1 = S2
    return f1
    
def F2(t,S1,S2,S3,S4):  #dx/dt
    f2 = -c*np.sqrt((S2**2)+(S4**2))*S2
    return f2

def F3(t,S1,S2,S3,S4):  #y(t)
    f3 = S4
    return f3

def F4(t,S1,S2,S3,S4):  #dy/dt
    f4 = -g-c*np.sqrt((S2**2)+(S4**2))*S4
    return f4

#Interval limits
a = 0
b = 10

#constants
g = 9.8
c = .1

#Specifying step size and initial conditions - given in problem statement
h = 0.1 
t_0 = 0
S_01 = 0
S_02 = 20
S_03 = 0
S_04 = 0

#Creating arrays to hold step sizes and results
N = math.floor((b-a)/h)
t = np.zeros((N+1,1))
S1 = np.zeros((N+1,1))
S2 = np.zeros((N+1,1))
S3 = np.zeros((N+1,1))
S4 = np.zeros((N+1,1))


#Filling arrays with initial values
for i in range(0,N+1):
    t[i] = t_0 + i*h
S1[0] = S_01
S2[0] = S_02
S3[0] = S_03
S4[0] = S_04

#Calculating/filling in rest of table values
for i in range(1,N+1):
    f01 = F1(t[i-1],S1[i-1],S2[i-1],S3[i-1],S4[i-1])
    f02 = F2(t[i-1],S1[i-1],S2[i-1],S3[i-1],S4[i-1])
    f03 = F3(t[i-1],S1[i-1],S2[i-1],S3[i-1],S4[i-1])
    f04 = F4(t[i-1],S1[i-1],S2[i-1],S3[i-1],S4[i-1])
    
    f11 = F1((t[i-1]+(h/2)),(S1[i-1]+(h/2)*f01),(S2[i-1]+(h/2)*f02),(S3[i-1]+(h/2)*f03),(S4[i-1]+(h/2)*f04))
    f12 = F2((t[i-1]+(h/2)),(S1[i-1]+(h/2)*f01),(S2[i-1]+(h/2)*f02),(S3[i-1]+(h/2)*f03),(S4[i-1]+(h/2)*f04))
    f13 = F3((t[i-1]+(h/2)),(S1[i-1]+(h/2)*f01),(S2[i-1]+(h/2)*f02),(S3[i-1]+(h/2)*f03),(S4[i-1]+(h/2)*f04))
    f14 = F4((t[i-1]+(h/2)),(S1[i-1]+(h/2)*f01),(S2[i-1]+(h/2)*f02),(S3[i-1]+(h/2)*f03),(S4[i-1]+(h/2)*f04))

    f21 = F1((t[i-1]+(h/2)),(S1[i-1]+(h/2)*f11),(S2[i-1]+(h/2)*f12),(S3[i-1]+(h/2)*f13),(S4[i-1]+(h/2)*f14))
    f22 = F2((t[i-1]+(h/2)),(S1[i-1]+(h/2)*f11),(S2[i-1]+(h/2)*f12),(S3[i-1]+(h/2)*f13),(S4[i-1]+(h/2)*f14))
    f23 = F3((t[i-1]+(h/2)),(S1[i-1]+(h/2)*f11),(S2[i-1]+(h/2)*f12),(S3[i-1]+(h/2)*f13),(S4[i-1]+(h/2)*f14))
    f24 = F4((t[i-1]+(h/2)),(S1[i-1]+(h/2)*f11),(S2[i-1]+(h/2)*f12),(S3[i-1]+(h/2)*f13),(S4[i-1]+(h/2)*f14))

    f31 = F1((t[i-1]+h),(S1[i-1]+h*f21),(S2[i-1]+h*f22),(S3[i-1]+h*f23),(S4[i-1]+h*f24))
    f32 = F2((t[i-1]+h),(S1[i-1]+h*f21),(S2[i-1]+h*f22),(S3[i-1]+h*f23),(S4[i-1]+h*f24))
    f33 = F3((t[i-1]+h),(S1[i-1]+h*f21),(S2[i-1]+h*f22),(S3[i-1]+h*f23),(S4[i-1]+h*f24))
    f34 = F4((t[i-1]+h),(S1[i-1]+h*f21),(S2[i-1]+h*f22),(S3[i-1]+h*f23),(S4[i-1]+h*f24))

    #RK4 alorithm
    S1[i] = S1[i-1] + (h/6)*(f01 + (2*f11) + (2*f21) + f31) 
    S2[i] = S2[i-1] + (h/6)*(f02 + (2*f12) + (2*f22) + f32) 
    S3[i] = S3[i-1] + (h/6)*(f03 + (2*f13) + (2*f23) + f33) 
    S4[i] = S4[i-1] + (h/6)*(f04 + (2*f14) + (2*f24) + f34) 

    
print('h = ',h)

#Outputting table - formatting 
T = np.concatenate((t,S1,S2,S3,S4),axis=1)
T = pd.DataFrame(T)
T.columns = ['t','x(t)','dx/dt','y','dy/dt']
print(T.to_string(index=False))

#Plotting results
myFigSize = (12,12)
plt.figure(figsize=myFigSize)
plt.subplot(1,1,1)
plt.plot(t,S1)
plt.plot(t,S2)
plt.grid(True)
plt.legend(['x(t)','dx/dt'])
plt.ylabel('x, v_x(t) [m]/[m/s]')
plt.xlabel('t [sec]')
plt.title('Plot of results')

myFigSize = (12,12)
plt.figure(figsize=myFigSize)
plt.subplot(1,1,1)
plt.plot(t,S3)
plt.plot(t,S4)
plt.grid(True)
plt.legend(['y(t)','dy/dt'])
plt.ylabel('y, v_y(t) [m]/[m/s]')
plt.xlabel('t [sec]')
plt.title('Plot of results')
plt.show()

print('The results are expected. The velocity in the x direction is shown to',
      'decrease from 20 m/s to 0. The plots show this happens and the position',
      'in the x direction is shown to increase and eventually come to a stop.')
print('The velocity in the y direction is shown to be negative and level off',
      'at around -9.8 m/s. The position data is shown to go negative and keep',
      'decreasing at a constant rate.')

