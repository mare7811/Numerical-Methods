# -*- coding: utf-8 -*-
"""

"""

import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt

###################################################################
#Global variables and functions
#
#This section defines the globally used variables and functions.
#
#input:
#       none
#
#ouput:
#       none
#
#local:
#       g - constant value for acceleration of gravity in m/s^2
#       k - drag coefficient for a particular sphere in kg/m
#       m - mass of particular sphere in kg
#       h - step sized used for each Runge Kutta method
#       t_0 - initial time value
#       y_0 - initial height value
#       v_0 - initial velocity value
#       a - lower interval boundary limit
#       b - upper interval boundary limit
#       N - number of intervals for the specified step size and interval
#
#functions:
#       F1() - Function describing derivative of x(t) 
#       
#       F2() - Function describing second derivative of x(t)
#       
###################################################################

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

#Specifying step size and initial conditions
h = 0.05
t_0 = 0
y_0 = 0
v_0 = 0

#Interval limits and number of intervals
a = 0
b = 10
N = math.floor((b-a)/h)


####################################################################
#Modified Euler Method
#
#This section applies the Modified Euler Method to solve the second 
#order differential equation x'' = g - (k/m)(x')^2.
#
#input:
#       none
#
#output:
#       T - 2-dimensional array that holds all values for outputting 
#
#local:
#       t_m - array to hold time values
#       y_m - array to hold height values
#       v_m - array to hold velocity values
#       y_mid - midpoint value of y
#       v_mid - midpoint value of v
#
####################################################################

#Defining arrays to hold values
t_m = np.zeros((N+1,1))
y_m = np.zeros((N+1,1))
v_m = np.zeros((N+1,1))

#Filling time array with step values and other arrays with initial values
for i in range(0,N+1):
    t_m[i] = t_0 + i*h
y_m[0] = y_0
v_m[0] = v_0


#Filling/calculating rest of array values using Modified Euler 
for i in range(1,N+1):
    y_mid = y_m[i-1] + (h/2)*F1(t_m[i-1],y_m[i-1],v_m[i-1])
    v_mid = v_m[i-1] + (h/2)*F2(t_m[i-1],y_m[i-1],v_m[i-1])
    
    y_m[i] = y_m[i-1] + h*F1(t_m[i-1],y_mid,v_mid)
    v_m[i] = v_m[i-1] + h*F2(t_m[i-1],y_mid,v_mid)


#Outputting table
T = np.concatenate((t_m,y_m,v_m),axis=1)
T = pd.DataFrame(T)
T.columns = ['t','y(t)','v(t)']
print('Modified Euler Method \n','h = ',h)
print(T.to_string(index=False))

#Plotting results
myFigSize = (15,15)
plt.figure(figsize=myFigSize)

plt.subplot(2,1,1)
plt.scatter(t_m,-y_m)
plt.grid(True)
plt.ylabel('y(t) [m]')
plt.xlabel('t [sec]')
plt.title('Height vs Time')

plt.subplot(2,1,2)
plt.scatter(t_m,-v_m)
plt.grid(True)
plt.ylabel('v(t) [m/s]')
plt.xlabel('t [sec]')
plt.title('Velocity vs Time')
plt.show()


#####################################################################
#RK4 Method
#
#This section applies the RK4 Method to solve the same differential
#equation as before.
#
#input:
#       none
#
#output:
#       T - 2-dimensional array that holds all values for outputting 
#
#local:
#       t - array to hold time values
#       y - array to hold height values
#       v - array to hold velocity values
#       f01 - value used in RK4 algorithm calulation
#       f02 - value used in RK4 algorithm calculation
#       f11 - value used in RK4 algorithm calculation
#       f12 - value used in RK4 algorithm calculation
#       f21 - value used in RK4 algorithm calculation
#       f22 - value used in RK4 algorithm calculation
#       f31 - value used in RK4 algorithm calculation
#       f32 - value used in RK4 algorithm calculation
#
#####################################################################
#Creating arrays to hold step sizes and results
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
print('RK4 Method \n','h = ',h)
print(T.to_string(index=False))

#Plotting results
myFigSize = (15,15)
plt.figure(figsize=myFigSize)

plt.subplot(2,1,1)
plt.scatter(t,-y)
plt.grid(True)
plt.ylabel('y(t) [m]')
plt.xlabel('t [sec]')
plt.title('Height vs Time')

plt.subplot(2,1,2)
plt.scatter(t,-v)
plt.grid(True)
plt.ylabel('v(t) [m/s]')
plt.xlabel('t [sec]')
plt.title('Velocity vs Time')
plt.show()






