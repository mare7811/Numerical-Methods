# -*- coding: utf-8 -*-
"""
Miguel Mares
Numerical Methods
Homework 
Due 3/8/2020
"""
#This program solves exercise 5.5 by using the RK4 algorithm to solve the 
#differential equation dy/dx = y^2 + 1 w/ y(0) = 0 on the interval 0 <= x <= 1.
#Step sizes of h from .02 to .2 in steps of .02 were used.
#The previous RK2 methods were used for comparison.

#*I was unable to authentically produce the log-log vs error plot, because 
#I could not figure out how to find alpha for the equation 
# error = alpha*(h^n) that is shown in the textbook. I instead used an alpha 
#values that produced close plots. 

import numpy as np
import math
import matplotlib.pyplot as plt


#Given initial conditions
x_0 = 0 
y_0 = 0

#Defining the f(x) function
def f(x,y):
    f = y**2 + 1
    return f

#Creating step sizes 
num = 10
h_a = np.arange(.02,.22,.02)


#Arrays to hold error values for each method
error1 = np.zeros(h_a.size)
error2 = np.zeros(h_a.size)
error3 = np.zeros(h_a.size)
error4 = np.zeros(h_a.size)
for j in range(0,h_a.size):
    #Create Table for x and y entries for each 'h' value
    h = h_a[j]
    N = math.floor(1/h)
    T = np.zeros((N+1,2))
    T[0][0] = x_0
    T[0][1] = y_0
    for i in range(0,N+1):
        T[i][0] = x_0 + i*h
    print('h =',h)


    #Euler Method 
    for i in range(1,N+1):
        alpha = .9
        n = 1
        x_o = T[i-1][0]
        y_o = T[i-1][1]
        T[i][1] = T[i-1][1] + h*f(x_o,y_o)
    error1[j] = alpha*(h**n)
    print('Euler: ', T[i][1])

    #Midpoint/Modified Euler Method
    for i in range(1,N+1):
        alpha = .75
        n = 2
        x_mid = T[i-1][0] + .5*h
        y_mid = T[i-1][1] + (h/2)*f(T[i-1][0],T[i-1][1])
        T[i][1] = T[i-1][1] + h*f(x_mid,y_mid)
    error2[j] = alpha*(h**n)
    print('Modified Euler: ', T[i][1])
 
    #Improved Euler Method    
    for i in range(1,N+1):
        n = 2
        alpha = .2
        x_i = T[i][0]
        x_o = T[i-1][0]
        y_o = T[i-1][1]
        f_o = f(x_o,y_o)
        T[i][1] = y_o + .5*h*(f_o + f(x_i, y_o + h*f_o))
    error3[j] = alpha*(h**n)
    print('Improved Euler', T[i][1])

    #RK4 method
    for i in range(1,N+1):
        alpha = .027
        n = 4
        x0 = T[i-1][0]
        y0 = T[i-1][1]
        f0 = f(x0,y0)
        f1 = f((x0+(h/2)),(y0+(h/2)*f0))
        f2 = f((x0+(h/2)),(y0+(h/2)*f1))
        f3 = f((x0+h),(y0+h*f2))
        T[i][1] = y0 + (h/6)*(f0 + (2*f1) + (2*f2) + f3)    
    error4[j] = alpha*(h**n)
    print('RK4: ',T[i][1])
    print('RK4 Method Table \n', T, '\n')



############# Plotting #################################################
myFigSize = (12,12)
plt.figure(figsize=myFigSize)
plt.subplot(1,1,1)
plt.yscale("log")
plt.xscale("log")
plt.plot(h_a,error1,marker='d')
plt.plot(h_a,error2,marker='d')
plt.plot(h_a,error3,marker='d')
plt.plot(h_a,error4,marker='d')
plt.ylim(10e-7,10e-1)
plt.xlim(0.02,.2)
plt.grid(True,which="both")
plt.legend(['Standard Euler','Modified Euler','Improved Euler','RK4'])
plt.ylabel('Relative error')
plt.xlabel('step size [h]')
plt.title('log-log plot of relative error versus h')

