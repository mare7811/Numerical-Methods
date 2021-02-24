# -*- coding: utf-8 -*-
"""
Miguel Mares
Numerical Methods
Homework 
Due 3/11/2020
"""
#This program solves exercise 5.8. The Runge-Kutta_Fehlberg method is used
#to solve the differential equation: dy/dx = 1/(x^2) with y(-1) = 1 on the 
#interval -1 <= x <= 1. The size of h was to be between .001 and .1
#with an error of at most 5e-5 (epsilon). Additionally, h was to be kept from
#changing to much: no less than .1h and no more than 4h. 

import numpy as np

#f(x,y) function from dy/dx = f(x,y)
def f(x,y):
    f = 1/(x**2)
    return f

#interval bounds
a = -1
b = 1

#error and max and min h size
epsilon = 5e-5
h_max = .1
h_min = .001

#Initial conditions
x_0 = -1
y_0 = 1


n = 4   #n'th order method that is used 
h = h_max   #Starting step size is h max
x0 = np.zeros(50)   #Create array large enough to calculate across entire interval
y0 = np.zeros(50)
x0[0] = x_0     #Set initial values in table
y0[0] = y_0


i = 1   #iteration counter
while(x0[i-1] != b):
    f0 = f(x0[i-1], y0[i-1])
    f1 = f(x0[i-1]+(h/4), y0[i-1]+(h/4)*f0)
    f2 = f(x0[i-1]+(3*h/8), y0[i-1]+(3*h/32)*f0+(9*h/32)*f1)
    f3 = f(x0[i-1]+(12*h/13), y0[i-1]+(1932*h/2197)*f0-(7200*h/2197)*f1+(7296*h/2197)*f2)
    f4 = f(x0[i-1]+h, y0[i-1]+(439*h/216)*f0-(8*h)*f1+(3680*h/513)*f2-(845*h/4104)*f3)
    f5 = f(x0[i-1]+(h/2), y0[i-1]-(8*h/27)*f0+(2*h)*f1-(3544*h/2565)*f2+(1859*h/4104)*f3-(11*h/40)*f4)

    y = y0[i-1] + h*((25/216)*f0 + (1408/2565)*f2 + (2197/4104)*f3 - (1/5)*f4)
    y_hat = y0[i-1] + h*((16/135)*f0 + (6656/12825)*f2 + (28561/56430)*f3 - (9/50)*f4 + (2/55)*f5)

    h_new = .9*h*(epsilon/np.abs((y-y_hat)/y_hat))**(1/(n+1))
   
    #Check if h_new is between h_min and h_max
    if((h_new) < .001):
        h_new = .001
    if(h_new > .1):
        h_new = .1
        
    #Check if h_new is smaller or larger than h    
    if(h_new < h):
        if(h_new < .1*h):   #Check if h changes more than .1h
            h_new = .1*h
        h = h_new           #Correct h without accepting previous h value
    else:
        if(h_new > 4*h):    #Check if h changes more than 4h
            h_new = 4*h
        x0[i] = x0[i-1] + h     #Accept h value and enter in table
        y0[i] = y               #Accept y value
        if(x0[i] > b):          #Check if x0+h value is beyond interval 
            h = b - x0[i-1]     #If beyond interval, adjust h to land exactly on interval boundary value
        else:                   #If h is not beyond interval, set h to h_new and continue loop
            h = h_new
            i += 1
        
        
#Print entries in a table
x0 = np.reshape(x0,(50,1))
y0 = np.reshape(y0,(50,1))
T = np.concatenate((x0,y0),axis=1)
print(T)
      
