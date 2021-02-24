# -*- coding: utf-8 -*-
"""
Miguel Mares
Numerical Methods 
Homework 
Due 3/6/2020
"""
#This program solves the given differential equation: dy/dx = y^2 + 1
#using the Euler, modified Euler, and improved Euler methods for a specified
#step size, h, on the interval 0<x<1. The initial condition of y(0)=0 was given.
#The program was ran three times for step sizes h = .2, .1, .05. 

import numpy as np
import math
import matplotlib.pyplot as plt

#Specified step size, h and number of intervals N
h = 0.1
N = math.floor(1/h)
print('h =',h)

#Given initial conditions
x_0 = 0 
y_0 = 0

#Defining the f(x) function
def f(x,y):
    f = y**2 + 1
    return f

#Creating the table of calculated values
T = np.zeros((N+1,2))
T[0][0] = x_0
T[0][1] = y_0

#Euler Method
for i in range(1,N+1):
    T[i][0] = x_0 + i*h
    x_o = T[i-1][0]
    y_o = T[i-1][1]
    T[i][1] = T[i-1][1] + h*f(x_o,y_o)
print('Euler Method \n',T,'\n')

#Midpoint/Modified Euler Method
for i in range(1,N+1):
    x_mid = T[i-1][0] + .5*h
    y_mid = T[i-1][1] + (h/2)*f(T[i-1][0],T[i-1][1])
    T[i][1] = T[i-1][1] + h*f(x_mid,y_mid)
print('Modified Euler Method \n',T,'\n')
    
#Improved Euler Method
for i in range(1,N+1):
    x_i = T[i][0]
    x_o = T[i-1][0]
    y_o = T[i-1][1]
    f_o = f(x_o,y_o)
    T[i][1] = y_o + .5*h*(f_o + f(x_i, y_o + h*f_o))
print('Improved Euler Method \n',T)
