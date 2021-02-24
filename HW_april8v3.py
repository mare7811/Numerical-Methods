# -*- coding: utf-8 -*-
"""
Miguel Mares
Numerical Methods
Homework
Due 4/8/2020
"""
#Exercise 5.37
#This program solves exercise 5.37 of the textbook with the give conditions:
#y'' = -((w^2)*u_0/T)*y; y(0) = 0, y(1) = 0, u_0 = .954 grams, T = 1000 Newtons
#for the fundamental frequency ,w1, of a vibrating string.

import numpy as np
import math
import matplotlib.pyplot as plt

#RK4 method - returns final value of solved second order ODE 
def RK4():
    global S1
    S2 = np.zeros((N+1,1))
    S1 = np.zeros((N+1,1))

    #Filling arrays with initial values
    for i in range(0,N+1):
        x[i] = x_0 + i*h
    S1[0] = S1_0
    S2[0] = S2_0

    #Calculating/filling in rest of table values
    for i in range(1,N+1):
        f01 = F1(x[i-1],S1[i-1],S2[i-1])
        f02 = F2(x[i-1],S1[i-1],S2[i-1])        
        f11 = F1((x[i-1]+(h/2)),(S1[i-1]+(h/2)*f01),(S2[i-1]+(h/2)*f02))
        f12 = F2((x[i-1]+(h/2)),(S1[i-1]+(h/2)*f01),(S2[i-1]+(h/2)*f02))        
        f21 = F1((x[i-1]+(h/2)),(S1[i-1]+(h/2)*f11),(S2[i-1]+(h/2)*f12))
        f22 = F2((x[i-1]+(h/2)),(S1[i-1]+(h/2)*f11),(S2[i-1]+(h/2)*f12))    
        f31 = F1((x[i-1]+h),(S1[i-1]+h*f21),(S2[i-1]+h*f22))
        f32 = F2((x[i-1]+h),(S1[i-1]+h*f21),(S2[i-1]+h*f22))
    
        S1[i] = (S1[i-1] + (h/6)*(f01 + (2*f11) + (2*f21) + f31)) 
        S2[i] = (S2[i-1] + (h/6)*(f02 + (2*f12) + (2*f22) + f32))
    print(S2[N])
    return(S1[N])

#Function that finds root using secant method
def secant(x1,x2,y1,y2):
    global w
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
            w = (r[i-2]*y[i-1] - r[i-1]*y[i-2])/(y[i-1]-y[i-2])
            r.append(w)
            c = RK4()
            y.append(c)
        else:
            r.append(0)
        i += 1
    return(r[i-2])
mu0 = .000954
L = 1
delta = .0005

def u(x):
    u = mu0 + (x-(L/2))*delta
    return u
#Defining functions for RK4 method
def F1(x,S1,S2):
    return S2

def F2(x,S1,S2):
    f2 = -((w**2)*u_0/T)*S1
    return f2

#Defining constants
L = 1
T = 1000
u_0 = .000954
w = 0

#step size, Interval limits, and number of intervals
h = 0.05
a = 0
b = L
y_lower = 0
y_upper = 0
N = math.floor((b-a)/h)

#Defining arrays to hold values of RK4 results and steps
x = np.zeros((N+1,1))
S1 = np.zeros((N+1,1))
S2 = np.zeros((N+1,1))

#Initial conditions for RK4 method
x_0 = a
S1_0 = y_lower
S2_0 = 3

#initial guesses of w, frequency
w_guess1 = 2000
w_guess2 = 3000

w = w_guess1
g1 = RK4()

w = w_guess2
g2 = RK4()

#Finding root, using initial w guesses and corresponding RK4 result
root = secant(w_guess1,w_guess2,g1,g2)
print('calculated w = ', root)

w = (np.pi/L)*np.sqrt(T/mu0)
print('Real w (for n=1):', w)

print('The calculated vs real w is very close, with error mostly coming from'
      ' the secant method function.')
print((np.abs(w-root)/w)*100)
w = (np.pi/L)*np.sqrt(T/mu0)
S = np.sin(w*np.sqrt(mu0/T)*x)

myFigSize = (15,15)
plt.figure(figsize=myFigSize)

plt.subplot(2,1,1)
plt.plot(x,S1)
plt.plot(x,S)
plt.grid(True)
plt.ylabel('y(t) [m]')
plt.xlabel('t [sec]')
plt.title('Height vs Time')