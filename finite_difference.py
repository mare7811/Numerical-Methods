# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 15:05:49 2020

@author: admin
"""
import math
import numpy as np
import matplotlib.pyplot as plt


h = .05
mu0 = .000954
T = 1000
L = 1
delta = .0005

nx = math.floor((1-0)/h) - 1
x = np.zeros(nx)
for i in range(1,nx+1):
    x[i-1] = i*h


mu = mu0+(x-L/2)*delta
a = -T/(h*h*mu)
b = -2*a
c = a

def deter(lamb):
    Det = np.zeros(nx)
    Det[0] = b[0] - lamb
    Det[1] = (b[1] - lamb)*Det[0] - a[1]*c[0]
    for i in range(2,nx):
        Det[i] = (b[i]-lamb)*Det[i-1] - a[i]*c[i-1]*Det[i-2]
    return Det[nx-1]

def secant(x1,x2,y1,y2):
    global A_l
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
        if((y[i-1]-y[i-2]) != 0):
            l_new = (r[i-2]*y[i-1] - r[i-1]*y[i-2])/(y[i-1]-y[i-2])
            r.append(l_new)
            y.append(deter(l_new))
        else:
            r.append(0)
        i += 1
    return(r[i-2])


lambda1 = 100
lambda2 = 200
determinant1 = deter(lambda1)
determinant2 = deter(lambda2)


root = secant(lambda1,lambda2,determinant1,determinant2)
omega = np.sqrt(root)
print('calculated w = ',omega)
############################################################################
N = math.floor((1-0)/h)

#Defining arrays to hold values of RK4 results and steps
x = np.zeros((N+1,1))
S1 = np.zeros((N+1,1))
S2 = np.zeros((N+1,1))

#Initial conditions for RK4 method
x_0 = 0
S1_0 = 0
#S2_0 = 3

#Defining functions for RK4 method
def F1(x,S1,S2):
    return S2

def F2(x,S1,S2):
    f2 = -((w**2)*u(x)/T)*S1
    return f2



def RK4(S2_0):
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


###########################################################################
def u(x):
    u = mu0 + (x-(L/2))*delta
    return u

x = np.zeros(nx+2)
for i in range(0,x.size):
    x[i] = i*h

w = omega
RK4(3.41856)

def w(x):
    w = (np.pi/L)*np.sqrt(T/u(x))
    return w 

z = w(x)
print(z)

print(np.average(w(x)))
S = np.sin(w(x)*np.sqrt(u(x)/T)*x)
############################################################################
myFigSize = (15,15)
plt.figure(figsize=myFigSize)

plt.subplot(2,1,1)
plt.plot(x,S1)
plt.plot(x,S)
plt.grid(True)
plt.ylabel('')
plt.xlabel('')
plt.title('')




