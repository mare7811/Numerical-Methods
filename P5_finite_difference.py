# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 15:05:49 2020

@author: admin
"""
import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#############################################################################
#Defining Variables
#
##This section defines the constants as global variables to be used by the
#functions in the program.
#
#input:
#       none
#
#ouput:
#       none
#
#global:
#       h - step size for RK4 method and defining x-axis plot points
#       mu0 - average mass density [kg/m]
#       T - tension of piano wire [N]
#       L - length of piano wire [m]
#       delta - variation of mass density per unit length [kg/m^2]
#       nx - number of interior points on string
#       x - array to hold string position values
#       mu - array of values of mass density of string at each interior point
#       a - 'a' element of 'A-lambdaI' matrix
#       b - 'b' element of 'A-lambdaI' matrix
#       c - 'c' element of 'A-lambdaI' matrix
#       w - analytical root from textbook
#
############################################################################
h = .05 #step size
mu0 = .000954   #average mass in kilograms
T = 1000    #Tension in Newtons
L = 1   #Length of string
delta = .0005   #

#Create points along string based on step size
nx = math.floor(L/h) - 1
x = np.zeros(nx)
for i in range(1,nx+1):
    x[i-1] = i*h

#'A' matrix elements
mu = mu0+(x-L/2)*delta
a = -T/(h*h*mu)
b = -2*a
c = a

w = (np.pi/L)*np.sqrt(T/mu0)


#############################################################################
#Defining Functions
#
##This section defines the functions that are used to solve the piano wire 
#problem.
#
#input:
#       x - 
#       lamb - 
#       x1 - first guess point for secant method
#       x2 - second guess point for secant method
#       y1 - integral point from RK4 method using first guess point
#       y2 - integral point from RK4 method using second guess point
#
#ouput:
#       u - calculated mass density at given x point
#       Det[nx-1] - 
#       r[i-2] - calculated root from secant method
#
#local:
#       Det - array that holds successive determinant values
#       r - array that holds roots
#       y - array that holds secant points
#       i - iteration counter
#       l_new - new calculated lambda from secant algorithm
#
#functions:
#       u() - Function that describes the mass density of the piano string at
#              at each x position
#
#       deter() - Function that calculates the determinant of the 'A-lambdaI' 
#                 matrix for the given lambda (lamb).
#
#       secant() - Function that performs the secant algorithm and finds the 
#                  root with the given initial points.
#       
############################################################################
def u(x):
    u = mu0 + (x-(L/2))*delta
    return u


def deter(lamb):
    Det = np.zeros(nx)
    Det[0] = b[0] - lamb
    Det[1] = (b[1] - lamb)*Det[0] - a[1]*c[0]
    for i in range(2,nx):
        Det[i] = (b[i]-lamb)*Det[i-1] - a[i]*c[i-1]*Det[i-2]
    return Det[nx-1]

#Function that finds the zero using the secant method and two given points
def secant(x1,x2,y1,y2):
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


#############################################################################
#FINITE DIFFERENCE METHOD
#
##This section uses the finite difference method to solve for the fundamental 
#frequency of the "piano wire" presented in the textbook. 
#
#input:
#       none
#
#ouput:
#       none
#
#local:
#       
#       lambda1 - first guess for lambda
#       lambda2 - second guess for lambda
#       determinant1 - determinant of 'A-lambdaI' matrix from lambda1
#       determinant2 - determinatn of 'A-lambdaI' matrix from lambda2
#       root - calculated root from secant method (lambda where determinant=0)
#       omega - squareroot of root (lambda) : w = lambda^2
#       x - array that holds points along string
#       w - analytical root from textbook
#       y_linear - analytical function of linear string from textbook
#       y_f - function of nonlinear string using finite difference method
#       p_dif - percent difference between linear and nonlinear fundamental 
#               frequency
#
##############################################################################
lambda1 = 100   #first guess for lambda
lambda2 = 200   #second guess for lambda
determinant1 = deter(lambda1)
determinant2 = deter(lambda2)

#root found from secant method
root = secant(lambda1,lambda2,determinant1,determinant2)
omega = np.sqrt(root)


###################### String Equations #####################################
#The original x matrix does not include endpoints so it is redefined to include 
#them
x = np.zeros(nx+2)
for i in range(0,x.size):
    x[i] = i*h

y_linear = np.sin(w*np.sqrt(mu0/T)*x)
y_f = np.sin(omega*np.sqrt(mu0/T)*x)

p_dif = (np.abs(w-omega)/w)*100

print("Finite Difference Method")
print('[nonlinear]  w = ',omega)
print('[linear] w = ', w)
print('Percent change:',p_dif,'%')

############################## Plotting ####################################
myFigSize = (15,15)
plt.figure(figsize=myFigSize)
plt.subplot(2,1,1)
plt.plot(x,y_f)
plt.plot(x,y_linear)
plt.grid(True)
plt.ylabel('y(x) [m]')
plt.xlabel('x [m]')
plt.title('Plot of strings')
plt.show()

########################## Outputting Table ####################################################
x.resize((x.size,1))
y_f.resize((y_f.size,1))
y_linear.resize((y_linear.size,1))
T = np.concatenate((x,y_f,y_linear),axis=1)
T = pd.DataFrame(T)
T.columns = ['x','y(x) [nonlinear]','y(x) [linear]']
print(T)

