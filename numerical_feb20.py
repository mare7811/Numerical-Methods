# -*- coding: utf-8 -*-
"""
Miguel Mares
Engr 428
Homework 
Due 2/12/2020
Exercise 4.1 and 4.2
"""
#This program uses linear, quadratic, and quartic integral approximation 
#methods for values of N (the number of intervals/sections) from 4 to 1028

import math
import numpy as np

#integral limits
a = 0
b = np.pi

#Defining the function to be integrated
def f(x):
    f = np.sin(x)
    return f

#Defining the first derivative
def f1(x):
    f1 = np.cos(x)
    return f1

#Defining the third derivative
def f3(x):
    f3 = -np.sin(x)
    return f3


#Creating the values of N (4,8,16,...,1024)
N = np.zeros(9)
i = 2
for j in range(0,9):
    N[j] = 2**i
    i += 1
    
#Calculating the integral approximation for each value of N    
print('Using the linear approximation: ')
for j in range(0,9):
    c = math.floor(N[j])        #Creates integer N values
    h = np.linspace(a,b,c+1)
    s = 0
    for i in range(1,c):
        s += f(h[i])*h[1]
    z = (h[1]/2)*f(h[0]) + s + (h[1]/2)*f(h[c])
    print('N = ', c,': ',z)
    

print('\nUsing the quadratic approximation: ')
for j in range(0,9):
    c = math.floor(N[j])
    h = np.linspace(a,b,c+1)
    s1 = 0
    s2 = 0
    for i in range(1,c,2):
        s1 += (h[1]/3)*4*f(h[i])
    for i in range(2,c-1,2):
        s2 += (h[1]/3)*2*f(h[i])
    z = (h[1]/3)*f(h[0]) + s1 + s2 + (h[1]/3)*f(h[c])
    print('N = ', c, ': ',z)
    
    
print('\nUsing the quartic approximation with one correction term: ')
for j in range(0,9):
    c = math.floor(N[j])        #Creates integer N values
    h = np.linspace(a,b,c+1)
    s = 0
    for i in range(1,c):
        s += f(h[i])*h[1]
    z = ((h[1]/2)*f(h[0]) + s + (h[1]/2)*f(h[c])) + (h[1]**2/12)*(f1(h[0])-f1(h[c]))
    print('N = ', c,': ',z)

print('\nUsing the quartic approximation with two correction terms: ')
for j in range(0,9):
    c = math.floor(N[j])        #Creates integer N values
    h = np.linspace(a,b,c+1)
    s = 0
    for i in range(1,c):
        s += f(h[i])*h[1]
    z = ((h[1]/2)*f(h[0]) + s + (h[1]/2)*f(h[c])) + (h[1]**2/12)*(f1(h[0])-f1(h[c])) - (h[1]**4/720)*(f3(h[0]-f3(h[c])))
    print('N = ', c,': ',z)


print('It is obvious that as the number of intervals increases (as N increases), \n') 
print('the more accurate the integral approximation is. It is also clear that \n')
print('the approximation is increasingly more accurate from the linear to quartic method \n')
print('However, it is not clear that the quartic approximation with one versus \n')
print('two corrections terms is any more accurate. I am guessing that the computer \n')
print('is only calculating to the nth decimal place and thus makes no difference \n')
print('as the corrections terms get small enough.')