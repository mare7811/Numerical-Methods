# -*- coding: utf-8 -*-
"""
Miguel Mares
Numerical Methods
Exercise 4.3
Due 2/14/2020
"""
#This program performs Romberg integration for the given function:y=x^-2*sin(x)
#but changes the variable x to x=1/t. The function then becomes y=t^2sin(1/t)
#and the limits become a=1 and b=0. Reversing the limits gives the equation 
#integrated from 0 to 1 for y=-t^2*sin(1/t).

import math
import numpy as np
import pandas as pd

#for output display settings
pd.set_option('display.max_colwidth',100)
pd.set_option('display.max_columns',100)

#Set limits of integral
a = 1e-9            #x=0, t=1      
b = 1               #x=infinity,t=0

#Define function to be integrated
def f(t):
    f = -((t**2)*np.sin(1/t))
    return f

#Set number of rows (values of m) to calculate up to
r = 5       #rows
c = r+3      #columns

#Create matrix with r rows and c columns
A = np.zeros((r,c))
for m in range(0,r):
    A[m][0] = m                 #m's
    A[m][1] = 2**m              #N's
    A[m][2] = (b-a)/A[m][1]     #h's


#This loop calculates the values for k=1 using the trapezoid rule
for m in range(0,r):
    c = math.floor(A[m][1])         #Creates integer N values
    h = np.linspace(a,b,c+1)        #Creates points for intervals
    s = 0
    for i in range(1,c):
        s += f(h[i])*h[1]
    A[m][3] = (h[1]/2)*f(h[0]) + s + (h[1]/2)*f(h[c])   #Calculates final approx.

#This loop calculates the rest of the values in the table
k = 1
for k in range(1,c-k-1):    
    for m in range(0,r-k):
        A[m+k][k+3] = ((4**k)*A[m+k][k+2] - A[m+k-1][k+2])/((4**k)-1)
 
#Prints out filled out matrix
print('       m       N       h_m      Tm,0      Tm,1      Tm,2      Tm,3 Tm,4 Tm,5')
b = pd.DataFrame(A)
print(b)

print('The integral is calculated to be approximately',A[r-1][r+2],)
