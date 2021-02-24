# -*- coding: utf-8 -*-
"""
Miguel Mares
Numerical Methods
Homework 
Due 2/20/2020
"""
#This program performs Romberg integration for the complete eliptic integral of
#the first kind, then calculates the period of a simple pendulum for initial 
#angles from 0 to 360 degrees in steps of 5 degrees. A value of l is arbitrarily
#set.

import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#for output display settings
pd.set_option('display.max_colwidth',100)
pd.set_option('display.max_columns',100)
pd.set_option('display.max_rows',100)

#Set limits of integral
a = 0                            
b = np.pi/2

#Set arbitrary length for pendulum and value of g
l = 10
g = 9.81

#Create array of angles and corresponding k values for K(k) equation.
steps = 5
angle = np.arange(0,360+steps,steps)
k1 = np.sin((angle*np.pi)/180)


#Set number of rows (values of m) to calculate up to for Romberg Table
r = 10     #rows
c = r+3      #columns


#Create matrix with r rows and c columns (Romberg Table)
A = np.zeros((r,c))
for m in range(0,r):
    A[m][0] = m                 #m's
    A[m][1] = 2**m              #N's
    A[m][2] = (b-a)/A[m][1]     #h's

#Create table for angles, K's, period, etc
Table = np.zeros((angle.size,5))


#This loops through each angle 0 to 360 and creates a Romberg table. The table 
#is calculated up to m=10 to maximize accuracy and thus uses the last entry in
#the table as the integral calculation. The K(k) value and corresponding period 
#value are then placed into the table.
for w in range(0,angle.size):
    def f(z):
        f = 1/np.sqrt(1-(k1[w]**2)*np.sin(z)**2)
        return f
    
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
    
    Table[w][2] = A[r-1][r+2]       #last entry in table - k's
    Table[w][3] = 4*np.sqrt(l/g)*Table[w][2]*np.sin((angle[w]*np.pi/180)/2) #T's
    Table[w][4] = l*np.sin((angle[w]*2)*np.pi/180)  #Amplitudes

#Create output table and format
Ts = np.zeros(angle.size)
for i in range(0,angle.size):
    Table[i][0] = angle[i]
    Table[i][1] = angle[i]*2
    Ts[i] = Table[i][3]
T = pd.DataFrame(Table)
T.columns = ['sin^-1(k)','angle','K(k)','T - period','Amplitude']
print(T)

#Print Plot of period vs angle
myFigSize = (12,12)
plt.figure(figsize=myFigSize)
plt.subplot(1,1,1)
plt.plot(angle*2,Ts)
plt.grid(True)
plt.ylabel('Period [s]')
plt.xlabel('Initial Angle [deg]')
plt.title('Period of Pendulum vs Angle')
print('2*pi*sqrt(l/g) = ',2*np.pi*np.sqrt(l/g))
print('At 180 deg and 360 deg there are discrepencies.')
print('At these angles, the period is infinite, and therefore incalculable.')