# -*- coding: utf-8 -*-
"""
Miguel Mares
Numerical Methods
Homework 
Due 4/6/2020
"""
#This program solves exercise 5.34 using the LU decomp, Gauss-Seidel, and SOR 
#method.


import numpy as np
import math
import pandas as pd
pd.options.display.float_format='{:,.2f}'.format    #format table for two decimal places
pd.set_option('display.max_columns',500)    

#Set step size, boundary limits, number of intervals
h = .1
a = 0
b = 1
N = math.floor((b-a)/h)

#Set tolerance and boundary values and alpha 
tolerance = 2e-4
y_lower = 0
y_upper = 100
alpha = .5

#Define 'initial guess' straight line function for first iteration
def y0(x):
    y0 = ((y_upper-y_lower)/(b-a))*x
    return y0

#Define x values
x = np.zeros(N+1)
for i in range(0,N+1):
    x[i] = a + i*h
p = x
    
#Fill in first iteration array
Y0 = np.zeros(N+1)
for i in range(0,N+1):
    Y0[i] = y0(x[i])


########################################################################
#Gauss Seidel Method
########################################################################
#Create another array for second iteration and Table array
Y = np.zeros(N+1)
Y1 = np.zeros(N+1)
T1 = np.vstack((Y0,Y1))
j = 1
done = False
while(done == False):   #loop until tolerance is reached
    T1[j][0] = y_lower   #Set boundary values in array
    T1[j][N] = y_upper
    for i in range(1,N):    #Fill in iteration values
        T1[j][i] = (1/(2-10*h**2))*((1-(5*h/2))*T1[j-1][i+1] + (1+(5*h/2))*T1[j][i-1] - 10*(h**2)*x[i])
        if(np.abs((T1[j][i]-T1[j-1][i])/T1[j][i]) < tolerance):    #check tolerance
            done = True   
    if(done == False):  #Add another array to table if another iteration is needed
        T1 = np.vstack((T1,Y))
    j += 1
    
#output table
T = pd.DataFrame(T1)
T.columns = ['y(0)','y(.1)','y(.2)','y(.3)','y(.4)','y(.5)','y(.6)','y(.7)','y(.8)','y(.9)','y(1)']
print('\nGauss-Seidel Method')
print(T)

    
########################################################################
#SOR Method
#######################################################################
T2 = np.vstack((Y0,Y1))
j = 1
done = False
while(done == False):   #loop until tolerance is reached
    T2[j][0] = y_lower   #Set boundary values in array
    T2[j][N] = y_upper
    for i in range(1,N):    #Fill in iteration values 
        T2[j][i] = (1/(2-10*h**2))*((1-(5*h/2))*T2[j-1][i+1] + (1+(5*h/2))*T2[j][i-1] - 10*(h**2)*x[i])
        T2[j][i] = T2[j][i] + alpha*(T2[j][i] - T2[j-1][i])
        if(np.abs((T2[j][i]-T2[j-1][i])/T2[j][i]) < tolerance):    #check tolerance
            done = True 
    if(done == False):
        T2 = np.vstack((T2,Y))
    j += 1


T = pd.DataFrame(T2)
T.columns = ['y(0)','y(.1)','y(.2)','y(.3)','y(.4)','y(.5)','y(.6)','y(.7)','y(.8)','y(.9)','y(1)']
print('\nSOR Method')
print('\n',T)

##################### Finite Difference Method################################
N = 9
y_a = (2+5*h)/(2*h**2)
y_b = (-2+(10*(h**2)))/h**2                                                                                                    
y_c = (2-5*h)/(2*h**2)
t = np.arange(.1,1,.1)

A = [[y_b,y_c,0,0,0,0,0,0,0],
     [y_a,y_b,y_c,0,0,0,0,0,0],
     [0,y_a,y_b,y_c,0,0,0,0,0],
     [0,0,y_a,y_b,y_c,0,0,0,0],
     [0,0,0,y_a,y_b,y_c,0,0,0],
     [0,0,0,0,y_a,y_b,y_c,0,0],
     [0,0,0,0,0,y_a,y_b,y_c,0],
     [0,0,0,0,0,0,y_a,y_b,y_c],
     [0,0,0,0,0,0,0,y_a,y_b]]

b = np.zeros((N,1))
for i in range(0,N):
    b[i] = 10*t[i]
b[0] = b[0] - y_a*y_lower
b[N-1] = b[N-1] - y_c*y_upper

############# LU DECOMP #############################################
#Creates the L matrix and sets known entries.
L = np.zeros((N,N))
for i in range(0,N):
    L[i][0] = A[i][0]
    

#Creates the U matrix and sets known entries.
U = np.zeros((N,N))
for i in range(0,N):
    U[0][i] = A[0][i]/L[0][0]

#The following functions calculate the sums for the specified entries. 
#For example, l_sum() sums up the appropriate entries for L[i][k].
def l_sum(i,k):
    l_sum = 0
    for j in range(0,k):
        l_sum += L[i][j]*U[j][k]
    return l_sum

def u_sum(k,j):
    u_sum = 0
    for i in range(0,k):
        u_sum += L[k][i]*U[i][j]
    return u_sum

def z_sum(i):
    z_sum = 0
    for k in range(0,i):
        z_sum += L[i][k]*z[k][0]
    return z_sum

def x_sum(i):
    x_sum = 0;
    for k in range(N-i+1,N):
        x_sum += U[N-i][k]*x[k][0]
    return x_sum



#This section calculates the unkown L and U entries using the known values and
#sum functions, then prints the L and U matrices. 
for k in range(1,N):
    for z in range(0,N):
        L[z][k] = A[z][k] - l_sum(z,k)
        if(L[k][k] == 0):
            U[k][z] = 0
        else:
            U[k][z] = (A[k][z] - u_sum(k,z))/L[k][k]
#print('L = \n',L,'\n')
#print('U = \n',U,'\n')

#This section calculates the z matrix using the sum function and b matrix.
z = np.zeros((N,1))
for i in range(0,N):
    z[i][0] = (b[i][0] - z_sum(i))/L[i][i]   
#print('z = \n',z,'\n')

#This section calculates the x values (solution) using the sum function and
#z matrix and prints the result.
x = np.zeros((N,1))
x[N-1][0]=z[N-1][0]
for i in range(2,N+1):
    x[N-i][0] = z[N-i][0] - x_sum(i)
#print('x = \n',x,'\n')

#This section checks the solution by multiplying the A matrix with the x matrix
r = np.zeros((N,1))
for j in range(0,N):
    for i in range(0,N):
        r[j][0] += x[i][0]*A[j][i] 
#print('solution check: \n',r)

##############################################################################
        
x0 = [[y_lower]]
xN = [[y_upper]]
T3 = np.vstack((x0,x))
T3 = np.vstack((T3,xN))
T3.resize(1,11)
T = pd.DataFrame(T3)
T.columns = ['y(0)','y(.1)','y(.2)','y(.3)','y(.4)','y(.5)','y(.6)','y(.7)','y(.8)','y(.9)','y(1)']
print('\nLU DECOMP Method')
print('\n',T)

print('\nAfter testing several values for alpha, I found that alpha=.75 was the',
      'the most efficient value which resulted in only 11 iterations to converge',
      'The SOR method seems to be more accurate than the Gauss-Seidel method',
      'but it is more difficult to tell how the accuracy compares to the LU',
      'Decomp Method.')






