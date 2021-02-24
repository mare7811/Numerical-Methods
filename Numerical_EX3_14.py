############################################################
#Miguel Mares
#Engr 428
#2/2/2020
#Exercise 3.14
############################################################

import numpy as np


#Original function
def f(x):
    f = x*np.exp(x)
    return f


x = 2           #Set x value

#Second derivative 3 point expression
def f_2(h):
    f_2 = (f(x+h)-2*f(x)+f(x-h))/h**2
    return f_2


#Create headings for table
B = ['  i','       h_i','        D1','      D2','       D3','       D4']


N = 4                       #Number of h entries
A = np.zeros((N,6))         #Create Richardson Extrapolation Table
h = .4

A[0][1] = h                 #Set first h entry

#Create index
for i in range(0,N):            
    A[i][0] = i

#Enter h values into appropriate boxes
for i in range(1,N):   
    A[i][1] = A[i-1][1]/2
    
#Calculate D1 values
for i in range(0,N):
    A[i][2] = f_2(A[i][1])
    
#Calculate D2 to DN values
for j in range (3,N+2):
    for i in range (j-2,N):
       A[i][j] = (2**(2*(j-2))*A[i][j-1]-A[i][j-1])/(2**(2*(j-2))-1)
    
print(B)    
print(A)

    