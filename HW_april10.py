"""

Miguel Mares
Numerical Methods
Homework
Due 4/10/2020
"""
#This function uses the finite difference method to solve the problem:
#y'' + 2y' + (k^2)y = 0, y(0) = 0 and y(1) = 0 for acceptable values of k.
#First different lambda values are swept from -40 to 0, and a determinant
#for each lambda tested is found from the matrix (A-lI). These values are used
#to find lambda's and determinants that are close to a root. The secant method

#is then used to on those values to find a root (lambda). This lambda is 
#where the determinant is zero and the k can be found from these. The k that
#is found is an allowable value of k.
import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt


#secant method function accepts two lambda values and two determinant values
#and then finds the closest root at those entries.
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
            A_lambda(l_new)
            y.append(determinant())
        else:
            r.append(0)
        i += 1
    return(r[i-2])

#This functions finds the determinant of the (A-lI) matrix
def determinant():
    global A_l
    global A
    Det = np.zeros(N)
    Det[0] = A_l[0][0]
    Det[1] = A_l[1][1]*Det[0] - A_l[1][0]*A_l[0][1]
    for i in range(2,N):
        Det[i] = A_l[i][i]*Det[i-1] - A_l[i][i-1]*A_l[i-1][i]*Det[i-2]
    return Det[N-1]

#This function updates the (A-lI) matrix to be used by the determinant() function
def A_lambda(l):
    global A_l
    global A
    for i in range(0,N):
        A_l[i][i] = A[i][i] - l
    
#specifying step size and boundary conditions
h = .2
a = 0
b = 1
y_upper = 0
y_lower = 0

#Number of intervals, minus 1
N = math.floor((b-a)/h) - 1

y_a = (1/h**2) + (1/h)
y_b = -2/(h**2)                                                                                                   
y_c = (1/h**2) - (1/h)

#matrix, A, from second order ODE, finite difference equation
A = np.zeros((N,N))
A[0][0] = y_b
A[0][1] = y_c
A[N-1][N-2] = y_a
A[N-1][N-1] = y_b
for i in range(1,N-1):
    A[i][i-1] = y_a
    A[i][i] = y_b
    A[i][i+1] = y_c

#This section sweeps through different lambda values from -40 to 0 and plots results
M = 40
l = np.zeros(M)
D = np.zeros(M)
for i in range(-M,0):
    l[i] = i
    A_lambda(l[i])
    D[i] = determinant()

#Plotting
myFigSize = (15,15)
plt.figure(figsize=myFigSize)
plt.subplot(1,1,1)
plt.plot(l,D)
plt.grid(True)
plt.show()

#This section creates a table of lambda and determinant values for outputting
l.resize((l.size,1))
D.resize((D.size,1))
T = np.concatenate((l,D),axis=1)
T = pd.DataFrame(T)
T.columns = ['lambda','determinant']
#print(T.to_string(index=False))

#These statements use the secant() function to find the root, given the lambda and determinant values
root1 = secant(l[30],l[31],D[30],D[31])
print('root1 = ',root1)

root2 = secant(l[5],l[6],D[5],D[6])
print('root2 = ',root2)

#finds the appropriate k for the root that is found (lambda = -k^2)
k1 = np.sqrt(-root1)
k2 = np.sqrt(-root2)
print('k1 = ',k1)
print('k2 = ',k2,'\n')

print('These k\'s can be plugged into the original problem and then be solved',
      'using previous methods.')
print('This method can be extended to find more roots and thus more acceptable',
      'values of k')





      
    