





import numpy as np
import matplotlib.pyplot as plt

##########################################################################
#The following section uses the given data to create an NxN matrix (A) and 
#Nx1 matrix (b) with least squares fit (either linear, quadratic or cubic) 
#based on the specified N where N = 2: linear, N = 3: quadratic,N = 4: cubic.
##########################################################################


N = 4       #For a nxn matrix, N=n
N_entries = 11
t_d = np.arange(0,1.1,.1)
y_d = [1.67203,1.79792,2.37791,2.66408,2.11245,2.43969,1.88843,1.59447,1.79643,1.07810,.21066]

def A_sums(m):
    ti_sum = 0
    for i in range(0,N_entries):
        ti_sum += t_d[i]**m
    return ti_sum

def b_sums(m):
    fi_ti_sum = 0 
    for i in range(0,N_entries):
        fi_ti_sum += y_d[i]*(t_d[i]**m)
    return fi_ti_sum

A = np.zeros((N,N))
for j in range(0,N):
    for i in range(0,N):
        A[i][j] = A_sums(i+j)
A[0][0] = N_entries
A = 10*A
print('A = \n',A,'\n')


b = np.zeros((N,1))
for i in range(0,N):
    b[i][0] = b_sums(i)
b = 10*b
print('b = \n',b,'\n')

#########################################################################
#The following section uses LU decompositon to solve the matrix A to satisfy 
#the b matrix. It first calculates and outputs the L and U matrix then 
#calculates the x matrix. The solved x matrix elements are the coefficents
#for the linear, quadratic, or cubic polynomial that is the best fit line
#for the data.
#########################################################################

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
print('L = \n',L,'\n')
print('U = \n',U,'\n')

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
print('x = \n',x,'\n')

#This section checks the solution by multiplying the A matrix with the x matrix
r = np.zeros((N,1))
for j in range(0,N):
    for i in range(0,N):
        r[j][0] += x[i][0]*A[j][i] 
print('solution check: \n',r)

###########################################################################
#The following section constructs the best fit line function and plots 
#the data and function on the same plot. It then calculates the error
#for each data entry and plots it on a seperate graph.
###########################################################################
steps = .01
t = np.arange(0,1+steps,steps)
if(N == 2):
    y = x[0][0] + x[1][0]*t
    def y1(t):
        y1 = x[0][0] + x[1][0]*t
        return y1
    print('y = ',x[0][0],'+',x[1][0],'t')

    
if(N==3):
    y = x[0][0] + x[1][0]*t + x[2][0]*t**2
    def y1(t):
        y1 = x[0][0] + x[1][0]*t + x[2][0]*t**2
        return y1
    print('y = ',x[0][0],'+',x[1][0],'t +',x[2][0],'t^2')


if(N==4):
    y = x[0][0] + x[1][0]*t + x[2][0]*t**2 + x[3][0]*t**3
    def y1(t):
        y1 = x[0][0] + x[1][0]*t + x[2][0]*t**2 + x[3][0]*t**3
        return y1
    print('y = ',x[0][0],'+',x[1][0],'t +',x[2][0],'t^2 +',x[3][0],'t^3')

#Plotting polynomial and given data
myFigSize = (12,12)
plt.figure(figsize=myFigSize)
plt.subplot(1,1,1)
plt.scatter(t_d,y_d)
plt.plot(t,y)
plt.grid(True)
plt.ylabel('height [m]')
plt.xlabel('time [s]')
plt.title('Least Squares Fit and Data [quad]')

residual = np.zeros(N_entries)
for i in range(0,N_entries):
    t_i = t_d[i]
    residual[i] = y_d[i] - y1(t_i) 
    
#Plotting the residual error vs time
myFigSize = (12,12)
plt.figure(figsize=myFigSize)
plt.subplot(1,1,1)
plt.scatter(t_d,residual)
plt.grid(True)
plt.ylabel('residual error [m]')
plt.xlabel('time [s]')
plt.title('Residual error vs time[quad]')
