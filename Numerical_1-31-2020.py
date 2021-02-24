############################################################
#Miguel Mares
#Engr 428
#1/30/2020
#LU Factorization
############################################################

import numpy as np

A = [[2,-1, 0,0],    #[b1,c1,0, 0]
     [-1,2,-1,0],    #[a2,b2,c2,0]
     [0,-1,2,-1],    #[0,a3,b3,c3]
     [0, 0,-1,2]]    #[0, 0,a4,b4]

r = [[1],
     [0],
     [0],
     [1]]

#Create matrix and set appropriate entries as 1's.
L = np.zeros((4,4))
L[0][0] = 1
L[1][1] = 1
L[2][2] = 1
L[3][3] = 1

#Create matrix and set appropriate known entries.
U = np.zeros((4,4))
U[0][1] = A[0][1]   #c1
U[1][2] = A[1][2]   #c2    
U[2][3] = A[2][3]   #c3
U[0][0] = A[0][0]   #B1

N = 3   #Number of diagonal entries minus 1. ()

#Calculate l's and B's.
for i in range(0,N):
    L[i+1][i] = A[i+1][i]/U[i][i]                              #l2 to lN+1
    U[i+1][i+1] = A[i+1][i+1]-(L[i+1][i]*U[i][i+1])            #B2 to BN+1

print (L,'\n')
print (U)

#Create p matrix
p = np.zeros((4,1))
p[0][0] = r[0][0]

#Calculate rho values p1 to pN
for i in range(1,N+1):
    p[i][0] = r[i][0] - L[i][i-1]*p[i-1][0]
    
print ('\n',p)

#Create x matrix
x = np.zeros((4,1))
x[N][0] = p[N][0]/U[N][N]

#Calculate x's
for i in range(N-1,-1,-1):
    x[i][0] = (p[i][0]-U[i][i+1]*x[i+1][0])/U[i][i]
    
print ('\n', x)

#Check answers
r_check = A*x
print('\n', r_check)
    