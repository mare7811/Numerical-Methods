############################################################
#Miguel Mares
#Engr 428
#2/4/2020
#Exercise 3.16
############################################################
#This program finds the least squares best linear fit for 
#the data in Table 3.2 in the textbook using LU decomposition
#to solve for the minimum solution.
#############################################################

import matplotlib.pyplot as plt
import numpy as np

#Table 3.2 entries
f = [.242,.616,1,1.881,11.86,29.33]
t = [.388,.724,1,1.524,5.2,9.51]
N = 6   #number of entries

A = np.zeros((2,2))

#Calculate sum of t_i and t_i^2
ti_sum = 0
for i in range(0,6):
    ti_sum += t[i]

ti2_sum = 0
for i in range(0,6):
    ti2_sum += t[i]**2

#Set appropriate matrix entries
A[0][0] = N
A[0][1] = ti_sum
A[1][0] = ti_sum
A[1][1] = ti2_sum

print('A = ', '\n', A, '\n')

#Calculate sum of f_i and f_i*t_i
fi_sum = 0
for i in range(0,6):
    fi_sum += f[i]
    
fi_ti_sum = 0
for i in range(0,6):
    fi_ti_sum += f[i]*t[i]

#set r entries
r = np.zeros((2,1))
r[0][0] = fi_sum
r[1][0] = fi_ti_sum

print('r = ', '\n', r, '\n')

#Create L matrix and set appropriate entries as 1's.
L = np.zeros((2,2))
L[0][0] = 1
L[1][1] = 1


#Create U matrix and set appropriate known entries.
U = np.zeros((2,2))
U[0][1] = A[0][1]   #c1
U[0][0] = A[0][0]   #B1

N = 1   #Number of diagonal entries minus 1. ()

#Calculate l's and B's.
for i in range(0,N):
    L[i+1][i] = A[i+1][i]/U[i][i]                              #l2 to lN+1
    U[i+1][i+1] = A[i+1][i+1]-(L[i+1][i]*U[i][i+1])            #B2 to BN+1

#print (L,'\n')
#print (U)

#Create rho matrix
p = np.zeros((4,1))
p[0][0] = r[0][0]

#Calculate rho values p1 to pN
for i in range(1,N+1):
    p[i][0] = r[i][0] - L[i][i-1]*p[i-1][0]
    
#print ('\n',p)

#Create x matrix
x = np.zeros((2,1))
x[N][0] = p[N][0]/U[N][N]

#Calculate x's
for i in range(N-1,-1,-1):
    x[i][0] = (p[i][0]-U[i][i+1]*x[i+1][0])/U[i][i]
print ('x = (a,b) ','\n', x,'\n')


r_check = np.zeros((2,1))

#Check answers
r_check[0][0] = A[0][0]*x[0][0] + A[0][1]*x[1][0]
r_check[1][0] = A[1][0]*x[0][0] + A[1][1]*x[1][0]
#print(r_check)

z = np.arange(0,10+.01,.01)

#Create linear fit equation using calculated a and b values.
y = x[0][0] + x[1][0]*z

#Plot data and linear fit line
myFigSize = (12,12)
plt.figure(figsize=myFigSize)
plt.subplot(1,1,1)
plt.scatter(t,f)
plt.plot(z,y)
plt.grid(True)
plt.ylabel('Period of Orbit [years]')
plt.xlabel('Relative Distance [AU]')
plt.title('Best Linear Fit')
plt.show()

print('The slope of this line is 3.176 and this should be the value of n. \n'
      'For the result to be consistent with Keplers Third Law, n should equal 1.5 \n'
      'The line is plotted by y=a+bz.')
