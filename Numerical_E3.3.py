#**************************************************************************#
#Miguel Mares
#Numerical Methods
#January 26, 2020
#Exercis 3.3
#**************************************************************************#

##The following is my attempt at plotting the Lagrange 
##interpolation of the Bessel function table values: plotting the 
##relative intensity as a function of p in increments of 0.1 for p.


###Lagrange - J1(p) values were used###

import numpy as np
import matplotlib.pyplot as plt

steps = 0.1
p = np.arange(.1,10+steps,steps)
p1 = np.arange(0.1,2+steps,steps)
p2 = np.arange(2,3+steps,steps)
p3 = np.arange(4,5+steps,steps)
p4 = np.arange(5,6+steps,steps)
p5 = np.arange(6,7+steps,steps)
p6 = np.arange(7,8+steps,steps)
p7 = np.arange(8,9+steps,steps)
p8 = np.arange(9,10+steps,steps)


x = [0,1,2,3,4,5,6,7,8,9,10]
f = [0,.4400505857,.5767248078,.3390589585,-.066043328,-.3275791376,
     -.2766838581,-.0046828235,.2346464469,.2453117866,.0434727462]

l = 1
J = 0
for j in range(0,3):
    for i in range(0,3):
        if i==j:
            l=l
        else:
            l = l*((p1-x[i])/(x[j]-x[i]))
    J = J + l*f[j]
I1 = ((2*J)/p1)**2


l=1
J=0
for j in range(1,4):
    for i in range(1,4):
        if i==j:
            l=l
        else:
            l = l*((p2-x[i])/(x[j]-x[i]))
    J = J + l*f[j]
I2 = ((2*J)/p2)**2    


l=1
J=0
for j in range(2,5):
    for i in range(2,5):
        if i==j:
            l=l
        else:
            l = l*((p3-x[i])/(x[j]-x[i]))
    J = J + l*f[j]
I3 = ((2*J)/p3)**2


l=1
J=0
for j in range(3,6):
    for i in range(3,6):
        if i==j:
            l=l
        else:
            l = l*((p4-x[i])/(x[j]-x[i]))
    J = J + l*f[j]
I4 = ((2*J)/p4)**2


l=1
J=0
for j in range(4,7):
    for i in range(4,7):
        if i==j:
            l=l
        else:
            l = l*((p5-x[i])/(x[j]-x[i]))            
        J = J + l*f[j]
I5 = ((2*J)/p5)**2


l=1
J=0
for j in range(5,8):
    for i in range(5,8):
        if i==j:
            l=l
        else:
            l = l*((p6-x[i])/(x[j]-x[i]))       
    J = J + l*f[j]
I6 = ((2*J)/p6)**2


l=1
J=0
for j in range(6,9):
    for i in range(6,9):
        if i==j:
            l=l
        else:
            l = l*((p7-x[i])/(x[j]-x[i]))  
    J = J+ l*f[j]
I7 = ((2*J)/p7)**2


l=1
J=0
for j in range(7,10):
    for i in range(7,10):
        if i==j:
            l=l
        else:
            l = l*((p8-x[i])/(x[j]-x[i]))  
    J = J+ l*f[j]
I8 = ((2*J)/p8)**2

I9 = [0,0,0]


I = np.concatenate([I1,I2,I3,I4,I5,I6,I7,I8,I9])


myFigSize = (12,12)
plt.figure(figsize=myFigSize)

plt.subplot(1,1,1)
plt.plot(p,I)
plt.grid(True)
plt.ylabel('I/I_0 ')
plt.xlabel('p')
plt.title('Relative Intensity vs p')




