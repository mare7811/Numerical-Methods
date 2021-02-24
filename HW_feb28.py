# -*- coding: utf-8 -*-
"""
Miguel Mares
Engr 428
Homework
Due February 2/28/2020
"""
#This program evaluates the integral of e^(-x^2) from 0 to 1 using quadratures 
#with N = 2,3,4,5. First the integral is put into the form with boundaries from
#-1 to 1 with the change of variable x=(y+1)/2. This gives the function 
#f(y)=.5*e^((-(y+1)/2)^2) integrated from -1 to 1.


import numpy as np

#N values 
N = [2,3,4,5]

#Table x_m values for each N value
x = [[.5773502691896258,-.5773502691896258],
     [.77445966692414834,0,-.77445966692414834],
     [.8611363115940526,.3399810435848563,-.8611363115940526,-.3399810435848563],
     [.9061798459386640,.5384693101056831,0,-.9061798459386640,-.5384693101056831]]

#Table W_m value for each N value
W = [[1,-1],
     [.55555555555555556,.8888888888888889,-.55555555555555556],
     [.3478548451374539,.6521451548625461,-.3478548451374539,-.6521451548625461],
     [.2369268850561891,.4786286704993665,.5688888888888889,-.2369268850561891,-.4786286704993665]]

#Definition of function
def f(y):
    f = .5*np.exp((-(y+1)/2)**2)
    return f

#Calculate sum
for i in range(0,len(N)):
    solution = 0
    for m in range(0,N[i]):
        solution += W[i][m]*f(x[i][m])
    print('solution for N =',N[i],':',solution)

print('Actual answer: .746824')