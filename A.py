# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 20:00:48 2020

@author: Newuser
"""

#Given Data    
N_entries = 11        #Number of data entries
t_d = np.arange(0,1.1,.1)
y_d = [1.67203,1.79792,2.37791,2.66408,2.11245,2.43969,1.88843,1.59447,1.79643,1.07810,.21066]

ti_sum = 0
ti2_sum = 0
ti3_sum = 0
ti4_sum = 0
fi_sum = 0
fi_ti_sum = 0
fi_ti2_sum = 0
for i in range(0,N_entries):
    ti_sum += t_d[i]
    ti2_sum += t_d[i]**2
    ti3_sum += t_d[i]**3
    ti4_sum += t_d[i]**4
    fi_sum += y_d[i]
    fi_ti_sum += y_d[i]*t_d[i]
    fi_ti2_sum += y_d[i]*t_d[i]**2

A = np.zeros((N,N))
A[0][0] = N_entries
A[0][1] = ti_sum
A[0][2] = ti2_sum
A[1][0] = ti_sum
A[1][1] = ti2_sum
A[1][2] = ti3_sum
A[2][0] = ti2_sum
A[2][1] = ti3_sum
A[2][2] = ti4_sum
print('A = \n',A,'\n')

b = np.zeros((N,1))
b[0][0] = fi_sum
b[1][0] = fi_ti_sum
b[2][0] = fi_ti2_sum
print('b = \n',b,'\n')
