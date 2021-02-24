# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 23:47:22 2020

@author: Newuser
"""

#######################################################################
#Miguel Mares
#
#Numerical Methods
#
#Homework 3
#
#Jan 24, 2020
#######################################################################


import numpy as np    
import matplotlib.pyplot as plt



a = .3
m = .511
v_0 = 10
h = .076199682

z_0 = np.sqrt((2*m/h**2)*v_0*a**2)

steps = .01
Z = np.arange(-5,10,steps)

x = np.tan(Z)
y = np.sqrt(z_0**2-Z**2)/Z


myFigSize = (12,12)
plt.figure(figsize=myFigSize)

plt.subplot(1,1,1)
plt.plot(Z,x)
plt.plot(Z,y)
plt.ylim(0,50)
plt.xlim(0,10)
plt.grid(True)
plt.ylabel('f(E)')
plt.xlabel('E')
plt.title('Plot of f(E) vs E')
