#######################################################################
#Miguel Mares
#Casey Kleinkopf
#Numerical Methods
#Jan 24, 2020
#Project 1
#######################################################################


import numpy as np    
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

#Set constants
a = .8
v_0 = 10
h2me = .0760


#Calculate z_0 based on constants
z_0 = np.sqrt((2*a**2*v_0)/h2me)
print('z_0 = ',z_0)

#Create inputs of equal step sizes.
steps = .01
E = np.arange(0.01,10+steps,steps)
EV = np.arange(0.01,10+steps,steps)/v_0
F = np.arange(.01,1+steps,steps)

#Define f(E) function
def f(E):
    f = z_0*np.sqrt(1-(E/v_0))*np.cos(z_0*np.sqrt(E/v_0))-z_0*np.sqrt(E/v_0)*np.sin(z_0*np.sqrt(E/v_0))
    return f

#Create entries to plot
z = a*f(E)
y = np.zeros(1000) #Line on x axis


#Plotting a*f(E) vs E
myFigSize = (12,12)
plt.figure(figsize=myFigSize)
plt.subplot(1,1,1)
plt.plot(EV,z)
plt.plot(E,y) #Line on x axis
plt.xlim(0,1)
plt.grid(True)
plt.ylabel('a*f(E)')
plt.xlabel('E/v_0')
plt.title('Plot of af(E) vs E/v_0')
plt.show()


root = fsolve(f,.5)
print(root,'\n')

roots = np.zeros(1001)

k = 0
n = np.zeros(1001)
for a in F:
    z_0 = np.sqrt((2*a**2*v_0)/h2me)
    print('a: ', a)
    roots[k] = fsolve(f,.5)
    print(roots[k])
    k += 1
    print('\n')
    

