
import numpy as np
import math
import matplotlib.pyplot as plt
import pandas as pd

#############################################################################
#Defining Variables
#
#This section defines the constants as global variables to be used by the
#functions in the program.
#
#input:
#       none
#
#ouput:
#       none
#
#global:
#       h - step size for RK4 method and defining x-axis plot points
#       mu0 - average mass density [kg/m]
#       T - tension of piano wire [N]
#       L - length of piano wire [m]
#       delta - variation of mass density per unit length [kg/m^2]
#       a - lower boundary limit
#       b - upper boundary limit
#       y_lower - lower boundary value
#       y_upper - upper boundary value
#       N - Number of intervals w/ the given step size on boundary limits
#       x - array to hold string position values
#       S1 - array to hold string height values
#       S2 - array to hold string velocity values
#       x_0 - initial string position
#       S1_0 - initial string height (at x=0)
#       S2_0 - initial string velocity - this is an arbitrary number, and not
#              necessarily correct
#       omega - analytical root from textbook
#
############################################################################
h = .05 #step size
mu0 = .000954   #average mass density in kg/m
T = 1000    #Tension in Newtons
L = 1   #Length of string
delta = .0005   #
omega = (np.pi/L)*np.sqrt(T/mu0)    #(eq. 5.135) from textbook - linear string 


# Interval limits, and number of intervals
a = 0
b = L
y_lower = 0
y_upper = 0
N = math.floor((b-a)/h)

#Defining arrays to hold values of RK4 results and steps
x = np.zeros((N+1,1))
S1 = np.zeros((N+1,1))
S2 = np.zeros((N+1,1))

#Defining initial conditions
x_0 = a
S1_0 = y_lower
S2_0 = 3


#############################################################################
#Defining Functions
#
#This section defines the functions that are used to solve the piano wire 
#problem.
#
#input:
#       x1 - first guess point for secant method
#       x2 - second guess point for secant method
#       y1 - integral point from RK4 method using first guess point
#       y2 - integral point from RK4 method using second guess point
#
#ouput:
#       u - calculated mass density at given x point
#       S1[N] - final point from evaluating DE using RK4 method
#       r[i-2] - calculated root from secant method
#
#local:
#       f01 - value used in RK4 algorithm calulation
#       f02 - value used in RK4 algorithm calculation
#       f11 - value used in RK4 algorithm calculation
#       f12 - value used in RK4 algorithm calculation
#       f21 - value used in RK4 algorithm calculation
#       f22 - value used in RK4 algorithm calculation
#       f31 - value used in RK4 algorithm calculation
#       f32 - value used in RK4 algorithm calculation
#       r - array that holds roots
#       y - array that holds secant points
#       i - iteration counter       
#       w - new calculated w from secant algorithm
#
#functions:
#       u() - Function that describes the mass density of the piano string at
#              at each x position
#
#       F1() - Function describing derivative of y(t) 
#       
#       F2() - Function describing second derivative of y(t)
#
#       RK4() - Function that performs the RK4 method algorithm and returns the 
#               final point of the integral.
#
#       secant() - Function that performs the secant algorithm and finds the 
#                  root with the given initial points.
#       
############################################################################
def u(x):
    u = mu0 + (x-(L/2))*delta
    return u

#Defining functions for RK4 method
def F1(x,S1,S2):
    return S2

def F2(x,S1,S2):
    f2 = -((w**2)*u(x)/T)*S1
    return f2

#RK4 method - returns final value of solved second order ODE 
def RK4():
    global S1
    S2 = np.zeros((N+1,1))
    S1 = np.zeros((N+1,1))

    #Filling arrays with initial values
    for i in range(0,N+1):
        x[i] = x_0 + i*h
    S1[0] = S1_0
    S2[0] = S2_0

    #Calculating/filling in rest of table values
    for i in range(1,N+1):
        f01 = F1(x[i-1],S1[i-1],S2[i-1])
        f02 = F2(x[i-1],S1[i-1],S2[i-1])        
        f11 = F1((x[i-1]+(h/2)),(S1[i-1]+(h/2)*f01),(S2[i-1]+(h/2)*f02))
        f12 = F2((x[i-1]+(h/2)),(S1[i-1]+(h/2)*f01),(S2[i-1]+(h/2)*f02))        
        f21 = F1((x[i-1]+(h/2)),(S1[i-1]+(h/2)*f11),(S2[i-1]+(h/2)*f12))
        f22 = F2((x[i-1]+(h/2)),(S1[i-1]+(h/2)*f11),(S2[i-1]+(h/2)*f12))    
        f31 = F1((x[i-1]+h),(S1[i-1]+h*f21),(S2[i-1]+h*f22))
        f32 = F2((x[i-1]+h),(S1[i-1]+h*f21),(S2[i-1]+h*f22))
    
        S1[i] = (S1[i-1] + (h/6)*(f01 + (2*f11) + (2*f21) + f31)) 
        S2[i] = (S2[i-1] + (h/6)*(f02 + (2*f12) + (2*f22) + f32))
    return(S1[N])

#Function that finds root using secant method
def secant(x1,x2,y1,y2):
    global w
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
        if((y[i-1]-y[i-2]) != 0 ):
            w = (r[i-2]*y[i-1] - r[i-1]*y[i-2])/(y[i-1]-y[i-2])
            r.append(w)
            y.append(RK4())
        else:
            r.append(0)
        i += 1
    return(r[i-2])


#############################################################################
#Shooting method
#
#This section uses the shooting method to solve for the fundamental frequency
#of the "piano wire" presented in the textbook. 
#
#input:
#       none
#
#ouput:
#       none
#
#local:
#       w_guess1 - first guess for fundamental frequency
#       w_guess2 - second guess for fundamental frequency
#       g1 - end point of integral for first guess 
#       g2 - end point of integral for second guess
#       root - calculated root from secant method
#       omega - analytical root from textbook
#       y_linear - analytical function of linear string from textbook
#       y_s - function of nonlinear string using shooting method
#       T -  2 dimensional array of x, linear y(x), and nonlinear y(x) values
#       p_dif - percent difference between linear and nonlinear fundamental 
#               frequency
#
############################################################################
#initial fundamental frequency guesses
w_guess1 = 2000
w_guess2 = 3000

w = w_guess1    #sets global variable to initial guess 1
g1 = RK4()      #end point of integral from RK4 method w/ initial guess 1

w = w_guess2    #sets global variable to initial guess 2
g2 = RK4()      #end point of integral from RK4 method w/ initial guess 2
 
#Finding root, using initial w guesses and corresponding RK4 result
root = secant(w_guess1,w_guess2,g1,g2)

################## String equations - y(x) ###################################
y_linear = np.sin(omega*np.sqrt(mu0/T)*x)   #(eq 5.132)
y_s = np.sin(root*np.sqrt(mu0/T)*x)

p_dif = (np.abs(omega-w)/omega)*100

print('Shooting Method')
print('[nonlinear]  w =', root)
print('[linear]  w =', omega)
print('Percent difference:',p_dif,'%')

########### Plotting #########################################################
myFigSize = (15,15)
plt.figure(figsize=myFigSize)

plt.subplot(2,1,1)
plt.plot(x,y_s)
plt.plot(x,y_linear)
plt.grid(True)
plt.ylabel('y(x) [m]')
plt.xlabel('x [m]')
plt.title('Plot of strings')
plt.show()

######################### Outputting Table ####################################
x.resize((x.size,1))
y_s.resize((y_s.size,1))
y_linear.resize((y_linear.size,1))
T = np.concatenate((x,y_s,y_linear),axis=1)
T = pd.DataFrame(T)
T.columns = ['x', 'y(x) [nonlinear]','y(x) [linear]']
print(T.to_string(index=False))

