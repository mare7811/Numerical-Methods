# -*- coding: utf-8 -*-
"""
Miguel Mares
Engr 428
Homework 
Due 3/4/2020
"""
#This program uses the Monte Carlo method to estimate the integral of the 
#given function f(x)=sin(x) from 0 to pi including the standard deviation.
#The random numbers and final answer are outputed.
#N is chosen so that the answer is fairly accurate and the program is ran
#3 times with different random numbers each time.
#
#I was not able to figure how to use importance sampling because I could not
#rearrange the y function to get x in terms of y.

import numpy as np
import random

N = 100

a = 0
b = np.pi

x = np.zeros(N)
for i in range(0,N):
    x[i] = random.randint(0,314)/100

def f(x):
    f = np.sin(x)
    return f

def sum2(N):
    sum2 = 0
    for i in range(0,N):
        sum2 += f(x[i])**2
    return sum2

def SD(N):
    SD = np.sqrt((((1/N)*sum2(N))-((1/N)*sum(N))**2)/(N-1))
    return SD

def sum(N):
    sum = 0
    for i in range(0,N):
        sum += f(x[i])
    return sum

monte = (b-a)*((1/N)*sum(N)+SD(N))

print(x)
print('solution = ',monte)