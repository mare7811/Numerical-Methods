# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 19:17:04 2020

@author: Newuser
"""


import numpy as np

def get_A(degree, t_vals):
    
    T_matrix = []
    
    for i in range(0, degree + (degree - 1)):
        
        T_matrix.append((t_vals ** i).sum())
        
    final_matrix = []
    
    for j in range(0, len(T_matrix) - degree + 1):
        row = []
        for k in range(0, degree):
            row.append(T_matrix[j + k])
        final_matrix.append(row)
        
    print("matrix A: ",final_matrix)
    return final_matrix


def get_b(degree, t_vals, f_vals):
    
    final_list = []
    
    for i in range(0, degree):
        entry = (f_vals * (t_vals ** i)).sum()
        final_list.append(entry)
        
    print("vector b: ", final_list)
    return final_list
        
    
        
        
T = np.arange(0, 1.1, 0.1)
temp = [1.67203, 1.79792, 2.37791, 2.66408, 2.11245, 2.43969, 1.88843, 1.59447, 1.79634, 1.07810, 0.21066]
F = np.array(temp)

get_A(2,T)
get_b(2,T,F)