
# coding: utf-8

# Alexander Hebert
# ECE 6390
# Matrix Sign Function


# Tested using Python v3.4

import numpy as np

def msf(A,eps,maxiter):
    """msf: Matrix Sign Function
    A matrix is assumed square and invertible
    eps = threshold value for iterations
    maxiter = maximum number of iterations
    
    returns A_j+1 and value of j+1
    flag = 1 if eps breaks loop
    flag = 0 if max iterations is reached"""
    
    A_j = A  # A_0 = A
    
    for j in range(0,maxiter):
        
        A_jp1 = 0.5*(A_j + np.linalg.inv(A_j))
        
        if (np.abs(np.trace(A_jp1) - np.trace(A_j)) < eps):
            flag = 1
            jp1 = j + 1
            return A_jp1, jp1, flag
            
        A_j = A_jp1
        
    flag = 0
    jp1 = j + 1
    return A_jp1, jp1, flag

