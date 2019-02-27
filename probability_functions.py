"""
probability filters that compute associated probability for properties of model
and observation SEDs
"""
import numpy as np

def P1(x):
    a = -0.332
    b = 53.1
    c = 15.5
    d = 0.5
    y = a*np.arctan(b*x-c)+d
    if y<0:
        return 0
    else:
        return y
        
def P2(x):
    a = -0.325
    b = 53.1
    c = 26.2
    d = 0.5
    y = a*np.arctan(b*x-c)+d
    if y<0:
        return 0
    else:
        return y
        
def P3(x):
    a = 12.6
    b = 17.6
    c = 37.6
    denom = 1+a*np.e**(b*x-c)
    y = 1-1/denom
    return y
    
        

        