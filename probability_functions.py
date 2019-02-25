import numpy as np

def piecewise(x):
    a = -1/0.35
    b = 1+1/3.5
    if x<=0.1:
        return 1
    elif x>=0.45:
        return 0
    else:
        return a*x+b
        
def logistic(x):
    a = 12.6
    b = 20
    c = 8.5
    denom = 1+a*np.e**(b*x-c)
    y = 1/denom
    if y<0:
        return 0
    else:
        return y
        
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