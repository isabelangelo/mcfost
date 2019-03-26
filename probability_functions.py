"""
probability filters that compute associated probability for properties of model
and observation SEDs
"""
import numpy as np

def P1(x):
    """
    Arctan filter used to test brightness ratio at 0.75um.
    Designed so Fdisk/F*<=0.1 has P=1 and Fdisk/F*>=0.5 has P=0.
    Arg:
        x(float):Fdisk/F* output of brightness_test1
    Returns:
        y(float): associated probability
    """
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
    """
    Arctan filter used to test brightness ratio at 4.5um.
    Designed so Fdisk/F*<=0.3 has P=1 and Fdisk/F*>=0.6 has P=0
    Arg:
        x(float):Fdisk/F* output of brightness_test2
    Returns:
        y(float): associated probability
    """
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
    """
    Logisitic filter used to test slope increase of SED between 2-10um.
    delta_m<1.5 have P=0 and delta_m>=2.5 have P=1
    Arg:
        x(float):m2-m1 output of slope_test
    Returns:
        y(float): associated probability
    """
    a = 12.6
    b = 17.6
    c = 37.6
    denom = 1+a*np.e**(b*x-c)
    y = 1-1/denom
    return y
    
def gaussian_kernel():
    """ 
    generates 2D gaussian kernel to convolve with mcfost image
    """
    
    # generate 1D gaussian
    t = np.linspace(-5, 5, 15)
    bump = np.exp(-0.75*t**2)
    bump /= np.trapz(bump) # normalize
    
    # make 2D kernel
    kernel = bump[:, np.newaxis] * bump[np.newaxis, :]
    
    return kernel


def rebin(arr, new_shape):
    """
    rebin an array into a desired shape
    """
    shape = (new_shape[0], arr.shape[0] // new_shape[0],
             new_shape[1], arr.shape[1] // new_shape[1])
    return arr.reshape(shape).mean(-1).mean(1)
        

        