"""
probability filters that compute associated probability for properties of model
and observation SEDs and model images
"""
import numpy as np
from scipy import optimize

### SED Probability Filters ###
def P1(x):
    """
    Arctan filter used to test brightness ratio at 0.75um.
    Designed so Fdisk/F*<=0.1 has P=1 and Fdisk/F*>=0.25 has P=0.
    Arg:
        x(float):Fdisk/F* output of brightness_test1
    Returns:
        y(float): associated probability
    """
#     a = 26
#     b = 99
#     c = 20
#     denom = 1+a*np.e**(b*x-c)
#     y=1/denom
#     return y  

    x0,x1 = [0.1,0.25]; xvals = [x0,x1]; yvals = [1,0]
    a,b = np.polyfit(xvals,yvals,deg=1)
    conds = [x < x0, (x0<=x)&(x<=x1), x > x1]
    funcs = 1, lambda x:a*x+b, 0
    y = np.piecewise(x, conds, funcs)
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
#     a = 12.6
#     b = 17.6
#     c = 37.6
#     denom = 1+a*np.e**(b*x-c)
#     y = 1-1/denom
#     return y
    x0,x1 = [1,2]; xvals = [x0,x1]; yvals = [0,1] #prev. x from 1.5>2.5
    a,b = np.polyfit(xvals,yvals,deg=1)
    conds = [x < x0, (x0<=x)&(x<=x1), x > x1]
    funcs = 0, lambda x:a*x+b, 1
    y = np.piecewise(x, conds, funcs)
    return y
    
def Pcolor(x):
    """
    Piecewise logistic filter used to test Class II slope at NIR 2-8um.
    P=0 alpha>0, P=1 for alpha<-0.05
    Arg:
        x(float):alpha computed in color_test
    Returns:
        y(float) associated probability
    """
#     a = 1
#     b = 1
#     c = 27
#     d = 0.25
#     num = a
#     denom = 1+b*np.e**(c*(x+d))
#     y=num/denom
#     return y

    x0,x1 = [-0.5,0]; xvals = [x0,x1]; yvals = [1,0]
    a,b = np.polyfit(xvals,yvals,deg=1)
    conds = [x < x0, (x0<=x)&(x<=x1), x > x1]
    funcs = 1, lambda x:a*x+b, 0
    y = np.piecewise(x, conds, funcs)
    return y
    
### Image Functions and Probability Filters ###

# these functions are needed to produce convolved images
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
    
# these are functions needed to perform second image test
def gaussian(height, center_x, center_y, width_x, width_y):
    """Returns a gaussian function with the given parameters"""
    width_x = float(width_x)
    width_y = float(width_y)
    return lambda x,y: height*np.exp(
                -(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)

def moments(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution by calculating its
    moments """
    total = data.sum()
    X, Y = np.indices(data.shape)
    x = (X*data).sum()/total
    y = (Y*data).sum()/total
    col = data[:, int(y)]
    width_x = np.sqrt(np.abs((np.arange(col.size)-y)**2*col).sum()/col.sum())
    row = data[int(x), :]
    width_y = np.sqrt(np.abs((np.arange(row.size)-x)**2*row).sum()/row.sum())
    height = data.max()
    return height, x, y, width_x, width_y

def fitgaussian(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution found by a fit"""
    params = moments(data)
    errorfunction = lambda p: np.ravel(gaussian(*p)(*np.indices(data.shape)) -
                                 data)
    p, success = optimize.leastsq(errorfunction, params)
    return p

# these functions are needed to compute associated probabilities 
def image_P1(x):
    """
    Logistic filter used to test used to test shape of star/disk,
    or ratio of sma_y to PSF sma_y. sma ratio 1.5< corresponds to P=0,
    and sma ratio>2 corresponds to P=1
    Arg:
        x(float):sma_ratio output of image_shape_test
    Returns:
        y(float): associated probability
    """
#     a = 1
#     b = 0.1
#     # old parameters (1.5-2)
#     c = 23
#     d = -1.8
#     # new parameters (1.25-1.5)
# #    c = 40
# #    d = -1.4
#     denom = 1 + b*np.e**(-c*(x+d))
#     y = a/denom
#     return y
    x0,x1 = [1.5,2]; xvals = [x0,x1]; yvals = [0,1]
    a,b = np.polyfit(xvals,yvals,deg=1)
    conds = [x < x0, (x0<=x)&(x<=x1), x > x1]
    funcs = 0, lambda x:a*x+b, 1
    y = np.piecewise(x, conds, funcs)
    return y
       
def image_P2(x):
    """
    Logistic filter used to test shape of star/disk,
    or ratio of sma to PSF sma. sma ratio <2 corresponds to P=0,
    and sma ratio>3 corresponds to P=1
    Arg:
        x(float):sma_ratio output of image_shape_test
    Returns:
        y(float): associated probability
    """
    a = 1
    b = 0.3
    c = 9.8
    d = -2.6
    denom = 1 + b*np.e**(-c*(x+d))
    y = a/denom
    return y

def image_P3(x):
    """
    *used when pixel shift=6*
    Logistic filter used to test flux ratio of disk/star,
    flux ratio <0.1 corresponds to P=0, and sma ratio>0.2 corresponds to P=1
    Arg:
        x(float):sma_ratio output of image_shape_test
    Returns:
        y(float): associated probability
    """
#     a = 1
#     b = 0.1
#     c = 110
#     d = -0.17
#     denom = 1 + b*np.e**(-c*(x+d))
#     y = a/denom
#     return y
    x0,x1 = [0.1,0.2]; xvals = [x0,x1]; yvals = [0,1]
    a,b = np.polyfit(xvals,yvals,deg=1)
    conds = [x < x0, (x0<=x)&(x<=x1), x > x1]
    funcs = 0, lambda x:a*x+b, 1
    y = np.piecewise(x, conds, funcs)
    return y
    
def image_P4(x):
    """
    *used when pixel shift=8*
    Logistic filter used to test flux ratio of disk/star,
    flux ratio <0.02 corresponds to P=0, and sma ratio>0.1 corresponds to P=1
    Arg:
        x(float):sma_ratio output of image_shape_test
    Returns:
        y(float): associated probability
    """
#     a = 1
#     b = 0.1
#     c = 150
#     d = -0.08
#     denom = 1 + b*np.e**(-c*(x+d))
#     y = a/denom
#     return y
    x0,x1 = [0.02,0.1]; xvals = [x0,x1]; yvals = [0,1]
    a,b = np.polyfit(xvals,yvals,deg=1)
    conds = [x < x0, (x0<=x)&(x<=x1), x > x1]
    funcs = 0, lambda x:a*x+b, 1
    y = np.piecewise(x, conds, funcs)
    return y

        

        