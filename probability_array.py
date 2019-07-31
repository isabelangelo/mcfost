"""
creates arrays to be used for statistical analysis of MCFOST grid.

1. generates an array that contains the associated probability for each grid model
array dimensions:
    
    dust_mass x Rc x f_exp x H0 x Rin x sd_exp x amax

so index of each dimension corresponds to the index of that parameter in 
grid_parameters dictions. For example, 

arr[0,0,0,0,0,0,0] corresponds to model_grid.get_model_index('1e-07,10,0.85,5,0.1,0,10')

2. generates masked array where excluded models have their value set to zero.

examples: arr[0] will be 0 if we exclude dust_mass=1e-7, 
         arr[:,0]=0 if we exclude Rc=10
    
written: Isabel Angelo (2019)
"""
import numpy as np
from astropy.io import fits
from generate_para import dust_mass, Rc, f_exp, H0, Rin, sd_exp, amax
from model_grid import *
from models import *
from test_edge_on import *

# define shape of grid arrays
if grid == 'new':
    arr_shape = (5,4,4,4,3,4,4,15)
if grid == 'old':
    arr_shape = (5,4,4,4,3,3,4,15)

def generate_probability_array():
    """
    generate hypercube of models where the value in the array corresponds to that
    model's probability of being edge-on
    
    example:
        
        arr[0,0,0,0,0,0,0] is the probability of Model(0) being edge-on
        
    Returns:
        probability_array.fits (FITS):file containing hypercube
    """
    # generate empty array to fill with probabilities
    arr = np.empty(arr_shape)
    for a in range(len(dust_mass)):
        a_s = str(dust_mass[a])
        for b in range(len(Rc)):
            b_s = str(Rc[b])
            for c in range(len(f_exp)):
                c_s = str(f_exp[c])
                for d in range(len(H0)):
                    d_s = str(H0[d])
                    for e in range(len(Rin)):
                        e_s = str(Rin[e])
                        for f in range(len(sd_exp)):
                            f_s = str(sd_exp[f])
                            for g in range(len(amax)):
                                g_s = str(amax[g])
                            
                                # get index with corresponding model                            
                                param_str = ','.join([a_s,b_s,c_s,d_s,e_s,f_s,g_s])
                                idx = get_model_index(param_str)
                                #p = compute_P(Model(idx))
                                p = image_compute_P(Model(idx))
                                arr[a,b,c,d,e,f,g,:] = p
                                print(idx)
                
    # write the large array to a text file
    hdu = fits.PrimaryHDU(data=arr)
    hdu.writeto('image_probability_array.fits')

    
def generate_masked_array(dust_mass=None, Rc=None, f_exp=None, H0=None,\
    Rin=None, sd_exp=None, amax=None):
    """
    Create array with masks on parameters that are NOT input.
    For example, Rc=[10,30,100] will return grid where Rc=300 are masked (set to 0)
    
    Args:
        dust_mass, Rc, etc. (list): list of values for a specified parameter to be
        included in masked array
    
    Returns:
        masked_arr(np.array): hypercube with masks as values
    """
    # generate array of ones to weight
    masked_arr = np.ones(arr_shape)
    # update ones to zeros for parameters that are masked
    if dust_mass is not None:
        for i in range(len(grid_parameters['dust_mass'])):
            if grid_parameters['dust_mass'][i] not in dust_mass:
                masked_arr[i] = np.zeros(masked_arr[i].shape)
    if Rc is not None:
        for i in range(len(grid_parameters['Rc'])):
            if grid_parameters['Rc'][i] not in Rc:
                masked_arr[:,i]=np.zeros(masked_arr[:,i].shape)
    if f_exp is not None:
        for i in range(len(grid_parameters['f_exp'])):
            if grid_parameters['f_exp'][i] not in f_exp:
                masked_arr[:,:,i]=np.zeros(masked_arr[:,:,i].shape)
    if H0 is not None:
        for i in range(len(grid_parameters['H0'])):
            if grid_parameters['H0'][i] not in H0:
                masked_arr[:,:,:,i] = np.zeros(masked_arr[:,:,:,i].shape)
    if Rin is not None:
        for i in range(len(grid_parameters['Rin'])):
            if grid_parameters['Rin'][i] not in Rin:
                masked_arr[:,:,:,:,i] = np.zeros(masked_arr[:,:,:,:,i].shape)
    if sd_exp is not None:
        for i in range(len(grid_parameters['sd_exp'])):
            if grid_parameters['sd_exp'][i] not in sd_exp:
                print(i)
                masked_arr[:,:,:,:,:,i] = np.zeros(masked_arr[:,:,:,:,:,i].shape)
    if amax is not None:
        for i in range(len(grid_parameters['amax'])):
            if grid_parameters['amax'][i] not in amax:
                masked_arr[:,:,:,:,:,:,i] = np.zeros(masked_arr[:,:,:,:,:,:,i].shape)                           
    return masked_arr
    
    
def generate_weighted_array(dust_mass_weights = [1,1,1,1,1],
                        Rc_weights = [1,1,1,1],
                        f_exp_weights = [1,1,1,1],
                        H0_weights = [1,1,1,1],
                        Rin_weights =  [1,1,1],
                        sd_exp_weights = [1,1,1,1],
                        amax_weights = [1,1,1,1]):
    """
    create an array with weights specified by input parameters. For example,
    Rc=[10,30,100,300] can have a corresponding Rc_weights=[0.15,0.2,0.3,0.35],
     with an output where all models with Rc=10 will be weighted by 0.15, 
     Rc=30 by 0.2, and so on
    
    Args:
        dust_mass_weights, etc. (list): list of weights where each index represents the
        weight of the corresponding parameter value
        
    Returns:
        weighted_arr(np.array): hypercube with weights as values
        
    """
    # generate empty array
    weighted_arr = np.ones(arr_shape)
    # weight dust mass
    for i in range(len(dust_mass)):
        weighted_arr[i] *= dust_mass_weights[i]
    # weight Rc
    for i in range(len(Rc)):
        weighted_arr[:,i] *= Rc_weights[i]
    # weight f_exp
    for i in range(len(f_exp)):
        weighted_arr[:,:,i] *= f_exp_weights[i]
    # weight H0
    for i in range(len(H0)):
        weighted_arr[:,:,:,i] *= H0_weights[i]
    # weight Rin
    for i in range(len(Rin)):
        weighted_arr[:,:,:,:,i] *= Rin_weights[i]
    # weight sd_exp
    for i in range(len(sd_exp)):
        weighted_arr[:,:,:,:,:,i] *= sd_exp_weights[i]
    # weight amax
    for i in range(len(amax)):
        weighted_arr[:,:,:,:,:,:,i] *= amax_weights[i]
            
    return weighted_arr
    
