"""
Defines a set of tests for determining if a Model or Obs object are observably
edge-on/optically thick. Probability for a given SED is generated based on weighted
probabilities computed by 3 tests: brightness1 (0.75um), brightness2 (4.5um),
and slope test. 

SED probability functions take a Model or Obs object as input and return
a list containing an associated probability (float) 0<P<1 for all inclinations.

Image probability functions only take models.

written: Isabel Angelo (2019)
"""
from models import *
from observations import *

def brightness_test1(obj):
    """
    Tests to see if SED is significantly blocked by dust at 0.75um.
    (most of edge-on disks should pass this test)
    
    Returns list that contains associated P for each model inclination
    """ 
    if type(obj)==Model:
    
        # set wavelength for brightness calculation
        # found by taking peak brightness of star with 1e-12 dust mass
        star_peak = 0.75 
        
        # compute peak brightness of the star
        star_sed = star_seds[0]
        F_star = star_sed[np.where(obj.wavelength==star_peak)][0]
        
        # store probabilities for each inclination
        P_i = []
        for inc_idx in range(len(obj.seds)):
            # compute brightness of SED at peak wavelength
            obj_sed = obj.seds[inc_idx]
            F = obj_sed[np.where(obj.wavelength==star_peak)][0]
            # store probability associated with brightness ratio
            P_i.append(P1(F/F_star))
    
        return P_i
    
    if type(obj)==Obs:
        
        return 1
    
    
def brightness_test2(obj):
    """
    Tests for brightness significantly blocked at 4.5um.
    (single-peaked ones should pass this test)
    
    Returns list that contains associated P for each model inclination
    """ 
    if type(obj)==Model:
        # wavelength where starlight is attenuated
        dim_wavelength = 4.5
        
        # compute corresponding brightness of the star
        star_sed = star_seds[0]
        F_star = star_sed[np.where(obj.wavelength==dim_wavelength)][0]
    
        # store probabilities for each inclination
        P_i = []
        for inc_idx in range(len(obj.seds)):
            # compute brightness of SED at dim wavelength
            obj_sed = obj.seds[inc_idx]
            F = obj_sed[np.where(obj.wavelength==dim_wavelength)][0]
            # store probability associated with brightness ratio
            P_i.append(P2(F/F_star))
        
        return P_i
        
    if type(obj)==Obs:
    
        return 1
    
    
def slope_test(obj):
    """
    Tests to see if slope goes from negative to positive i.e. "double-peaked" shape.
    (two-peaked disks and ones with large Rin should pass)
    
    Returns list that contains associated P for each model inclination
    """
    if type(obj)==Model:
        P_i = []
        for inc_idx in range(len(obj.seds)):
            s = obj.get_slopes(inc_idx)
            windows, slopes = 10**s[0], s[1]
    
            min = slopes[np.where((windows>=2) & (windows<=10))].min()
            max = slopes[np.where((windows>10) & (windows<=20))].max()
    
            dif = max-min
            #print('model'+str(obj.n_model), obj.inclinations[inc_idx], dif, P3(dif))
            P_i.append(P3(dif))
            
        return P_i
    
    if type(obj)==Obs:

        s = obj.get_slopes()
        windows, slopes, slope_errs = 10**s[0], s[1], s[2]
        
        min = slopes[np.where((windows>=2) & (windows<=10))].min()
        max = slopes[np.where((windows>10) & (windows<=20))].max()
        
        dif = max-min
        #print(self.obs_name, dif, P3(dif))
        return P3(dif)
        
def compute_P(obj):
    """
    Compute final weighted probability from 3 tests for each model inclination
    """
    sum = np.array(brightness_test1(obj))+np.array(brightness_test2(obj))+np.array(slope_test(obj))
    return sum/3.
        
        
        
        
        
        
