"""
Defines a set of tests for determining if a Model or Obs object are observably
edge-on/optically thick. Probability for a given SED is generated based on weighted
probabilities computed by 3 tests: brightness1 (0.75um), brightness2 (4.5um),
and slope test. Probability for a given image is generated based on brightest pixel
location, disk shape and outer disk flux.

SED probability functions take a Model or Obs object as input and return
a list containing an associated probability (float) 0<P<1 for all inclinations.
Image probability functions only take models.

written: Isabel Angelo (2019)
"""
from models import *
from observations import *

### SED Tests ###
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
            ## store probability associated with brightness ratio
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
        
def color_test(obj):
    """
    Tests for color of object.
    
    Returns list that contains color for each model inclination
    """
    # compute 0.75 zeroth magnitude flux
    x1,x2 = 6.38e-7*10**6.,7.97e-7*10**6. #micron
    dx1,dx2 = 1.38e-7, 1.49e-7 #meters
    y1,y2 = 2.2387e-02*dx1, 1.1482e-02*dx2 #W/m^2
    m = (y2-y1)/(x2-x1); b = y1-(m*x1) #W/m^2/micron, W/m^2
    F0_75 = m*0.75+b #W/m^2      
    # compute 4.5 zeroth magnitude flux
    F0_45_WmHz = 179.7*(1e-26) #converting from jansky
    F0_45_Wmm = F0_45_WmHz*2.998e8/(4.5e-6)**2. #F_nu *1e-26*c/lambda^2=F_lambda
    F0_45 = F0_45_Wmm*1.1e-6 #W/m^2
    
    # handle models
    if type(obj)==Model:
        #store probabilities
        P_i = []
        for inc_idx in range(len(obj.seds)):
            # compute flux at 0.75 and 4.5 um
            obj_sed = obj.seds[inc_idx]
            F_75 = obj_sed[np.where(obj.wavelength==0.75)][0]
            F_45 = obj_sed[np.where(obj.wavelength==4.5)][0]
            # compute magnitudes
            m_75 = -2.5*np.log10(F_75/F0_75)
            m_45 = -2.5*np.log10(F_45/F0_45) 
            color = m_75-m_45
            P_i.append(color)
            return P_i
            
    # handle observations
    if type(obj)==Obs:
        # determine index of fluxes to compute
        idx_75 = np.abs(obj.wavelength-0.75).argmin()
        idx_45 = np.abs(obj.wavelength-4.5).argmin()
        # compute fluxes
        F_75 = obj.sed[idx_75]
        F_45 = obj.sed[idx_45]
        # compute magnitudes
        m_75 = -2.5*np.log10(F_75/F0_75)
        m_45 = -2.5*np.log10(F_45/F0_45)
        color = m_75-m_45
        return color
        
        

        
def compute_P(obj):
    """
    Compute final weighted probability from 3 tests for each model inclination
    """
    #sum = np.array(brightness_test1(obj))+np.array(brightness_test2(obj))+np.array(slope_test(obj))
    sum = np.array(color_test(obj))
    #return sum/3.
    return sum
    
    
### Image Tests ###
# def image_brightpixel_test(model):
#     """
#     determines binary edge-on probability based on criteria of the pixel location
#     of the brightest pixel in the image
#     """
#     P_i = []
#     images = model.convolved_images
#     # get center pixel of images
#     center_pixel = int(images[0].shape[0]/2-0.5)
#     # determine location of brightest pixel + associated probability
#     for inc_idx in range(len(model.images)):
#         im=images[inc_idx]
#         if im[center_pixel,center_pixel]!=im.max():
#             P_i.append(1)
#         else:
#             P_i.append(0)
#     return P_i   
    
# def image_ysma_test(model):
#     """
#     determines edge-on probability (0<=P<=1) based on ratio of convolved image sma to 
#     TinyTim PSF sma
#     
#     note: width_x corresponds to the vertical axis (rows) and 
#     width_y corresponds to the horizontal axis (columns)
#     """
#     P_i = []
#     # compute sma ratio of image to PSF
#     images=model.convolved_images
#     for inc_idx in range(len(images)):
#         data=images[inc_idx]
#         params = fitgaussian(data)
#         fit = gaussian(*params)
#         (height, x, y, width_x, width_y) = params
#         sma_ratio = width_x/PSFwidth_x
#         # append associated probability
#         P_i.append(image_P1(sma_ratio))
#     return P_i  

def image_verticalfluxratio_test(model):
    """
    determines edge-on probability (0<=P<=1) based on the flux ratio of 
    an offset row brightest pixel to the image brightest pixel
    """
    # define vertical shift for flux ratio computation
    shift = 6
    # define flux ratio threshold
    P_i = []
    images = model.convolved_images
    for inc_idx in range(len(images)):
        # get image for inc
        im = images[inc_idx]
        # store location of maximum
        (rowmax,colmax) = np.unravel_index(np.argmax(im, axis=None), im.shape)
        # brightest pixels in row shifted above + below
        num1 = np.max(im[rowmax-shift,:]);num2 = np.max(im[rowmax+shift,:])
        # compute ratio
        num = max(num1,num2)
        denom = im[rowmax,colmax] # brightest pixel
        flux_ratio = num/denom
        P_i.append(image_P3(flux_ratio))
    return P_i

        
PSFparams = fitgaussian(tinytim_PSF)
PSFfit = gaussian(*PSFparams)
(PSFheight, PSFx, PSFy, PSFwidth_x, PSFwidth_y) = PSFparams


def image_shape_test(model):
    """
    determines edge-on probability (0<=P<=1) based on ratio of convolved image sma to 
    TinyTim PSF sma
    """
    P_i = []
    # compute sma ratio of image to PSF
    images=model.convolved_images
    for inc_idx in range(len(images)):
        data=images[inc_idx]
        params = fitgaussian(data)
        fit = gaussian(*params)
        (height, x, y, width_x, width_y) = params
        sma_ratio = width_y/PSFwidth_y
        # append associated probability
        P_i.append(image_P1(sma_ratio))
    return P_i
    
    
def image_horizontalfluxratio_test(model):
    """
    determines edge-on probability (0<=P<=1) based on the flux ratio of 
    an offset column brightest pixel to the image brightest pixel
    """
    # define horizontal shift for flux ratio computation
    shift = 8 
    # define flux ratio threshold
    P_i = []
    images = model.convolved_images
    for inc_idx in range(len(images)):
        # get image for inc
        im = images[inc_idx]
        # store location of maximum
        (rowmax,colmax) = np.unravel_index(np.argmax(im, axis=None), im.shape)
        # brightest pixels in column shifted left + right
        num1 = np.max(im[:,colmax-shift]); num2 = np.max(im[:,colmax+shift])
        # compute ratio
        num = max(num1,num2)
        denom = im[rowmax,colmax] # brightest pixel
        flux_ratio = num/denom
        P_i.append(image_P4(flux_ratio))
    return P_i
    
def image_compute_P(model):
    """
    Compute final weighted probability from 3 tests for each model image inclination
    """
    sum = np.array(image_verticalfluxratio_test(model))+\
            np.array(image_shape_test(model))+\
            np.array(image_horizontalfluxratio_test(model))
    return sum/3.       
        
        
        
        
