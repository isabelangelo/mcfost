"""
This code defines the Model class for handling model SED's from MCFOST.
The class contains the SEDs, their wavelengths, and methods
for plotting individual and multiple SEDs.

Written: Isabel Angelo (2019)
"""
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats.mstats as stats
from reduce_sed import parseline

###TO DO:#
###make input for static method list instead of just 2?
###but should we only normalize in overplot? For now I like that it defaults to it

# define path to models
model_path = 'grid_i15/' # add ../ for leo

# set up plot axes
def setup_plot():
    """
    puts a plot on a log scale and labels the axes for an SED
    """
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('$\lambda [\mu m]$')
    plt.ylabel(r'$\lambda F_{\lambda} [\frac{W}{m^2}]$')
    

class Model(object):
    """
    Model object
    
    Args:
        n_model(int): integer associated with corresponding MCFOST model
        norm(bool): 'True' will normalize according to the geometric mean 
                        between 2-10um of the highest inclination in the model, 
                        'False' plots raw data
    """    
    def __init__(self, n_model, norm=True):
        
        # retrieve sed and wavelength data from model
        self.n_model = n_model
        self.filepath = model_path + 'model'+str(n_model)+'/data_th/sed_rt.fits'
    
        hdulist = fits.open(self.filepath)
        sedlist = hdulist[0].data[0][0]
        self.wavelength = hdulist[1].data
        
        # normalize according to highest inclination
        if norm==True:
            mean_range = np.where((2<self.wavelength) & (self.wavelength<10))
            gmean = stats.gmean(sedlist[-1][mean_range])
            self.seds = [sed/gmean for sed in sedlist]
            
        # for non-normalized case
        else:
            self.seds = sedlist
        
    def plot(self):
        """
        plots an MCFOST model SED [W/m^2] versus wavelength [um]
        """
        # plot SEDs for all model inclinations
        n_i = len(self.seds)
        for i in range(n_i):                                                                                                       
            plt.plot(self.wavelength, self.seds[i], color='k', linewidth=0.5)
            
        setup_plot()
        plt.xlim(self.wavelength[0],self.wavelength[-1])
        plt.title('model '+str(self.n_model))
       # plt.show()
        
    @staticmethod
    def overplot(model1, model2):
        """
        overplots two MCFOST model SED's with multiple inclinations and creates a 
        corresponding legend
        """
        # plot SEDs for all model inclinations for each model
        n_i = len(model1.seds)
        for i in range(n_i):                                                                                                       
            lines1 = plt.plot(model1.wavelength, model1.seds[i], color='k', linewidth=0.5)
            lines2 = plt.plot(model2.wavelength, model2.seds[i], color='r', linewidth=0.5)
         
        setup_plot()
        # set up legend that only shows both models without individual inclinations
        plt.legend([lines1[-1],lines2[-1]], 
            ['model '+str(model1.n_model), 'model '+str(model2.n_model)])   
            
    def get_slopes(self, window):

        # take log for analysis
        logw = np.log10(self.wavelength)
        logs = np.log10(self.seds[-1])
        window = np.log10(window)
        
        # define bins
        window_start = np.arange(logw[0],logw[-1],0.01)
        window_center = window_start+window/2.
        
        # set up figure subplots
        ax0 = plt.subplot(211)
        ax0.plot(logw,logs,'k.')
        
        # store slopes for each window
        slopes = []       
        for i in window_start:
            fitpoints = np.where((i<=logw)&(logw<=i+window))[0]
            if len(fitpoints)>=2:
                # perform linear fit
                x = logw[fitpoints]
                y = logs[fitpoints]
                a,b=np.polyfit(x,y,1)
                # plot polyfit as a test
                ax0.plot(x,a*x+b, '-')
                # store slope
                slopes.append(a)
            else:
                slopes.append(np.nan)
                
        # store unique slope and corresponding window values (median index)
        unq_slopes = np.unique(slopes) 
        unq_window_median = []
        
        for s in unq_slopes:
            slope_idx = np.where(slopes==s)[0]
            if len(slope_idx)>0:
                idx_median = int(np.median(slope_idx))
                unq_window_median.append(window_center[idx_median])
            else:
                unq_window_median.append(np.nan)
        
        ax1=plt.subplot(212, sharex=ax0)       
        # plot all slopes
        ax1.plot(window_center, slopes, '.', color='LightBlue')
        # plot unique slope values        
        ax1.plot(unq_window_median,unq_slopes,'r.')
        
        plt.suptitle('model '+str(self.n_model))
        ax0.set_ylabel('$\lambda F_\lambda$')
        ax1.set_ylabel('slope')
        ax1.set_xlabel('$\lambda$')
        plt.setp(ax0.get_xticklabels(), visible=False)
        plt.subplots_adjust(hspace=.0)
        ## put the model parameter values below the title?
        plt.show()
        
        ##TO DO:
        # push to github, write for observations!!
        
        
        
        
               
