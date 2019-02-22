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
import find_model

# define path to models
#model_path = 'grid_i15/' # add ../ for leo
model_path='/Volumes/backup/grid_i15/'

# set up plot axes
def setup_plot(ax):
    """
    puts a plot on a log scale and labels the axes for an SED
    """
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('$\lambda [\mu m]$')
    ax.set_ylabel(r'$\lambda F_{\lambda} [\frac{W}{m^2}]$')

class Model(object):
    """
    Model object
    
    Args:
        n_model(int): integer associated with corresponding MCFOST model
        norm(bool): 'True' will normalize according to the geometric mean 
                        between 2-10um of the highest inclination in the model, 
                        'False' plots raw data
    """    
    def __init__(self, n_model, norm=False):
        
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
            
        # generate model inclinations- is it in the files?
        i_0, i_f = np.radians(45), np.radians(90)
        cosi = np.linspace(np.cos(i_0), np.cos(i_f), 15)
        self.inclinations = np.degrees(np.arccos(cosi))
        
        # store parameters
        self.parameters = find_model.get_parameters(n_model)
        
    def plot(self):
        """
        plots an MCFOST model SED [W/m^2] versus wavelength [um]
        """        
        # plot SEDs for all model inclinations
        f, ax = plt.subplots()
        n_i = len(self.seds)
        for i in range(n_i):                                                                                                       
            ax.plot(self.wavelength, self.seds[i], color='k', linewidth=0.5)
            
        setup_plot(ax)
        ax.set_xlim(self.wavelength[0],self.wavelength[-1])
        ax.set_title('model '+str(self.n_model))
        plt.show()
        
    @staticmethod
    def overplot(model1, model2):
        """
        overplots two MCFOST model SED's with multiple inclinations and creates a 
        corresponding legend
        """
        # plot SEDs for all model inclinations for each model
        f, ax = plt.subplots()
        n_i = len(model1.seds)
        for i in range(n_i):                                                                                                       
            lines1 = ax.plot(model1.wavelength, model1.seds[i], color='k', linewidth=0.5)
            lines2 = ax.plot(model2.wavelength, model2.seds[i], color='r', linewidth=0.5)
    
        # set up legend that only shows both models without individual inclinations
        ax.legend([lines1[-1],lines2[-1]], 
            ['model '+str(model1.n_model), 'model '+str(model2.n_model)])
            
        setup_plot(ax)   
        plt.show() 
            
    def get_slopes(self, inc_idx, window=4):

        # take log for analysis
        logw = np.log10(self.wavelength)
        logs = np.log10(self.seds[inc_idx])
        window = np.log10(window)
        
        # define bins
        window_start = np.arange(logw[0],logw[-1],0.01)
        window_center = window_start+window/2.
        
        # store slopes for each window
        slopes = []      
        for i in window_start:
            fitpoints = np.where((i<=logw)&(logw<=i+window))[0]
            if len(fitpoints)>=2:
                # perform linear fit
                x = logw[fitpoints]
                y = logs[fitpoints]
                a,b=np.polyfit(x,y,1)
                # store slope
                slopes.append(a)
            else:
                slopes.append(np.nan)
      
        # store unique slope and corresponding window values (median index)
        indexes = np.unique(slopes, return_index=True)[1]
        unq_slopes = np.array([slopes[index] for index in sorted(indexes)])
        unq_slopes = unq_slopes[~np.isnan(unq_slopes)]
        unq_window_median = []
        
        for s in unq_slopes:
            slope_idx = np.where(slopes==s)[0]
            idx_median = int(np.median(slope_idx))
            unq_window_median.append(window_center[idx_median])
        
        # store window+slope data as tuples       
        slopes_all = (window_center, slopes)
        slopes_unq = (unq_window_median, unq_slopes)
               
        return slopes_unq
        
    def plot_slopes(self):
    
        # set up figure subplots
        ax0 = plt.subplot(211)
        ax1=plt.subplot(212, sharex=ax0)
        
        # plot each SED and corresponding slopes
        for i in np.arange(0,15,7):
            slopes_unq = self.get_slopes(i)
            # plot sed in top panel
            ax0.plot(self.wavelength, self.seds[i], '-')
            # plot unique slope values
            ax1.plot(10**np.array(slopes_unq[0]), slopes_unq[1], '.', 
                    label='i='+str(np.round(self.inclinations[i])))
        
        # set labels and subplots            
        plt.suptitle('model '+str(self.n_model))
        ax0.set_title(str(self.parameters), fontsize=10)
        plt.setp(ax0.get_xticklabels(), visible=False)
        plt.subplots_adjust(hspace=.05)
        setup_plot(ax0)
        ax1.set_xscale('log') 
        ax1.set_ylabel('slope')
        ax1.set_xlabel('log($\lambda$)')
        ax1.legend()
        plt.show()
        
    def brightness_test1(self, inc_idx=-1):
        
        # set wavelength for brightness calculation
        star_peak = 0.75 # found by taking peak brightness of star with 1e-12 dust mass
        
        # compute brightness ratio at 45 and 90 degrees
        i_45, i_90 = self.seds[0],self.seds[inc_idx]
        F_45 = i_45[np.where(self.wavelength==star_peak)][0]
        F_90 = i_90[np.where(self.wavelength==star_peak)][0]
        
        def P(x):
            a,b,c,d = -0.332,53.1,15.5,0.5
            y = a*np.arctan(b*x-c)+d
            if y<0:
                return 0
            else:
                return y
                
        return P(F_90/F_45)
        
        
    def brightness_test2(self, inc_idx=-1):
        pass
        
        
        
               
