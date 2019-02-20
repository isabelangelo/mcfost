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
from find_model import *
from observations import Obs

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
            
        # generate model inclinations- is it in the files?
        i_0, i_f = np.radians(45), np.radians(90)
        cosi = np.linspace(np.cos(i_0), np.cos(i_f), 15)
        self.inclinations = np.degrees(np.arccos(cosi))
        
        # store parameters
        self.parameters = get_parameters(n_model)
        
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
                # plot polyfit as a test
                #ax0.plot(x,a*x+b, '-')
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
        
        print(unq_slopes)        
        #return slopes_all, slopes_unq
        
    def plot_slopes(self, templates=False):
    
        # set up figure subplots
        ax0 = plt.subplot(211)
        ax1=plt.subplot(212, sharex=ax0)
        
        # plot each SED and corresponding slopes
        for i in np.arange(0,15,7):
            slopes_all, slopes_unq = self.get_slopes(i)

            # plot sed in top panel
            ax0.plot(self.wavelength, self.seds[i], '-')
            # plot all slope values
            #ax1.plot(10**np.array(slopes_all[0]), slopes_all[1], '.', color='LightBlue')
            # plot unique slope values
            ax1.plot(10**np.array(slopes_unq[0]), slopes_unq[1], '.', 
                    label='i='+str(np.round(self.inclinations[i])))
                    
        if templates==True:
            template1_all, template1_unq = Obs('hh30').get_slopes() 
            template2_all, template2_unq = Obs('i04302').get_slopes()
            template3_all, template3_unq = Obs('SSTTAUJ042021.4+281349_low').get_slopes()
            
            ax1.plot(10**np.array(template1_unq[0]),template1_unq[1], 'k-', linewidth=0.4)
            ax1.plot(10**np.array(template2_unq[0]),template2_unq[1], 'k-', linewidth=0.4)
            ax1.plot(10**np.array(template3_unq[0]),template3_unq[1], 'k-', linewidth=0.4)
        
        # set labels and subplots            
        ax0.set_ylabel('log($\lambda F_\lambda$)')
        ax0.set_xscale('log')
        ax0.set_yscale('log')
        plt.setp(ax0.get_xticklabels(), visible=False)
        ax1.set_xscale('log') 
        ax1.set_ylabel('slope')
        ax1.set_xlabel('log($\lambda$)')
        ax1.legend()
        ax0.set_title(str(self.parameters), fontsize=10)
        plt.suptitle('model '+str(self.n_model))
        plt.subplots_adjust(hspace=.05)
        plt.show()
        


##TO DO: maybe we don't need to return all slopes? 
# every value in array for unique slopes ends in 13927?
# add a for loop to plot templates?
        
        
        
               
