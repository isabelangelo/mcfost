"""
This code defines the Model class for handling model SED's from MCFOST.
The class contains the SEDs, their wavelengths, and methods
for plotting individual and multiple SEDs.

Written: Isabel Angelo (2019)
"""
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm 
import scipy.stats.mstats as stats
from astropy.convolution import convolve_fft
from reduce_sed import parseline
from model_grid import *
from probability_functions import *

import imageio
import glob
import os

# define path to HST PSF
tinytim_PSF = fits.open('../result00_psf.fits')[0].data
tinytim_PSFnorm = tinytim_PSF/np.max(tinytim_PSF) # normalized PSF
PSF_20mas = rebin_image(tinytim_PSFnorm, 0.0158,0.02) # rebinned PSF to 20mas/pixel

# store model star SED
star_path = model_path + 'test_peak/data_th/sed_rt.fits'
star_seds = fits.open(star_path)[0].data[0][0]

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
            
        # store the model image
        impath = self.filepath[:-19]+'data_0.6/RT.fits'
        self.images = fits.open(impath)[0].data[0][0]
        
        # store convolved+rebinned images image
        # convolve images with rebinned PSF (20mas/pixel)
        convolved_prebinned_images = [convolve_fft(i, PSF_20mas) for i in self.images]
        # rebin to 40mas/pixel via averaging
        binpix = int((self.images[0].shape[0]-1)/2)
        self.convolved_images = [rebin(i[:-1,:-1], (binpix,binpix)) for i in convolved_prebinned_images]
            
        # generate model inclinations- is it in the files?
        i_0, i_f = np.radians(45), np.radians(90)
        cosi = np.linspace(np.cos(i_0), np.cos(i_f), 15)
        self.inclinations = np.degrees(np.arccos(cosi))
        
        # store parameters
        self.parameters = get_model_parameters(n_model)
        
    def plot(self):
        """
        plots an MCFOST model set of SEDs [W/m^2] versus wavelength [um]
        """        
        # plot SEDs for all model inclinations
        f, ax = plt.subplots()
        n_i = len(self.seds)
        for i in range(n_i):                                                                                                       
            ax.plot(self.wavelength, self.seds[i], color='k', linewidth=0.5)
        
        # plot host star
        ax.plot(self.wavelength, star_seds[0],color='midnightblue', linewidth=0.8, linestyle=':')
            
        setup_plot(ax)
        ax.set_xlim(self.wavelength[0],self.wavelength[-1])
        ax.set_title(str(self.parameters), fontsize=10)
        plt.suptitle('model '+str(self.n_model))
        plt.show()
       
    def create_gif(self):
        """
        plot a gif that shows both the SED and output image change
        as inclination increases.
        """
    
        # generate plots as pngs
        for inc_idx in range(0,15):
            fig,ax =plt.subplots(1,3,figsize=(17,5))
            
            # plot SED inclinations
            n_i = len(self.seds)
            for i in range(n_i): 
                if i==inc_idx: # plot input inclination in bold
                    ax[0].plot(self.wavelength, self.seds[i], 'k-') 
                else:                                                                                                   
                    ax[0].plot(self.wavelength, self.seds[i], color='gray', linewidth=0.5)
            setup_plot(ax[0])
            
            # plot host star
            ax[0].plot(self.wavelength, star_seds[0],color='midnightblue', linewidth=0.8, linestyle=':')
            # print inclination
            ax[0].annotate('i='+str(np.round(self.inclinations[inc_idx]))[:-2]+'$^\circ$',\
                            (0.87,0.95),xycoords='axes fraction')
    
            # get image and convulation data
            image = self.images[inc_idx] # image
            convolved_image = self.convolved_images[inc_idx]
    
            vmin = 10*image.mean();vmax = image.max() 
            cm = 'RdPu'
    
            # plot original image
            im1 = ax[1].imshow(image, norm=LogNorm(),\
                    vmin=vmin,vmax=vmax, aspect='auto', cmap=cm)
            #fig.colorbar(im1, ax=ax[1])#orientation='horizontal',ax=ax[1])
            #ax[1].set_title('MCFOST Image')
                
            # plot convolved with tinytim PSF
            im2 = ax[2].imshow(convolved_image, norm=LogNorm(),\
                    vmin=vmin,vmax=vmax, aspect='auto', cmap=cm)
            fig.colorbar(im2, ax=ax[2],fraction=0.05,pad=0.01)#orientation='horizontal',ax=ax[2])
            #ax[2].set_title('Convolved with HST PSF')
    
            fig.suptitle('Model '+str(self.n_model))
            fig.text(.5, .9, str(self.parameters)[1:-1], ha='center', fontsize=9)
            plt.savefig('../'+str(inc_idx)+'.png')
            plt.close()
        
        # combine pngs into .gif
        images=[]
        for i in range(0,15):
            filename='../'+str(i)+'.png'
            images.append(imageio.imread(filename))
            os.remove(filename)
        imageio.mimsave('../Model'+str(self.n_model)+'.gif', images)
         
        
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
        """
        computes the slope if the SED as a function of wavelength
    
        Args:
            inc_index(int): index of inclination that slopes are calculated for, 
            corresponding to index of self.inclinations
            window(int): size of the window for which the slope is calculated.
            window=4 corresponds to window spanning a factor of 4, i.e. 0.1-0.4
        
        Returns:
            slopes_unq(tuple): tuple (unq_window, unq_slopes) where unq_slopes
            represents each specific slope value computed for the SED and 
            unq_window represents their median corresponding index
        """    
        # take log for analysis
        logw = np.log10(self.wavelength)
        logs = np.log10(self.seds[inc_idx])
        window = np.log10(window)
        
        # define bins
        window_start = np.arange(logw[0],logw[-1],0.01)
        window_center = window_start+window/2.
        
        # store slopes for each window
        slopes = np.array([])      
        for i in window_start:
            fitpoints = np.where((i<=logw)&(logw<=i+window))[0]
            if len(fitpoints)>=2:
                # perform linear fit
                x = logw[fitpoints]
                y = logs[fitpoints]
                a,b=np.polyfit(x,y,1)
                # store slope
                slopes = np.append(slopes, a)
            else:
                slopes = np.append(slopes, np.nan)
      
        # store unique slope and corresponding window values (median index)
        indexes = np.unique(slopes, return_index=True)[1]
        unq_slopes = np.array([slopes[index] for index in sorted(indexes)])
        unq_slopes = unq_slopes[~np.isnan(unq_slopes)]
        unq_window = np.array([])
        
        for s in unq_slopes:
            slope_idx = np.where(slopes==s)[0]
            idx_median = int(np.median(slope_idx))
            unq_window = np.append(unq_window, window_center[idx_median])
        
        # store window+slope data as tuples       
        slopes_all = (window_center, slopes)
        slopes_unq = (unq_window, unq_slopes)
               
        return slopes_unq
        
    def plot_slopes(self):
        """
        plots the SED's of a model in top panel and their slopes as a function
        of wavelength in the bottom panel
        """
        # set up figure subplots
        ax0 = plt.subplot(211)
        ax1=plt.subplot(212, sharex=ax0)
        
        # plot each SED and corresponding slopes
        for i in np.arange(0,15,7):
            slopes_unq = self.get_slopes(i)
            # plot sed in top panel
            ax0.plot(self.wavelength, self.seds[i], '-')
            # plot unique slope values
            ax1.plot(10**slopes_unq[0], slopes_unq[1], '.', 
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
         
        
          
        
               
