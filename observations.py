"""
This code defines Obs, a class for handling observed SEDs. The class contains the SEDs, 
their wavelengths, and methods for plotting individual and multiple SEDs. 

Written: Isabel Angelo (2019)
"""
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import warnings
import scipy.stats.mstats as stats
from reduce_sed import parseline
from probability_functions import *

# define path to observation files
obs_path = '/Volumes/backup/disks/reduced_seds/'
#obs_path = 'seds/reduced_seds/'

# set up plot axes
def setup_plot(ax):
    """
    puts a plot on a log scale and labels the axes for an SED
    """
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('$\lambda [\mu m]$')
    ax.set_ylabel(r'$\lambda F_{\lambda} [\frac{W}{m^2}]$')
    
# dm/dy term for error calculations
def dmdy(j, x_arr, N):
    num = np.sum(x_arr)-N*x_arr[j]
    denom = (np.sum(x_arr))**2.-N*np.sum(x_arr**2.)
    return num/denom

# suppress polyfit warnings    
warnings.simplefilter('ignore', np.RankWarning)

class Obs(object):
    """
    Observation object
    
    Args:
        obs_name(str): string name of object that appears as the fileroot in the original
                        SED text file. For example saucer_reduced.sed will have
                        obs_name 'saucer'
        norm(bool): 'True' will normalize according to the geometric mean 
                    between 2-10um of observed SED, 'False' plots raw data
    """
    def __init__(self, obs_name, norm=False):
    
        # retrieve data from observation file
        self.obs_name = obs_name
        self.filepath = obs_path+obs_name+'_reduced.sed'
        lines = open(self.filepath).readlines()
        data = np.array([parseline(l) for l in lines[2:]])
        
        # compute wavelength, % error, and SED in W/m^s
        self.wavelength = np.array([float(i) for i in data[:,0]])
        flux_jy = [float(i) for i in data[:,1]]
        sed = np.array([j*1e-26*3e8/(l*1e-6) for j,l in zip(flux_jy,self.wavelength)])
        percent_err = [float(i) for i in data[:,2]]
        
        # normalize according observations
        if norm==True:
            mean_range = np.where((2<self.wavelength) & (self.wavelength<10))
            gmean = stats.gmean(sed[mean_range])
            self.sed = sed/gmean
        
        # for non-normalized case
        else:
            self.sed = sed
        
        # convert errors 
        self.err = np.array([i*j*0.01 for i,j in zip(self.sed,percent_err)])
        
        # determine points that are upper limits
        indices = [i for i, s in enumerate(self.err) if '-'==str(s)[0]]
        self.uplims = self.err<0
        self.uplims[indices]=True
        
        # fix upper limit arrow lengths to 30% of upper limit value
        self.err[self.uplims]=-0.3*self.sed[self.uplims]
        
        # set 0 error to default 10% of flux
        zero_err = np.where(self.err==0)
        self.err[zero_err]=0.1*self.sed[zero_err]
    
    def plot(self):
        """
        plots an observed SED [W/m^2] versus wavelength [um]
        """
        # plot SED with errorbars and upper limits
        f, ax = plt.subplots()
        ax.errorbar(self.wavelength, self.sed, yerr=np.abs(self.err), uplims = self.uplims,
        fmt='o', color='darkslateblue', markersize=3, alpha=0.8, label=self.obs_name)
        
        setup_plot(ax)
        ax.set_title(self.obs_name)
        plt.show()
        
    @staticmethod
    def overplot(obs1, obs2):
        """
        overplots two observed SED's and creates a corresponding legend.
        example: Obs.overplot(Obs('hh30'), Obs('hvtauc))
        
        Args: 
            obs1, obs2(Obs): observation objects to be plotted
        Returns:
            overplotted observations with legend
        """
        # plot each observation and label accordingly
        f, ax = plt.subplots()
        ax.errorbar(obs1.wavelength, obs1.sed, yerr=np.abs(obs1.err), uplims = obs1.uplims,
        fmt='o', markersize = 3,label=obs1.obs_name)
        
        ax.errorbar(obs2.wavelength, obs2.sed, yerr=np.abs(obs2.err), uplims = obs2.uplims,
        fmt='o', markersize = 3,label=obs2.obs_name)
         
        setup_plot(ax)
        ax.legend()
        plt.show()
        
    def get_slopes(self, window=4):
        """
        computes the slope if the SED as a function of wavelength
    
        Args:
            window(int): size of the window for which the slope is calculated.
            window=4 corresponds to window spanning a factor of 4, i.e. 0.1-0.4
        
        Returns:
            slopes_unq(tuple): tuple (unq_window, unq_slopes, unq_err)
            where unq_slopes represents each specific slope value computed for the SED 
            and unq_window, unq_err represents their median corresponding index and errors
        """  

        # take log for analysis
        logw = np.log10(self.wavelength)
        logs = np.log10(self.sed)
        logerr = 0.434*self.err/self.sed
        window = np.log10(window)
    
        # remove points that are upper limits
        logw = logw[self.uplims==False]
        logs = logs[self.uplims==False]
        
        # define bins
        window_start = np.arange(logw[0],logw[-1],0.01)
        window_center = window_start+window/2.
        
        # store slopes for each window
        slopes, slope_err = np.array([]), np.array([])     
        for i in window_start:
            fitpoints = np.where((i<=logw)&(logw<=i+window))[0]
            if len(fitpoints)>=2:
                # perform linear fit
                x = logw[fitpoints]
                y = logs[fitpoints]
                a,b=np.polyfit(x,y,1)
                # store slope
                slopes = np.append(slopes,a)
                # store slope error
                dmdy_arr = [dmdy(j,x,len(x))for j in range(len(x))]
                sigma_arr = [logerr[i] for i in fitpoints] #convert to log?
                dm_arr = [t1**2.*t2**2. for t1,t2 in zip(dmdy_arr,sigma_arr)]
                dm = np.sqrt(np.sum(dm_arr))
                slope_err = np.append(slope_err, dm)
                # plot polyfit as a test
                #ax0.plot(x,a*x+b, '-')
            else:
                slopes = np.append(slopes, np.nan)
                slope_err = np.append(slope_err, np.nan)
              
        # store unique slope and corresponding window values (median index)
        indexes = np.unique(slopes, return_index=True)[1]
        unq_slopes = np.array([slopes[index] for index in sorted(indexes)])
        unq_slopes = unq_slopes[~np.isnan(unq_slopes)]
        unq_window, unq_err= np.array([]), np.array([])
        
        for s in unq_slopes:
            slope_idx = np.where(slopes==s)[0]
            idx_median = int(np.median(slope_idx))
            unq_window = np.append(unq_window, window_center[idx_median])
            unq_err = np.append(unq_err, slope_err[idx_median])
        
        # store window+slope data as tuples       
        slopes_all = (window_center, slopes, slope_err)
        slopes_unq = (unq_window, unq_slopes, unq_err)
                
        return slopes_unq
        
    def plot_slopes(self):
        """
        plots the observed SED in top panel and its slope as a function
        of wavelength in the bottom panel
        """
    
        # set up figure subplots
        ax0 = plt.subplot(211)
        ax1=plt.subplot(212, sharex=ax0)            
        # plot sed in top panel
        slopes_unq = self.get_slopes()
        ax0.errorbar(self.wavelength, self.sed, yerr=self.err, 
                    fmt = '.--', label=self.obs_name)
        # plot unique slope values        
        ax1.errorbar(10**slopes_unq[0], slopes_unq[1], yerr=slopes_unq[2],
                     fmt='.', markersize=3)
            
        # set labels and subplots
        #plt.suptitle(self.obs_name)
        plt.setp(ax0.get_xticklabels(), visible=False)
        ax0.legend()
        plt.subplots_adjust(hspace=.05)
        setup_plot(ax0)
        ax1.set_xscale('log')
        #ax1.set_ylim(-5,5)
        ax1.set_ylabel('slope')
        ax1.set_xlabel('log($\lambda$)')
        plt.subplots_adjust(hspace=None)
        #plt.show()
        
    
    
    
    
    
    
        