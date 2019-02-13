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

###TO DO:#
###fix limits for observation x axis?
###make input for static method list instead of just 2?
###but should we only normalize in overplot? For now I like that it defaults to it

# define path to observation files
obs_path = 'seds/reduced_seds/'

# set up plot axes
def setup_plot():
    """
    puts a plot on a log scale and labels the axes for an SED
    """
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('$\lambda [\mu m]$')
    plt.ylabel(r'$\lambda F_{\lambda} [\frac{W}{m^2}]$')

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
    def __init__(self, obs_name, norm=True):
    
        # retrieve data from observation file
        self.obs_name = obs_name
        self.filepath = obs_path+obs_name+'_reduced.sed'
        lines = open(self.filepath).readlines()
        data = np.array([parseline(l) for l in lines[2:]])
        
        # compute wavelength, % error, and SED in W/m^s
        self.wavelength = np.array([float(i) for i in data[:,0]])
        flux_jy = [float(i) for i in data[:,1]]
        sed = np.array([j*1e-23*3e10/(l*1e-4) for j,l in zip(flux_jy,self.wavelength)])
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
        plt.errorbar(self.wavelength, self.sed, yerr=np.abs(self.err), uplims = self.uplims,
        fmt='o', color='darkslateblue', markersize=3, alpha=0.8, label=self.obs_name)
        
        setup_plot()
        plt.title(self.obs_name)
        #plt.show()
        
    @staticmethod
    def overplot(obs1, obs2):
        """
        overplots two observed SED's and creates a corresponding legend
        """
        # plot each observation and label accordingly
        plt.errorbar(obs1.wavelength, obs1.sed, yerr=np.abs(obs1.err), uplims = obs1.uplims,
        fmt='o', markersize = 3,label=obs1.obs_name)
        
        plt.errorbar(obs2.wavelength, obs2.sed, yerr=np.abs(obs2.err), uplims = obs2.uplims,
        fmt='o', markersize = 3,label=obs2.obs_name)
         
        setup_plot()
        plt.legend()
        plt.show()
        
    def get_slopes(self, window):

        # take log for analysis
        logw = np.log10(self.wavelength)
        logs = np.log10(self.sed)
        window = np.log10(window)
    
        # remove points that are upper limits
        logw = logw[self.uplims==False]
        logs = logs[self.uplims==False]
        
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
        
        # finish plot        
        ax1=plt.subplot(212, sharex=ax0)       
        # plot all slopes
        ax1.plot(window_center, slopes, '.', color='LightBlue')
        # plot unique slope values        
        ax1.plot(unq_window_median,unq_slopes,'r.')
        # set labels
        plt.suptitle(self.obs_name)
        ax0.set_ylabel('log($\lambda F_\lambda$)')
        ax1.set_ylabel('slope')
        ax1.set_xlabel('log($\lambda$)')
        ax1.set_ylim(-5,5)
        plt.setp(ax0.get_xticklabels(), visible=False)
        #plt.subplots_adjust(hspace=.0)
        ## put the model parameter values below the title?
        plt.show()