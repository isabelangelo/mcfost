"""
This code defines Obs, a class for handling observed SEDs. The class contains the SEDs, 
their wavelengths, and methods for plotting individual and multiple SEDs. 

Written: Isabel Angelo (2019)
"""
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
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
        
        # fix upper limit arrow lengths to 10% of upper limit value
        self.err[self.uplims]=-0.3*self.sed[self.uplims]
    
    def plot(self):
        """
        plots an observed SED [W/m^2] versus wavelength [um]
        """
        # plot SED with errorbars and upper limits
        plt.errorbar(self.wavelength, self.sed, yerr=np.abs(self.err), uplims = self.uplims,
        fmt='o', color='darkslateblue', markersize=3, alpha=0.8, label=self.obs_name)
        
        setup_plot()
        plt.title(self.obs_name)
        
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