"""
This code contains a variety of methods for plotting outputs form MCFOST
"""
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import glob

def plot_SED(filename, label=None):
    """
    plot SED [W/m^2] versus wavelength [um]
    
    Args:
        filename(str): name of path to sed file, example: 'model1/data_th/sed_rt.fits'
        label(str): label of SED to be put in a legend
    """
    hdulist = fits.open(filename)
    sed = hdulist[0].data[0][0][0]
    _lambda = hdulist[1].data
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('$\lambda$')
    plt.ylabel('$\lambda F_{\lambda}$')
    
    if label:
        plt.plot(_lambda, sed, label=label)
        plt.legend()
    else: 
        plt.plot(_lambda, sed)

    
def plot_2d(filename, inc=0):
    """
    plot image output from MCFOST
    
    Args:
        filename(str): name of path to image file, example: 'model1/data_0.6/RT.fits'
        inc(int): index of model inclination to display
    """
    hdulist = fits.open(filename)
    sed = hdulist[0].data[0][0][inc]
    plt.imshow(sed, vmax=20*sed.mean())
    plt.xlim(40,60)
    plt.ylim(40,60)
    #plt.colorbar()
    
def SED_error(datapath):
    """
    compute error for multiple SEDs generated for given set of parameters 
    
    Args:
        datapath(str): path containing all data_th_* directories
    Returns:
        computed error of all SED's generated
    """
    sedlist = glob.glob(datapath+'/data*/sed_rt.fits')
    seds = []
    for sed in sedlist:
        d = fits.open(sed)[0].data[0][0]
        seds.append(d)
    stds = np.std(seds, axis=0)
    return np.mean(stds)
        
def plot_SED_irange(fileroot):
    """
    plot SEDs [W/m^2] versus wavlenght [um] for all inclinations of model
    
    Args:
        fileroot(str): directory where model is stored, example: 'model1'
    """
    hdulist = fits.open(fileroot+'/data_th_1/sed_rt.fits')
    sed = hdulist[0].data[0][0]
    n_i = len(sed)
    _lambda = hdulist[1].data
    
    for i in range(n_i):
        plt.plot(_lambda, sed[i], linewidth=0.5)
    
    plt.ylim(1e-19,1e-11)    
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('$\lambda$')
    plt.ylabel('$\lambda F_{\lambda}$')
    plt.text(100,1e-12, str(SED_error(fileroot)))
    plt.title(fileroot)
    
def plot_SED_errorbars(fileroot, t):
    """
    plot average SEDs [W/m^2] versus wavelength [um] with errorbars
    for set of SEDs for a given parameter set
    
    Args:
        fileroot(str): path containing all data_th_* directories
        t(str): float representing average runtime of models
    """
    sedpaths = glob.glob(fileroot+'/data*/sed_rt.fits')
    _lambda = fits.open(sedpaths[0])[1].data

    seds = []
    edge_on_seds = []
    
    for path in sedpaths:
        d = fits.open(path)[0].data
        seds.append(d[0][0][:-1])
        edge_on_seds.append(d[0][0][-1])
        
    sed_avg = np.mean(seds, axis=0)
    errors = np.std(edge_on_seds, axis=0)
    edge_on_sed = np.mean(edge_on_seds, axis=0)
    
    plt.errorbar(_lambda, edge_on_sed, yerr=errors, color='k', ecolor='r', elinewidth = 2)
    for sed in sed_avg:
        plt.plot(_lambda, sed, color='grey', linewidth=0.5)
        
    plt.ylim(1e-19,1e-11)    
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('$\lambda$')
    plt.ylabel('$\lambda F_{\lambda}$')
    plt.text(100,1e-12, 't~'+t)
    plt.title(fileroot)
    
    
def plot_lambda(fileroot, n):
    """
    plot a MCFOST wavelength file
    
    Args:
        fileroot(str): name of wavelength file, example: 'Lambda/lambda3'
        n(int): determines where on the y axis the lambda values will lie
    """
    lambdalist = list(open(fileroot+'.lambda'))
    lambdax = [float(x) for x in lambdalist]
    plt.plot(lambdax, n*np.ones(len(lambdax)), '.', label=fileroot)
    plt.xlabel('$\lambda [\mu m]$')
    plt.xscale('log')
    plt.yticks([])
    plt.axvline(x=0.4,color='k')
    plt.axvline(x=0.7,color='k')