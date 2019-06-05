"""
This code computes statistics on model grid from MCFOST 
and defines for plotting them.

written: Isabel Angelo (2019)
"""

###NOTE : this code can't really go between both grids yet,
# but that's okay because you can't make the arrays for the new one yet anyways

from models import *
from model_grid import grid_parameters
from probability_array import generate_masked_array

# load binary array with probabilities
binary_arr_path = model_path + 'binary_array.fits'
binary_array = fits.open(binary_arr_path)[0].data

# generate array with masked values
# this gets changed when you want to mask different values
masked_array = generate_masked_array(
    dust_mass = [1e-07, 1e-06, 1e-05, 0.0001, 0.001],
    Rc = [10, 30, 100, 300],
    f_exp = [0.85, 1.0, 1.15, 1.3],
    H0 = [5, 10, 15, 20],
    Rin =  [0.1, 1, 10],
    sd_exp = [0, -0.5, -1, -1.5],
    amax = [10, 100, 1000, 10000]
    )

# generate weighted array here

# update array with masks and weights    
binary_array *= weighted_array # might not need masked array??


# generate array to multiply by binary that maps

# update dictionary to inclue inclinations
param_dict = model_grid.grid_parameters.copy()
param_dict['inc']=[str(i)[:2] for i in list(Model(0).inclinations)]

def plotP_1D(ax, paramstr):
    """
    plot histogram for parameter showing probability count
    """
    # get list of parameter values
    param_list = param_dict[paramstr]
    
    # collapse axes
    keys = list(param_dict.keys())
    axes = [7,6,5,4,3,2,1,0]
    axes.remove(keys.index(paramstr))
    values = np.sum(binary_array,axis=tuple(axes)) # do we want to normalize?
    max = binary_array.size/(len(param_list))
    arr = values/max
    arr2=np.concatenate(([arr[0]],arr,[arr[-1]])) # looks better in plot
    
    # plot values in steps
    a = np.arange(0,len(param_list))
    a2 = np.arange(0,len(param_list)+2)
    ax.step(a2,arr2,where='mid',color='darkslateblue')
    ax.fill_between(a2,arr2,step='mid',alpha=0.6,color='darkslateblue')
    
    ax.set_xticks(a+1)
    ax.set_xticklabels(param_list,rotation=-45)
    #ax.tick_params(direction='in')
    
    ax.set_xlim(a2[0]+0.5,a2[-1]-0.5)
    ax.set_ylim(0,1)
    
    ax.set_xlabel(paramstr);ax.set_ylabel(paramstr)
    #plt.show()
    
    
def plotP_2D(fig, ax, paramstrx,paramstry,cmap='Purples',cbar=False):
    """
    plot correlation image for 2 input parameters
    """
    # get list of parameter values
    param_listx = param_dict[paramstrx]
    param_listy = param_dict[paramstry]
    
    # generate axes to collapse
    keys = list(param_dict.keys())
    axes = [7,6,5,4,3,2,1,0]
    axes.remove(keys.index(paramstrx))
    axes.remove(keys.index(paramstry))
    
    # collapse into 2D array and normalize
    values = np.sum(binary_array,axis=tuple(axes))
    max = binary_array.size/(len(param_listx)*len(param_listy))
    arr = values/max
    
    # determine tick labels
    labelsx = [str(i) for i in param_listx]
    labelsy = [str(i) for i in param_listy]
    
    # flip array so that labels will match values
    if keys.index(paramstrx)<keys.index(paramstry):
        arr = arr.T
        
    ax.set_xticks(np.arange(0,len(param_listx)));ax.set_xticklabels(labelsx,rotation=-45)
    ax.set_yticks(np.arange(0,len(param_listy)));ax.set_yticklabels(labelsy)
    
    ax.set_xlabel(paramstrx);ax.set_ylabel(paramstry)
    im = ax.imshow(arr, origin='lower',vmin=0.1, vmax=1, norm=LogNorm(),aspect='auto',cmap=cmap)
    #plt.colorbar()#;plt.show()
    
    # plot colorbar
    if cbar==True:
        cb = fig.colorbar(im, cax=plt.axes([0.1,0,1,0.02]), orientation='horizontal')
    
def plot_corner():
    """
    generates corner plot where all histogram and 2D correlation plots are shown
    """
    # list of parameter names and corresponding indexes
    pnames = ['inc']+list(grid_parameters.keys())
    ilist = np.arange(0,len(pnames),1)
    
    # generate combinations for correlation plots
    pcomb = [(pnames[a],pnames[b]) for a in ilist for b in ilist if b>=a]
    icomb = [(a,b) for a in ilist for b in ilist if b>=a]
    allcomb = [(a,b) for a in ilist for b in ilist]
    
    # create figure
    fig, axes = plt.subplots(8,8, figsize=(6,6))
    for comb in allcomb:
        ax = axes[comb[1],comb[0]]
        if comb in icomb:
            x,y=pnames[comb[1]],pnames[comb[0]]
            if x==y:
                plotP_1D(ax, x)
            elif comb==(0,1):
                plotP_2D(fig,ax,y,x,cbar=True)
            else:
                plotP_2D(fig,ax,y,x,cmap='Purples')
        else:
            ax.axis('off')
        # remove unnecessary labels/ticks
        if comb[1]!=7:
            ax.set_xlabel('')
            #ax.set_xticks([])
        if comb[0]!=0:
            ax.set_ylabel('')
            ax.set_yticks([])
        if comb[0]==0 and comb[1]==7:
            ax.set_xticks([2,7,12])
            ax.set_xticklabels(['52','70','84'])
            
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.show()
    