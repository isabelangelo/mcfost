from models import *
from model_grid import grid_parameters

d = fits.open('binary_array.fits')[0].data

def plotP_1D(ax, paramstr, ticks=True):
    """
    plot histogram for parameter showing probability count
    """
    # update dictionary to inclue inclinations
    param_dict = model_grid.grid_parameters.copy()
    param_dict['inc']=[str(i)[:2] for i in list(Model(0).inclinations)]
    param_dict['sd_exp']=[0, -0.5, -1.5] #REMOVE THIS LATER
    # get list of parameter values
    param_list = param_dict[paramstr]
    
    # collapse axes
    keys = list(param_dict.keys())
    axes = [7,6,5,4,3,2,1,0]
    axes.remove(keys.index(paramstr))
    values = np.sum(d,axis=tuple(axes)) # do we want to normalize?
    
    ## spacing based on distance
    #plt.step(param_list,values,where='mid')
    # forced even spacing
    a = np.arange(0,len(param_list))
    ax.step(a,values,where='mid',color='darkslateblue')
    
    if ticks==False:
        ax.set_xticks([]);ax.set_yticks([])
    else:
        ax.set_xticks(a)
        ax.set_xticklabels(param_list)
        
    ax.set_ylim(0,values.max()+values.min())
    
    #ax.set_xlabel(paramstr);ax.set_ylabel(paramstr)
    # set y limits to 0 so you can see true variance
    #plt.show()
    
    
def plotP_2D(fig, ax, paramstrx,paramstry,cbar=False, ticks=True):
    """
    plot correlation image for 2 input parameters
    """
    # update dictionary to inclue inclinations
    param_dict = model_grid.grid_parameters.copy()
    param_dict['inc']=[str(i)[:2] for i in list(Model(0).inclinations)]

    # get list of parameter values
    param_listx = param_dict[paramstrx]
    param_listy = param_dict[paramstry]
    
    # generate axes to collapse
    keys = list(param_dict.keys())
    axes = [7,6,5,4,3,2,1,0]
    axes.remove(keys.index(paramstrx))
    axes.remove(keys.index(paramstry))
    
    # collapse into 2D array and normalize
    values = np.sum(d,axis=tuple(axes))
    max = d.size/(len(param_listx)*len(param_listy))
    arr = values/max
    
    # determine tick labels
    labelsx = [str(i) for i in param_listx]
    labelsy = [str(i) for i in param_listy]
    
    # flip array so that labels will match values
    if keys.index(paramstrx)<keys.index(paramstry):
        arr = arr.T
        
    if ticks==False:
        ax.set_xticks([]);ax.set_yticks([])
    else:
        ax.set_xticks(np.arange(0,len(param_listx)));ax.set_xticklabels(labelsx)
        ax.set_yticks(np.arange(0,len(param_listy)));ax.set_yticklabels(labelsy)
    
    ax.set_xlabel(paramstrx);ax.set_ylabel(paramstry)
    im = ax.imshow(arr, vmin=0.1, vmax=1, norm=LogNorm(),aspect='auto')
    #plt.colorbar()#;plt.show()
    
    # plot colorbar
    #if cbar==True:
    #    cb = fig.colorbar(im, cax=plt.axes([0.1,0,1,0.02]), orientation='horizontal')
    
def plot_corner():
    ilist = np.arange(0,7,1)
    allcomb = [(a,b) for a in ilist for b in ilist]
    icomb = [(a,b) for a in ilist for b in ilist if b>=a]
    
    pnames = list(grid_parameters.keys())
    pcomb = [(pnames[a],pnames[b]) for a in ilist for b in ilist if b>=a]
    
    fig, axes = plt.subplots(7,7, figsize=(6,6))
    for comb in allcomb:
        ax = axes[comb[1],comb[0]]
        if comb in icomb:
            x,y=pnames[comb[1]],pnames[comb[0]]
            if x==y:
                plotP_1D(ax, x)
            #elif comb==(0,1):
            #    plotP_2D(fig,ax,y,x,cbar=True)
            else:
                plotP_2D(fig,ax,y,x)
        else:
            ax.axis('off')
    plt.show()

    
    
#REMOVE LINE 13
    