from models import *
from model_grid import grid_parameters

d = fits.open('binary_array.fits')[0].data

def plotP_1D(paramstr):
    """
    plot histogram for parameter showing probability count
    """
    # update dictionary to inclue inclinations
    param_dict = model_grid.grid_parameters.copy()
    param_dict['inc']=[str(i)[:2] for i in list(Model(0).inclinations)]
    
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
    plt.step(a,values,where='mid',color='darkslateblue')
    plt.xticks(a,param_list)
    
    plt.xlabel(paramstr);plt.ylabel('$N_{edge-on}$')
    plt.show()
    
    
def plotP_2D(paramstrx,paramstry):
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
 
    plt.xticks(np.arange(0,len(param_listx)),labelsx)
    plt.yticks(np.arange(0,len(param_listy)),labelsy)
    plt.xlabel(paramstrx);plt.ylabel(paramstry)
    
    plt.imshow(arr, vmin=0, vmax=1, aspect='auto')
    plt.colorbar();plt.show()
    
    # shared colorbar for subplots
    #cax = plt.axes([0.85, 0.1, 0.075, 0.8])
    #plt.colorbar(cax=cax)
    