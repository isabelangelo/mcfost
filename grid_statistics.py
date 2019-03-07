from models import *
from model_grid import grid_parameters

d = fits.open('binary_array.fits')[0].data
n_models = d.size/15.

def plot_parameter_probabilities(paramstr, cmap='inferno'):

    # get list of parameter values
    param_list = model_grid.grid_parameters[paramstr]
    
    # generate axes to collapse
    keys = list(grid_parameters.keys()) # you can import this instead
    axes = [6,5,4,3,2,1,0]
    axes.remove(keys.index(paramstr))

    # collapse into 2D array and normalize
    values = np.sum(d,axis=tuple(axes))
    max = n_models/len(param_list)
    arr = values/max

    # determine tick labels
    inc = [str(i)[:2] for i in Model(0).inclinations]
    dm = [str(d) for d in param_list]
    
    plt.imshow(arr, vmin=0, vmax=1, cmap=cmap, aspect='auto')
    plt.xticks(np.arange(0,15,1),inc);plt.yticks(np.arange(0,len(param_list)),dm)
    plt.ylabel(paramstr);plt.xlabel('inclination')
    plt.title(paramstr+' edge-on probabilities')
    plt.colorbar();plt.show()

    # shared colorbar for subplots
    #cax = plt.axes([0.85, 0.1, 0.075, 0.8])
    #plt.colorbar(cax=cax)
    
    
#TO DO: make this so you can correlate any input two?
#maybe the default can be inclination?