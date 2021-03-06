"""
This code computes statistics on model grid from MCFOST 
and defines for plotting them.

written: Isabel Angelo (2019)
"""

###NOTE : this code can't really go between both grids yet,
# but that's okay because you can't make the arrays for the new one yet anyways

from models import *
from model_grid import *
from probability_array import *
from weights import *

# load binary array with probabilities
binary_arr_path = model_path + 'SED_and_image.fits'
binary_array_unflipped = fits.open(binary_arr_path)[0].data

# flip binary array along flaring exp axis
binary_array = np.flip(binary_array_unflipped, axis=5)

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

# generate weights array here
unweighted_values = [[1,1,1,1,1], # dust_mass_weights
                [1,1,1,1], # Rc_weights
                [1,1,1,1], # f_exp_weights
                [1,1,1,1], # H0_weights
                [1,1,1], # Rin_weights
                [1,1,1,1], # sd_exp_weights
                [1,1,1,1]] # amax weights


weighted_values = [mass_taurus_truncated, #[1,1,1,1,1], # dust_mass_weights
                [1,1,1,0], # Rc_weights
                [0,1,1,1], # f_exp_weights
                H0_gaussian_norm, # [1,1,1,1] # H0_weights
                [1,1,1], # Rin_weights
                [1,1,1,1], # sd_exp_weights
                [1,1,1,1]] # amax weights


# CHANGE THIS LINE TO CHOOSE WHICH WEIGHT SET
grid_weights = dict(zip(grid_keys, unweighted_values))
    

weighted_array = generate_weighted_array(
    dust_mass_weights = grid_weights['dust_mass'],
    Rc_weights = grid_weights['Rc'],
    f_exp_weights = grid_weights['f_exp'],
    H0_weights = grid_weights['H0'],
    Rin_weights =  grid_weights['Rin'],
    sd_exp_weights = grid_weights['sd_exp'],
    amax_weights = grid_weights['amax']
    )


# generate array to multiply by binary that maps

# update dictionary to inclue inclinations
param_dict = grid_parameters.copy()
param_dict['inc']=[str(i)[:2] for i in list(Model(0).inclinations)]

# dictionary of parameter names + label strings
pnames = ['inc']+list(grid_parameters.keys())
pvalues = ['$i$ [$^o$]', '$M_d$ [$M_\odot$]', '$R_c$ [AU]', r'$\beta$', '$H_0$ [AU]', '$R_{in}$ [AU]', r'$\gamma$', '$a_{max}$ [$\mu$m]']
pnames_dict=dict(zip(pnames, pvalues))

#dictionary with axis labels
labels_dict = {'dust_mass':['$10^{-7}$','$10^{-6}$','$10^{-5}$','$10^{-4}$','$10^{-3}$'], 
              'Rc':[str(i) for i in param_dict['Rc']],
              'f_exp':[str(i) for i in param_dict['f_exp']],
               'H0':[str(i) for i in param_dict['H0']],
              'Rin':[str(i) for i in param_dict['Rin']],
              'sd_exp': [str(i) for i in param_dict['sd_exp']][::-1],
              'amax': ['$10$','$10^{2}$','$10^{3}$','$10^{4}$'],
              'inc': [str(i) for i in param_dict['inc']]}

def plotP_1D(ax, paramstr,linecolor='black'):
    """
    plot histogram for parameter showing probability count
    """
    # get list of parameter values
    param_list = param_dict[paramstr]
    
    # collapse axes
    keys = list(param_dict.keys())
    axes = [7,6,5,4,3,2,1,0]
    axes.remove(keys.index(paramstr))
    
    weighted_binary_array = binary_array*weighted_array #bool*weight
    numerator = np.sum(weighted_binary_array,axis=tuple(axes))#take sum
    denominator = np.sum(weighted_array, axis=tuple(axes))#sum weights
    arr = numerator/denominator

    arr2=np.concatenate(([arr[0]],arr,[arr[-1]])) # looks better in plot
    
    # plot values in steps
    a = np.arange(0,len(param_list))
    a2 = np.arange(0,len(param_list)+2)
    ax.step(a2,arr2,where='mid',color=linecolor)
#    ax.fill_between(a2,arr2,step='mid',alpha=0.6,color='darkslateblue')
    
    ax.set_xticks(a+1)
    ax.set_xticklabels(labels_dict[paramstr],rotation=40)
    ax.xaxis.set_tick_params(labelsize=7,rotation=40)
    
    ax.set_yticks([0,0.5,1])
    ax.set_yticklabels(labels=[0,0.5,1],rotation=40)
    ax.yaxis.set_tick_params(labelsize=7)
    #ax.tick_params(direction='in')
    
    ax.set_xlim(a2[0]+0.5,a2[-1]-0.5)
    ax.set_ylim(0,1)
    
    ax.set_xlabel(pnames_dict[paramstr]);ax.set_ylabel(pnames_dict[paramstr])
    #plt.show()
    
    
def plotP_2D(fig, ax, paramstrx,paramstry,cmap='Purples',cbar=False):
    """
    plot correlation image for 2 input parameters
    
    weighted array is computed as follows
    
    weighted_value = sum(bool*weight)/sum(weight)
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
    weighted_binary_array = binary_array*weighted_array #bool*weight
    numerator = np.sum(weighted_binary_array,axis=tuple(axes))#take sum
    denominator = np.sum(weighted_array, axis=tuple(axes))#sum weights
    arr = numerator/denominator
    
    # determine tick labels
    labelsx = labels_dict[paramstrx]
    labelsy = labels_dict[paramstry]
    
    # flip array so that labels will match values
    if keys.index(paramstrx)<keys.index(paramstry):
        arr = arr.T
        
    ax.set_xticks(np.arange(0,len(param_listx)));ax.set_xticklabels(labelsx,rotation=40)
    ax.xaxis.set_tick_params(labelsize=7)
    
    ax.set_yticks(np.arange(0,len(param_listy)));ax.set_yticklabels(labelsy,rotation=40)
    ax.yaxis.set_tick_params(labelsize=7)
    
    ax.set_xlabel(pnames_dict[paramstrx]);ax.set_ylabel(pnames_dict[paramstry])
    im = ax.imshow(arr, origin='lower',vmin=0.1, vmax=1, norm=LogNorm(),aspect='auto',cmap=cmap)
    #plt.colorbar()#;plt.show()
    
    # plot colorbar
    if cbar==True:
        cb = fig.colorbar(im, cax=plt.axes([0.1,0,1,0.02]), orientation='horizontal')
        
def plot_corner(cm='Purples'):
    """
    generates corner plot where all histogram and 2D correlation plots are shown
    """
    # list of parameter names and corresponding indexes
    ilist = np.arange(0,len(pnames),1)
    
    # generate combinations for correlation plots
    pcomb = [(pnames[a],pnames[b]) for a in ilist for b in ilist if b>=a]
    icomb = [(a,b) for a in ilist for b in ilist if b>=a]
    allcomb = [(a,b) for a in ilist for b in ilist]
    
    # create figure
    fig, axes = plt.subplots(8,8, figsize=(10,10))
    for comb in allcomb:
        ax = axes[comb[1],comb[0]]
        if comb in icomb:
            x,y=pnames[comb[1]],pnames[comb[0]]
            if x==y:
                plotP_1D(ax, x)
            elif comb==(0,1):
                plotP_2D(fig,ax,y,x,cmap=cm)#,cbar=True)
            else:
                plotP_2D(fig,ax,y,x,cmap=cm)
        else:
            ax.axis('off')
        # remove unnecessary labels/ticks
        if comb[1]!=7: #remove x labels from inner subplots
            ax.set_xlabel('')
            #ax.set_xticks([])
        if comb[0]!=0: # remove y label from inner subplots
            ax.set_ylabel('')
            ax.set_yticks([])
        if comb[0]==0 and comb[1]==7: # label inclination x axis
            ax.set_xticks([2,7,12])
            ax.set_xticklabels(['52','70','84'])
        if (comb[0]==0 and comb[1]==1): # label dust mass y axis 
            ax.set_yticks([0,2,4])
        if (comb[0]==1 and comb[1]==7): # label dust mass x axis 
            ax.set_xticks([0,2,4])
            ax.set_xticklabels(['$10^{-7}$','$10^{-5}$','$10^{-3}$'])
            
    plt.subplots_adjust(wspace=0, hspace=0)
    #plt.savefig('corner_plots/'+s+'_corner.png')
    plt.show()
    
def bar_graph(ax, paramstr):
    """
    plot histogram for parameter showing probability count
    """
    # get list of parameter values
    param_list = param_dict[paramstr]
    values = [float(i) for i in param_dict[paramstr]]
    
    # collapse axes
    keys = list(param_dict.keys())
    axes = [7,6,5,4,3,2,1,0]
    axes.remove(keys.index(paramstr))
    
    weighted_binary_array = binary_array*weighted_array #bool*weight
    numerator = np.sum(weighted_binary_array,axis=tuple(axes))#take sum
    denominator = np.sum(weighted_array, axis=tuple(axes))#sum weights
    arr = numerator/denominator

    plt.bar(values, arr, width=2, label='SED', color='purple', alpha=0.4)
    #plt.show()
    