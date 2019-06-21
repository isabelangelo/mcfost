"""
The purpose of this code is to retrieve models from an MCFOST grid with a set
of desired model parameters

Written: Isabel Angelo (2018)
"""

from ast import literal_eval
from generate_para import dust_mass, Rc, f_exp, H0, Rin, sd_exp, amax
import numpy as np
import sys

# define path to models
model_path='/Volumes/backup/disks/grid_new/' # new grid
#model_path='/Volumes/backup/disks/old_grid/' # old grid

# store which grid is being used
if 'new' in model_path:
    grid = 'new'
else:
    grid = 'old'

# retrieve list with indices corresponding to parameter combinations
model_idx_path = model_path+ 'model_indices.txt'
f = open(model_idx_path,'r')
model_list = [literal_eval(line.strip()) for line in f]

# store grid parameters in dictionary
grid_keys = ['dust_mass', 'Rc','f_exp', 'H0', 'Rin', 'sd_exp', 'amax']
grid_values = [dust_mass, Rc, f_exp, H0, Rin, sd_exp, amax]
grid_parameters = dict(zip(grid_keys, grid_values))


def get_model_index(parameter_input):
    """
    reads in list of parameter values and outputs the name of the grid model with 
    corresponding parameters
    
    Args:
        parameter_input(str): string containing comma separated values of parameters in the format
        'dust_mass,Rc,f_exp,H0,Rin,sd_exp,amax'
        example: '1e-07,10,0.85,5,0.1,-0.5,10'
        *note: for our models, porosity will always be 0
    
    Returns:
        model_index(int): index of model_list array corresponding to model with input parameters,
        for example model_index=4 corresponds to model4.para
    """

    parameter_list = [float(p) for p in parameter_input.split(',')]
    
    if grid=='old':
        parameter_list.insert(-1,0) # for porosity=0 in grid
        
    model_index = model_list.index(tuple(parameter_list))
    return model_index
    
def get_model_parameters(model_index):
    """
    takes in model number and returns dictionary of corresponding model parameters
    
    Args:
        model_index(int): model number corresponding to index in model_list
        
    Returns:
        parameter_dict(dict): dictionary with parameter names as keys with their
        corresponding MCFOST values
    """
    
    parameters = model_list[model_index]
    keys = ['dust_mass', 'Rc','f_exp', 'H0', 'Rin', 'sd_exp', 'amax']
    values = list(parameters)
    
    if grid == 'old':
        values.pop(-2) # remove porosity
        
    parameter_dict = dict(zip(keys,values))
    return parameter_dict
    
def get_grid_indexes(dust_mass=None, Rc=None, f_exp=None, H0=None,\
    Rin=None, sd_exp=None, amax=None,  n=None):
    """
    returns all model numbers of grid 
    satisfying parameter conditions specified by inputs
    
    Args:
        dust_mass, Rc, f_exp, etc, (float,int): desired values
        of corresponding parameters to be retrieved from grid
        
        n(int): number of randomly selected model indices to generate. 
        Returns all possible indices if None.
        
        example: get_grid_indexes(dust_mass=1e-3) 
        returns all model #'s with dustmass=1e-3
    """
    # generate lists of indexes for each individual condisiton
    indexes_all = []
    if dust_mass is not None:
        a = [i for i,x in enumerate(model_list) if x[0]==dust_mass]
        indexes_all.append(a)
    if Rc is not None:
        b = [i for i,x in enumerate(model_list) if x[1]==Rc]
        indexes_all.append(b)
    if f_exp is not None:
        c = [i for i,x in enumerate(model_list) if x[2]==f_exp]
        indexes_all.append(c)
    if H0 is not None:
        d = [i for i,x in enumerate(model_list) if x[3]==H0]
        indexes_all.append(d)
    if Rin is not None:
        e = [i for i,x in enumerate(model_list) if x[4]==Rin]
        indexes_all.append(e) 
    if sd_exp is not None:
        f = [i for i,x in enumerate(model_list) if x[5]==sd_exp]
        indexes_all.append(f)
    if amax is not None:
        g = [i for i,x in enumerate(model_list) if x[6]==amax]
        indexes_all.append(g)
    
    # find grid indexes that satisfy all conditions
    if len(indexes_all)>0: 
        grid_indexes = set(indexes_all[0])
        for s in indexes_all[1:]:
            grid_indexes.intersection_update(s)
        grid_indexes=list(grid_indexes)
        
    # for case where no conditions are input    
    else:
        grid_indexes = list(np.arange(0,len(model_list),1))
    
    # generate random model numbers from grid indexes   
    if n is not None:
        return np.random.choice(grid_indexes, n)       
    else:
        return grid_indexes

