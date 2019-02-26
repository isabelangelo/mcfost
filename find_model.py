"""
The purpose of this code is to retrieve models from an MCFOST grid with a set
of desired model parameters

Written: Isabel Angelo (2018)
"""

from ast import literal_eval
from generate_para import dust_mass, Rc, f_exp, H0, Rin, sd_exp, porosity, amax
import sys

# retrieve list with indices corresponding to parameter combinations
f = open('model_indices_15i.txt','r')
model_list = [literal_eval(line.strip()) for line in f]

# store grid parameters in dictionary
grid_keys = ['dust_mass', 'Rc','f_exp', 'H0', 'Rin', 'sd_exp', 'porosity', 'amax']
grid_values = [dust_mass, Rc, f_exp, H0, Rin, sd_exp, porosity, amax]
grid_parameters = dict(zip(grid_keys, grid_values))

# find index corresponding to set of input parameters
def get_model(parameter_input):
    """
    reads in list of parameter values and outputs the name of the grid model with 
    corresponding parameters
    
    Args:
        parameter_input(str): string containing comma separated values of parameters in the format
        'dust_mass,Rc,f_exp,H0,Rin,sd_exp,porosity,amax'
        example: '1e-07,10,0.85,5,0.1,-0.5,0,10'
        *note: for our models, porosity will always be 0
    
    Returns:
        model_index(int): index of model_list array corresponding to model with input parameters,
        for example model_index=4 corresponds to model4.para
    """

    parameter_list = [float(p) for p in parameter_input.split(',')]
    model_index = model_list.index(tuple(parameter_list))
    return model_index
    
def get_parameters(model_index):
    """
    takes in model number and returns dictionary of corresponding model parameters
    
    Args:
        model_index(int): model number corresponding to index in model_list
        
    Returns:
        parameter_dict(dict): dictionary with parameter names as keys with their
        corresponding MCFOST values
    """
    
    parameters = model_list[model_index]
    keys = ['dust_mass', 'Rc','f_exp', 'H0', 'Rin', 'sd_exp', 'porosity', 'amax']
    values = list(parameters)
    parameter_dict = dict(zip(keys,values))
    return parameter_dict
