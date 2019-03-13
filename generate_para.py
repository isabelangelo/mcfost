"""
This code will define and generate a grid of parameter files to be run on MCFOST. 
Lists contain the variable values of the grid, and the parameter files are 
stored as model*.para.

Written: Isabel Angelo (2018)
"""
# define grid
dust_mass = [1e-7, 1e-6, 1e-5, 1e-4, 1e-3]

Rc = [10, 30, 100, 300]

f_exp = [0.85, 1.0, 1.15, 1.30]

H0 = [5, 10, 15, 20]

Rin = [0.1, 1, 10]

#sd_exp = [0, -0.5, -1.5]
sd_exp = [0, -0.5, -1, -1.5]


amax = [10, 100, 1000, 10000]

# generate all possible combinations of grid parameters
combinations = [(a,b,c,d,e,f,g) for a in dust_mass for b in Rc for c in f_exp \
for d in H0 for e in Rin for f in sd_exp for g in amax]

# generate initial parameters
from pymcfost import parameters
para = parameters.McfostParams(filename='ref3.0.para')

# update fixed parameters
# para.phot.nphot_image = 1.28e4
# para.simu.use_default_wl = False
# #para.wavelengths.file = 'lambda3.lambda'
# para.wavelengths.file = 'lambda4.lambda'
# #para.map.nx = 251
# #para.map.ny = 251
# para.map.nx = 375
# para.map.ny = 375
# para.map.size=1050 #added for second grid
# para.map.RT_imin = 45
# para.map.RT_imax = 90
# para.map.RT_ntheta = 15
# para.simu.radial_migration = False
# para.simu.dust_sublimation = False
# para.simu.hydrostatic_eq = False
# para.simu.viscous_heating = False
# para.zones[0].Rout = 600 #added for second grid
# para.zones[0].geometry = 2
# para.zones[0].dust[0].amin = 0.01
# para.zones[0].dust[0].n_grains = 50
# para.zones[0].dust[0].amax = 10
# 
# # generate individual parameter files
# for i in range(len(combinations)):
#     new_para = para
#     comb = combinations[i]
#     new_para.zones[0].dust_mass = comb[0]
#     new_para.zones[0].Rc = comb[1]
#     new_para.zones[0].flaring_exp = comb[2]
#     new_para.zones[0].h0 = comb[3]
#     new_para.zones[0].Rin = comb[4]
#     new_para.zones[0].surface_density_exp = comb[5]
#     #new_para.zones[0].dust[0].porosity = comb[6]
#     new_para.zones[0].dust[0].amax = comb[6]
#     new_para.writeto('model'+str(i)+'.para')



