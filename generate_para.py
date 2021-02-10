from pymcfost import parameters

# define grid
host_mass = [0.15, 0.3, 0.6, 1.2]
dust_fraction = [0.001, 0.01, 0.1]
Rc = [10, 30, 100, 300]
f_exp = [1.05, 1.15, 1.25]
H0 = [5, 10, 15, 20, 25]
Rin = [0.1, 1, 10]
sd_exp = [-0.5, -1]
dust_settling = [1e-5,1e-4,1e-3,'none']

# store star parameters for each host mass
Teff_values = [3080, 3313, 3574, 4048]
R_values = [0.96, 1.30, 1.79, 2.36]
file_values = ['lte3100-4.5.NextGen.fits.gz', 'lte3300-4.5.NextGen.fits.gz', \
               'lte3600-4.5.NextGen.fits.gz', 'lte4000-4.5.NextGen.fits.gz']

# generate all possible combinations of grid parameters
combinations = [(a,b,c,d,e,f,g,h) for a in host_mass \
                for b in dust_fraction \
                for c in Rc \
                for d in f_exp \
                for e in H0 \
                for f in Rin \
                for g in sd_exp \
                for h in dust_settling]

# generate initial parameters
para = parameters.Params(filename='ref3.0.para')

# update fixed parameters 
para.phot.nphot_image = 1.28e4
para.simu.use_default_wl = False
para.wavelengths.file = 'lambda4.lambda'
para.map.nx = 375
para.map.ny = 375
para.map.size=1050
para.map.RT_imin = 45
para.map.RT_imax = 90
para.map.RT_ntheta = 15
para.simu.radial_migration = False
para.simu.dust_sublimation = False
para.simu.hydrostatic_eq = False
para.simu.viscous_heating = False
para.zones[0].Rout = 600
para.zones[0].geometry = 2
para.zones[0].dust[0].amin = 0.01
para.zones[0].dust[0].n_grains = 50
para.stars[0].is_bb = False
para.simu.dust_settling_type = 3 # Fromang settling
# amax already fixed @1mm

# generate individual parameter files
for i in range(len(combinations)):
    new_para = para
    comb = combinations[i]
    # set host mass and associated parameters
    M = comb[0]; i_M=host_mass.index(M)
    para.stars[0].Teff = Teff_values[i_M]
    para.stars[0].R = R_values[i_M]
    para.stars[0].file = file_values[i_M]
    # set dust mass according to disk mass fraction
    M_d = comb[1]*M/100 # assuming 1:100 gas-to-dust ratio
    # other parameters
    new_para.zones[0].dust_mass = M_d    
    new_para.zones[0].Rc = comb[2]
    new_para.zones[0].flaring_exp = comb[3]
    new_para.zones[0].h0 = comb[4]
    new_para.zones[0].Rin = comb[5]
    new_para.zones[0].surface_density_exp = comb[6]
    # set settling
    if comb[7]=='none':
        para.simu.dust_settling_type = 0
    else:
        para.simu.viscosity = comb[7]
    # write to new parameter file
    new_para.writeto('ParameterFiles/model'+str(i)+'.para')
