import numpy as np
from astropy.io import fits
from generate_para import dust_mass, Rc, f_exp, H0, Rin, sd_exp, porosity, amax
from model_grid import get_model_index
from models import *
from test_edge_on import *

arr = np.empty((5,4,4,4,3,3,4,15)) # put in lengths to be more straightforward?

for a in range(len(dust_mass)):
    a_s = str(dust_mass[a])
    for b in range(len(Rc)):
        b_s = str(Rc[b])
        for c in range(len(f_exp)):
            c_s = str(f_exp[c])
            for d in range(len(H0)):
                d_s = str(H0[d])
                for e in range(len(Rin)):
                    e_s = str(Rin[e])
                    for f in range(len(sd_exp)):
                        f_s = str(sd_exp[f])
                        for g in range(len(amax)):
                            g_s = str(amax[g])
                            
                            # get index with corresponding model                            
                            param_str = ','.join([a_s,b_s,c_s,d_s,e_s,f_s,'0',g_s])
                            idx = get_model_index(param_str)
                            p = compute_P(Model(idx))
                            arr[a,b,c,d,e,f,g,:] = p
                            print(idx)
                
# write the large array to a text file
hdu = fits.PrimaryHDU(data=arr)
hdu.writeto('probability_array.fits')