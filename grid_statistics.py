from models import *

d = fits.open('binary_array.fits')[0].data
n_models = d.size/15.

# plot dust mass versus inclination color map
dust_mass_values = np.sum(d,axis=(6,5,4,3,2,1))
dust_mass_max = n_models/len(model_grid.grid_parameters['dust_mass'])
dust_mass_arr = dust_mass_values/dust_mass_max

inc = [str(i)[:2] for i in Model(0).inclinations]
dm = [str(d) for d in model_grid.grid_parameters['dust_mass']]

fig, ax = plt.subplots()
ax.set_xticks(np.arange(0,14,1));ax.set_xticklabels(inc);ax.set_yticklabels(dm)
ax.set_xlabel('inclination');ax.set_ylabel('dust mass')
ax.set_title('dust mass edge-on probabilities')
plt.imshow(dust_mass_arr)
plt.colorbar();plt.show()

#bug: why does this code take so long
#also some ticks dont show up