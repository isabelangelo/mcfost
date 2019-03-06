from model_grid import get_grid_indexes
from models import *
from test_edge_on import *
import matplotlib.gridspec as gridspec

# extreme cases
high_dust_mass = get_grid_indexes(dust_mass=1e-3, n=6)
steep_surface_density = get_grid_indexes(sd_exp=-1.5,n=6)
small_Rin = get_grid_indexes(Rin=0.1,n=6)
small_flaring_exp = get_grid_indexes(f_exp=0.85,n=6)

# control for each extreme parameter
mid_dust_mass = get_grid_indexes(dust_mass=1e-5, n=6)
mid_surface_density = get_grid_indexes(sd_exp=-0.5,n=6)
mid_Rin = get_grid_indexes(Rin=1,n=6)
mid_flaring_exp = get_grid_indexes(f_exp=1,n=6)

# random sample from grid
random = get_grid_indexes(n=6)

def makeplot(gridlist, title):
    fig=plt.figure(figsize=(5, 6))
    gs = gridspec.GridSpec(4,2)

    ax = fig.add_subplot(gs[0, :])
    for i in range(len(gridlist)):
        m = Model(gridlist[i])

        p = compute_P(m)
        ax.plot(m.inclinations, p, label=str(i))
    
        ax2 = fig.add_subplot(gs[i+2])
        for j in range(len(m.seds)):                                                                                                       
            ax2.loglog(m.wavelength, m.seds[j], color='k', linewidth=0.5)
        ax2.set_xticks([])
        ax2.set_yticks([])
        ax2.annotate(str(i),(1e3,1e-14))


    ax.set_xlabel('$i$')
    ax.set_ylabel('P_avg')
    ax.set_title(title)
    ax.legend()
    plt.savefig(title+'.png')
    plt.close()
    print(title)

makeplot(random, 'randomly selected')
makeplot(high_dust_mass, 'high dust mass (1e-3)')
makeplot(steep_surface_density, 'steep surface density')
makeplot(small_Rin, 'small inner radius')
makeplot(small_flaring_exp, 'small flaring exponent')
makeplot(mid_dust_mass, 'medium dust mass (1e-5)')
makeplot(mid_surface_density, 'medium surface density')
makeplot(mid_Rin, 'medium inner radius')
makeplot(mid_flaring_exp, 'medium flaring exponent')



