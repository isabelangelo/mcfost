from model_grid import get_grid_indexes
from models import *
from test_edge_on import *
import matplotlib.gridspec as gridspec

# controlling for Rc
Rc10 = get_grid_indexes(Rc=10, n=6)
Rc30 = get_grid_indexes(Rc=30, n=6)
Rc100 = get_grid_indexes(Rc=100, n=6)
Rc300 = get_grid_indexes(Rc=300, n=6)

# high dust mass
high_dust_mass = get_grid_indexes(dust_mass=1e-3, n=6)

# control for flaring exponent
f_exp085 = get_grid_indexes(f_exp=0.85, n=6)
f_exp1 = get_grid_indexes(f_exp=1.0, n=6)
f_exp115 = get_grid_indexes(f_exp=1.15, n=6)
f_exp13 = get_grid_indexes(f_exp=1.3, n=6)

# random sample from grid
random1 = get_grid_indexes(n=6)
random2 = get_grid_indexes(n=6)

def makeplot(gridlist, title):
    fig=plt.figure(figsize=(6, 8))
    gs = gridspec.GridSpec(4,4)
    ax = fig.add_subplot(gs[0, :])
    
    # for images
    subplot_arr = np.arange(0,12,1)[0::2]
    num = 0
    for i in subplot_arr:
        m = Model(gridlist[num])

        p_threshold = 0.5
        p = image_compute_P(m)
        ax.plot(m.inclinations, p, label=str(num))
        ax.set_ylim(0,1)
    
        ax2 = fig.add_subplot(gs[i+4])
        ax3 = fig.add_subplot(gs[i+5])
        if len(np.where(p>p_threshold)[0])==0:
            image2 = m.convolved_images[-2]
            image3 = m.convolved_images[-1]
            ax3.annotate('never edge on',(0,100))
        else:
            first_edgeon_idx = np.where(p>p_threshold)[0][0]
            image2 = m.convolved_images[first_edgeon_idx-1]
            image3 = m.convolved_images[first_edgeon_idx]
            ax2.annotate(str(m.inclinations[first_edgeon_idx-1])[:2],(100,100))
            ax3.annotate(str(m.inclinations[first_edgeon_idx])[:2],(100,100))
              
        vmin2 = 10*image2.mean();vmax2 = image2.max() 
        vmin3 = 10*image3.mean();vmax3 = image3.max() 
                                                                                                           
        ax2.imshow(image2, norm=LogNorm(),\
                    vmin=vmin2,vmax=vmax2, cmap='RdPu')
        ax3.imshow(image3, norm=LogNorm(),\
                    vmin=vmin3,vmax=vmax3, cmap='RdPu')
                    
        ax2.yaxis.set_major_locator(plt.NullLocator())
        ax2.xaxis.set_major_formatter(plt.NullFormatter())
        ax3.yaxis.set_major_locator(plt.NullLocator())
        ax3.xaxis.set_major_formatter(plt.NullFormatter())
        ax2.set_ylabel(str(num), rotation=0)

        num += 1

    # for SEDS   
    #for i in range(len(gridlist)):
    #    m = Model(gridlist[i])
        #p = compute_P(m)
        #ax.plot(m.inclinations, p, label=str(i))

        #for j in range(len(m.seds)):                                                                                                       
        #    ax2.loglog(m.wavelength, m.seds[j], color='k', linewidth=0.5)
        #ax2.set_xticks([])
        #ax2.set_yticks([])
        #ax2.annotate(str(i),(1e3,1e-14))

    ax.set_xlabel('$i$')
    ax.set_ylabel('P_avg')
    ax.set_title(title)
    ax.legend()
    plt.savefig('../'+title+'.png')
    plt.close()
    print(title)

makeplot(Rc10, 'Rc = 10')
makeplot(Rc30, 'Rc = 30')
makeplot(Rc100, 'Rc = 100')
makeplot(Rc300, 'Rc = 300')
makeplot(high_dust_mass, 'high dust mass (1e-3)')
makeplot(f_exp085, 'flaring exponent = 0.85')
makeplot(f_exp1, 'flaring exponent = 1.0')
makeplot(f_exp115, 'flaring exponent = 1.15')
makeplot(f_exp13, 'flaring exponent = 1.3')
makeplot(random1, 'randomly selected models 1')
makeplot(random2, 'randomly selected models 2')



