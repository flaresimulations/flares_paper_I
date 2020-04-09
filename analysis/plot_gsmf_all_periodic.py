"""
Plot fitted GSMF and original data
"""
import sys
import numpy as np

import matplotlib as mpl
from matplotlib import cm
import matplotlib.pyplot as plt
# mpl.rcParams['text.usetex'] = True

import fitDF.fitDF as fitDF
import fitDF.models as models
import fitDF.analyse as analyse

model = models.DoubleSchechter()

# massBinLimits = np.linspace(7.45, 13.25, 30)
# massBins = np.logspace(7.55, 13.15, 29)
from methods import mass_bins
massBins, massBinLimits = mass_bins() 
# print(np.log10(massBins))
# print(massBinLimits)

import flares

## ---- Plot GSMF
fig, ax = plt.subplots(1,1,figsize=(6,6))
    



for sim,bsize,fname in zip(['agn','ref'],[50,100],['EAGLE_AGNdT9_sp_info','EAGLE_REF_sp_info']):
    print(sim,bsize)
    if sim == 'agn': continue

    fl = flares.flares(fname='../../flares/data/%s.hdf5'%fname, sim_type='PERIODIC')

    tags = fl.tags
    zeds = [float(tag[5:].replace('p','.')) for tag in tags]
    ticks = np.linspace(0.05, .95, len(tags))
    colors = [ cm.viridis(i) for i in ticks ]
    
    ## Load FLARES ## 
    import h5py
    mstar = {}
    with h5py.File('../../flares/data/%s.hdf5'%fname, 'r') as f:
        for tag in fl.tags:
            mstar[tag] = f['%s/Galaxy/Mstar_30'%tag][:]
    #mstar = fl.load_dataset('Mstar_30',arr_type='Subhalo')
    
    for tag,z,c in zip(tags,zeds,colors):
        
        phi, phi_sigma, hist = fl.calc_df(mstar[tag], tag, bsize**3, massBinLimits)
    
        ## ---- Get fit
        #sample_ID = '%s_gsmf_%s'%(sim,tag)
        #a = analyse.analyse(ID='samples', model=model, sample_save_ID=sample_ID, verbose=False)
        
        # from methods import switch_samples
        # _samples = switch_samples(a.samples)
        #a = analyse.analyse(ID='samples', model=model, sample_save_ID=sample_ID, verbose=False, samples=_samples)
    
        fl.plot_df(ax, phi, phi_sigma, hist, massBins=massBins, color=c, lines=False, label='')
        
        # model.update_params(a.median_fit)
        # ax.plot(np.log10(massBins), a.model.log10phi(np.log10(massBins)), color=c)



## FLARES ##
fl_zoom = flares.flares(fname='../../flares/data/flares.hdf5')
for tag,c in zip(fl_zoom.tags,colors):
    print(tag)
    sample_ID = 'flares_gsmf_%s'%(tag)
    a = analyse.analyse(ID='samples', model=model, sample_save_ID=sample_ID, verbose=False)
    model.update_params(a.median_fit)
    xvals = np.linspace(7,15,1000)
    ax.plot(xvals, a.model.log10phi(xvals), color=c, alpha=1, ls='dashed')



# ax.text(0.1, 0.1, '$z = %.1f$'%z, transform=ax.transAxes, size=12)
ax.set_xlim(8.0, 11.6)
ax.set_ylim(-8.5, -1.5) 
ax.grid(alpha=0.5)

# sm = plt.cm.ScalarMappable(cmap='viridis', norm=plt.Normalize(vmin=0, vmax=1))
# sm._A = []  # # fake up the array of the scalar mappable
# cbaxes = fig.add_axes([0.74, 0.46, 0.03, 0.4]) 
# cbar = plt.colorbar(sm, ticks=ticks, cax=cbaxes)
# cbar.ax.set_yticklabels(zeds)
# cbar.ax.set_ylabel('$z$', size=18, rotation=90)
    
ax.set_ylabel('$\mathrm{log_{10}}\,(\phi \,/\, \mathrm{Mpc^{-3}} \, \mathrm{dex^{-1}})$', size=16)
ax.set_xlabel('$\mathrm{log_{10}} \, (M_{*} \,/\, \mathrm{M_{\odot}})$', size=16)

# ax.legend(frameon=False, loc=1);

# plt.show()
imgf='images/gsmf_all_%s.png'%sim
print(imgf)
#plt.show()
fig.savefig(imgf, dpi=150, bbox_inches='tight')

