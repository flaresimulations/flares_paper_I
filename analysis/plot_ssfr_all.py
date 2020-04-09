"""
Plot fitted GSMF and original data
"""
import sys
import numpy as np

from scipy.stats import binned_statistic

from matplotlib import cm
import matplotlib.pyplot as plt

import fitDF.fitDF as fitDF
import fitDF.models as models
import fitDF.analyse as analyse

from methods import mass_bins, binned_weighted_quantile

import flares
fl = flares.flares(fname='../../flares/data/flares.hdf5')

tags = fl.tags
zeds = [float(tag[5:].replace('p','.')) for tag in tags]

# massBinLimits = np.linspace(7.5, 13.5, 16)
# massBins = np.logspace(7.7, 13.3, 15)
massBinLimits = np.linspace(7.5, 13.5, 13)
massBins = np.logspace(7.75, 13.25, 12)
print(massBinLimits)
print(np.log10(massBins))

# massBins, massBinLimits = mass_bins() 


## Load FLARES ## 
sfr = fl.load_dataset('SFR_inst_30', arr_type='Subhalo')
mstar = fl.load_dataset('Mstar_30', arr_type='Subhalo')
centrals = fl.load_dataset('Centrals', arr_type='Subhalo')
R = 14./0.6777

## ---- Overdensity Weights
dat = np.loadtxt(fl.weights, skiprows=1, delimiter=',')
weights = dat[:,8]
index = dat[:,0]

## ---- Plot GSMF
fig, ax = plt.subplots(1,1,figsize=(6,6))

ticks = np.linspace(0.05, .95, len(tags))
colors = [ cm.viridis(i) for i in ticks ]



for tag,z,c in zip(tags,zeds,colors):
    print(tag)
    
    mstar_all = []
    sfr_all = []
    w_all = []

    for i,halo in enumerate(fl.halos):
        w = weights[np.where(["%02d"%i == halo for i in index])[0]]
    
        mstar_temp = mstar[halo][tag]#[centrals[halo][tag]]
        sfr_temp = sfr[halo][tag]#[centrals[halo][tag]]

        # plot centrals individually
        plt.scatter(np.log10(mstar_temp[centrals[halo][tag]]),
                    np.log10(sfr_temp[centrals[halo][tag]] / mstar_temp[centrals[halo][tag]]),
                    s=5,color=c)

        mask = (mstar_temp > 0.) & (sfr_temp > 1e-3)

        w_all = np.append(w_all,np.repeat(w,np.sum(mask)))
        mstar_all = np.append(mstar_all,mstar_temp[mask])
        sfr_all = np.append(sfr_all,sfr_temp[mask])


    _sfs = binned_weighted_quantile(np.log10(mstar_all), sfr_all / mstar_all, w_all, bins=massBinLimits, quantiles=[0.16,0.5,0.84])

    plt.plot(np.log10(massBins), np.log10(_sfs[:,1]), color=c, lw=3)
    plt.fill_between(np.log10(massBins), np.log10(_sfs[:,0]), np.log10(_sfs[:,2]), 
                     color=c, alpha=0.4, linewidth=0)



ax.set_xlim(8,11.5)
ax.set_ylim(-9.5,-7.1) 
ax.grid(alpha=0.5)

# sm = plt.cm.ScalarMappable(cmap='viridis', norm=plt.Normalize(vmin=0, vmax=1))
# sm._A = []  # # fake up the array of the scalar mappable
# cbaxes = fig.add_axes([0.73, 0.62, 0.03, 0.24]) 
# cbar = plt.colorbar(sm, ticks=ticks, cax=cbaxes)
# cbar.ax.set_yticklabels(zeds)
# cbar.ax.set_ylabel('$z$', size=16, rotation=90)
   
ax.set_ylabel('$\mathrm{log_{10}}\,(\mathsf{sSFR} \,/\, \mathrm{yr^{-1}})$', size=15)
ax.set_xlabel('$\mathrm{log_{10}} \, (M_{\mathrm{*}} \,/\, \mathrm{M_{\odot}})$', size=15)


# plt.show()
imgf='images/ssfr_all.png'
print(imgf)
fig.savefig(imgf, dpi=150, bbox_inches='tight')

