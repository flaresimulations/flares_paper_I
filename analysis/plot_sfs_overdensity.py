"""
Plot fitted GSMF and original data
"""
import sys
import numpy as np

from scipy.stats import binned_statistic

from matplotlib import cm
import matplotlib.colors as colors
import matplotlib.pyplot as plt

import fitDF.fitDF as fitDF
import fitDF.models as models
import fitDF.analyse as analyse

import flares
fl = flares.flares(fname='../../flares/data/flares.hdf5')

tags = fl.tags
zeds = [float(tag[5:].replace('p','.')) for tag in tags]

massBinLimits = np.linspace(7.5, 13.5, 13)
massBins = np.logspace(7.75, 13.25, 12)
# print(massBinLimits)
# print(np.log10(massBins))
# from methods import mass_bins
# massBins, massBinLimits = mass_bins() 


## Load FLARES ## 
sfr = fl.load_dataset('SFR_inst_30', arr_type='Subhalo')
mstar = fl.load_dataset('Mstar_30', arr_type='Subhalo')
centrals = fl.load_dataset('Centrals', arr_type='Subhalo')
R = 14./0.6777

## ---- Overdensity Weights
dat = np.loadtxt(fl.weights, skiprows=1, delimiter=',')
log1pdelta = dat[:,5]
weights = dat[:,8]
index = dat[:,0]

## ---- Plot GSMF
fig, ax = plt.subplots(1,1,figsize=(6,6))

cNorm  = colors.Normalize(vmin=log1pdelta.min(), vmax=log1pdelta.max())
scalarMap = cm.ScalarMappable(norm=cNorm, cmap=cm.plasma)

# ticks = np.linspace(0.05, .95, len(tags))
# colors = [ cm.viridis(i) for i in ticks ]


# for tag,z,c in zip(tags,zeds,colors):
tag = fl.tags[5]
z = float(tag[5:].replace('p','.'))
print(tag)

for i,(halo,d) in enumerate(zip(fl.halos,log1pdelta)):

    mstar_temp = mstar[halo][tag]#[centrals[halo][tag]]
    sfr_temp = sfr[halo][tag]#[centrals[halo][tag]]
    
    # plot centrals individually
    # plt.scatter(np.log10(mstar_temp[centrals[halo][tag]]),
    #             np.log10(sfr_temp[centrals[halo][tag]]), s=5,color=c)

    
    mask = (mstar_temp > 0.) & (sfr_temp > 1e-3)
    mstar_temp = mstar_temp[mask]
    sfr_temp = sfr_temp[mask]
    
    # plt.scatter(np.log10(mstar_temp),np.log10(sfr_temp), s=5,color=scalarMap.to_rgba(d))        

    if np.sum(mask) > 0:

        _sfs, dummy, dummy = binned_statistic(np.log10(mstar_temp), sfr_temp, statistic='median', bins=massBinLimits)
    
        # _sfs_16, dummy, dummy = binned_statistic(np.log10(mstar_temp), sfr_temp,
        #                                          statistic=lambda y: np.percentile(y, 16), bins=massBinLimits)

        # _sfs_84, dummy, dummy = binned_statistic(np.log10(mstar_temp), sfr_temp,
        #                                          statistic=lambda y: np.percentile(y, 84), bins=massBinLimits)

        hist, dummy, dummy = binned_statistic(np.log10(mstar_temp), sfr_temp, statistic='count', bins=massBinLimits)

        mask = hist >= 0
        plt.plot(np.log10(massBins[mask]), np.log10(_sfs[mask]), color=scalarMap.to_rgba(d), lw=1)
        # plt.fill_between(np.log10(massBins[mask]), np.log10(_sfs_16[mask]), np.log10(_sfs_84[mask]),
        #                  color=scalarMap.to_rgba(d), alpha=0.2, linewidth=0)



ax.text(0.1, 0.1, '$z = %.1f$'%z, transform=ax.transAxes, size=16)
ax.set_xlim(8,11.5)
ax.set_ylim(-1,3) 
ax.grid(alpha=0.5)


dbinLims = [-0.3,-0.15,-0.075,0.0,0.1,0.2,0.28,0.3] # np.linspace(0.3,-0.3,8)
dbins = dbinLims[:-1] +  np.diff(dbinLims)/2
dbins = np.array(["%.2f"%db for db in dbins]).astype(float)

ticks = np.linspace(0.05,.95,len(dbins))
sm = plt.cm.ScalarMappable(cmap='plasma', norm=plt.Normalize(vmin=0, vmax=1))
sm._A = []  # # fake up the array of the scalar mappable
cbaxes = fig.add_axes([0.72, 0.17, 0.03, 0.28])
cbar = plt.colorbar(sm, ticks=ticks, cax=cbaxes)
cbar.ax.set_yticklabels(dbins)
cbar.ax.set_ylabel('$\mathrm{log_{10}}(1 \,+\,\delta)$', size=16, rotation=90)
 
ax.set_ylabel('$\mathrm{log_{10}}\,(\mathsf{SFR} \,/\, \mathrm{M_{\odot}} \, \mathrm{yr^{-1}})$', size=16)
ax.set_xlabel('$\mathrm{log_{10}} \, (M_{\mathrm{*}} \,/\, \mathrm{M_{\odot}})$', size=16)


#plt.show()
imgf='images/sfs_overdensity.png'
print(imgf)
fig.savefig(imgf, dpi=150, bbox_inches='tight')

