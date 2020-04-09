"""
Plot fitted GSMF and original data
"""
import json
import sys
import numpy as np

from scipy.stats import binned_statistic

from matplotlib import cm
import matplotlib.pyplot as plt

from methods import piecewise_linear

from methods import mass_bins, binned_weighted_quantile
exec(open("./obs_data_sfs.py").read())

import flares
fl = flares.flares(fname='../../flares/data/flares.hdf5')

tags = fl.tags
zeds = [float(tag[5:].replace('p','.')) for tag in tags]


## ---- Overdensity Weights
dat = np.loadtxt(fl.weights, skiprows=1, delimiter=',')
weights = dat[:,8]
index = dat[:,0]

## ---- Plot
ticks = np.linspace(0.05, .95, len(tags))
colors = [ cm.viridis(i) for i in ticks ]

fig, (ax1,ax2) = plt.subplots(2,1,figsize=(5.5,10.5))

plt.subplots_adjust(hspace=0.1)

axes = [ax1,ax2]#,ax4,ax5,ax6]

for c,ax,tag,z in zip(colors[4:],axes,tags[4:],zeds[4:]):
    print(tag)

    with open('samples/sfs_fit_%s.json'%tag) as f:
        p = json.load(f)
    
    x = np.linspace(8,12,int(1e3))

    x0,y0,m1,m2 = p['x0']['median'],p['y0']['median'],p['m1']['median'],p['m2']['median']
    ax.plot(x,piecewise_linear(x-9.7,*[x0,y0,m1,m2]),color=c,lw=4)

    ## Observations ##

    mask = (santini17['z_low'] < z-0.1) & (santini17['z_high'] > z-0.1)
    if np.sum(mask) != 0:
        s17_artist = ax.errorbar(np.log10(santini17['mstar'][mask]), np.log10(santini17['sfr'][mask]), 
            yerr=santini17['sigma_sfr'][mask], fmt='p', label='Santini+17', color='grey')


    mask = ((salmon15['z']-0.5) < z) & ((salmon15['z']+0.5) > z)
    if np.sum(mask) != 0: 
        print('salmon15')
        s15_artist = ax.errorbar(salmon15['logM'][mask], salmon15['logSFR'][mask], 
                yerr=salmon15['sigma_MC'][mask], fmt='s', label='Salmon+15', color='grey')


    ax.text(0.1, 0.8, '$z = %.1f$'%z, transform=ax.transAxes, size=15)
    ax.set_xlim(8,11.5)
    ax.set_ylim(-1,3) 
    ax.grid(alpha=0.5)


ax2.set_xlabel('$\mathrm{log_{10}} \, (M_{\mathrm{*}} \,/\, \mathrm{M_{\odot}})$', size=16)

for ax in axes:#[ax4,ax5,ax6]:
    ax.set_ylabel('$\mathrm{log_{10}}\,(\mathrm{SFR} \,/\, \mathrm{M_{\odot}} \, \mathrm{yr^{-1}})$', size=16)

#for ax in [ax1,ax4]:
#    ax.set_ylabel('$\mathrm{log_{10}}\,(\mathrm{SFR} \,/\, \mathrm{M_{\odot}} \, \mathrm{yr^{-1}})$', size=16)



for ax in [ax2]:#[ax5,ax6,ax2,ax3]:
    ax.set_yticklabels([])

#for ax in [ax1,ax2,ax3]:
#    ax.set_xticklabels([])


ax2.legend(frameon=False, loc=3);


#plt.show()
imgf='images/sfs_obs.png'
print(imgf)
fig.savefig(imgf, dpi=150, bbox_inches='tight')

