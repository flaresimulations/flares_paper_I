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

import flares
fl = flares.flares(fname='../../flares/data/flares.hdf5')

tags = fl.tags
zeds = [float(tag[5:].replace('p','.')) for tag in tags]

model = models.DoubleSchechter()

binLimits = np.linspace(-3.0, 3.0, 21)
bins = np.logspace(-2.85, 2.85, 20)
# binLimits = np.linspace(-0.2, 3.2, 18)
# bins = np.logspace(-0.1, 3.1, 17)

print(np.log10(bins))
print(binLimits)

## Load FLARES ## 
sfr = fl.load_dataset('SFR_inst_30',arr_type='Subhalo')
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
    ## ---- get data
    phi_all = np.zeros(len(bins))
    phi_sigma = np.zeros((len(fl.halos),len(bins)))
    hist_all = np.zeros(len(bins))
    
    for i,halo in enumerate(fl.halos):
    
        w = weights[np.where(["%02d"%i == halo for i in index])[0]]
    
        sfr_temp = sfr[halo][tag]
    
        V = (4./3) * np.pi * (R)**3
   
        hist, dummy = np.histogram(np.log10(sfr_temp), bins = binLimits)
        hist = np.float64(hist)
        phi = (hist / V) / (binLimits[1] - binLimits[0])

        phi_sigma[i] = (np.sqrt(hist) / V) / (binLimits[1] - binLimits[0])

        phi_all += np.array(phi) * w
        phi_sigma[i] *= w 
        hist_all += hist
    
    
    phi_sigma = np.sqrt(np.sum(np.square(phi_sigma), axis=0))
    fl.plot_df(ax, phi_all, phi_sigma, hist_all, massBins=bins, color=c, lines=False, label='', lw=5)

    ## ---- Get fit
    sample_ID = 'flares_sfrf_%s'%(tag)

    a = analyse.analyse(ID='samples', model=model, sample_save_ID=sample_ID, verbose=False)
    model.update_params(a.median_fit)
    
    xvals = np.linspace(-0.3,3,1000)
    ax.plot(xvals, a.model.log10phi(xvals), color=c)


# ax.text(0.1, 0.1, '$z = %.1f$'%z, transform=ax.transAxes, size=12)
ax.set_xlim(-2.5,3.2)
ax.set_ylim(-10, -1.2) 
ax.grid(alpha=0.5)

sm = plt.cm.ScalarMappable(cmap='viridis', norm=plt.Normalize(vmin=0, vmax=1))
sm._A = []  # # fake up the array of the scalar mappable
cbaxes = fig.add_axes([0.16, 0.16, 0.03, 0.33]) 
cbar = plt.colorbar(sm, ticks=ticks, cax=cbaxes)
cbar.ax.set_yticklabels(zeds)
cbar.ax.set_ylabel('$z$', size=18, rotation=90)
    
ax.set_ylabel('$\mathrm{log_{10}}\,(\phi \,/\, \mathrm{Mpc^{-3}} \, \mathrm{dex^{-1}})$', size=16)
ax.set_xlabel('$\mathrm{log_{10}} \, (\mathrm{SFR} \,/\, \mathrm{M_{\odot}} \, \mathrm{yr^{-1}})$', size=16)

# ax.legend(frameon=False, loc=1);

# plt.show()
imgf='images/sfrf_all.png'
print(imgf)
#plt.show()
fig.savefig(imgf, dpi=150, bbox_inches='tight')

