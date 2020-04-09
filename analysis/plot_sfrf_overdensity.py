"""
Plot the GSMF for the whole sample, fit with a single or double schechter function

Basic template for other scripts
"""

import numpy as np
import scipy

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
# mpl.rcParams['text.usetex'] = True

import fitDF.fitDF as fitDF
import fitDF.models as models
import fitDF.analyse as analyse

binLimits = np.linspace(-3.0, 3.0, 21)
bins = np.logspace(-2.85, 2.85, 20)
print(np.log10(bins))
print(binLimits)


import flares
fl = flares.flares('../../flares/data/flares.hdf5')

R = 14./0.6777

mstar = fl.load_dataset('SFR_inst_30',arr_type="Subhalo")

## Overdensity Weights ##
dat = np.loadtxt(fl.weights, skiprows=1, delimiter=',')
log1pdelta = dat[:,5]
weights = dat[:,8]
index = dat[:,0]

dbinLims = [-0.3,-0.15,-0.075,0.0,0.1,0.2,0.28,0.3] # np.linspace(0.3,-0.3,8)
dbins = dbinLims[:-1] +  np.diff(dbinLims)/2
dbins = np.array(["%.2f"%db for db in dbins]).astype(float)


dselect = np.digitize(log1pdelta, dbinLims) - 1
dindex = np.arange(0,np.max(dselect)+1)
print(dselect)
print(dindex)

ticks = np.linspace(0.05, .95, len(dindex))
colors = [ cm.plasma(i) for i in ticks ]


## save normalisation
od_norm = np.zeros((3,len(dindex)))
print(od_norm.shape)



fig, (ax1,ax2,ax3) = plt.subplots(3,1,figsize=(5,14))

for l,(tag_idx,ax) in enumerate(zip([5,3,1],[ax1,ax2,ax3])):
    tag = fl.tags[tag_idx]
    print(tag)
    z = float(tag[5:].replace('p','.'))

    for k,(j,c) in enumerate(zip(dindex,colors)):
        print("N:",np.sum(dselect == j))
    
        phi_all = np.zeros(len(bins))
        phi_sigma = np.zeros((len(fl.halos[dselect == j]),len(bins)))
        hist_all = np.zeros(len(bins))
        V_total = 0.
    
        mstar_temp = []
    
        for i,halo in enumerate(fl.halos[dselect == j]):
            mstar_temp.append(mstar[halo][tag][mstar[halo][tag] > 0.].copy())
    
        mstar_temp = np.concatenate(mstar_temp)
            
        V = (4./3) * np.pi * R**3 * np.sum(dselect == j)
    
        phi, phi_sigma, hist = fl.calc_df(mstar_temp, tag, V, binLimits)
    
        fl.plot_df(ax, phi, phi_sigma, hist_all, massBins=bins, color=c, lines=True, label=dbins[j])

        pivot = 0.5
        mask = (np.log10(bins) < pivot+.1) & (np.log10(bins) > pivot-.1)
        if np.sum(mask) > 0:
            od_norm[l,k] = np.log10(phi[mask])[0]


    ax.text(0.8, 0.85, '$z = %.1f$'%z, transform=ax.transAxes, size=15)
    ax.set_xlim(-3,3)
    ax.set_ylim(-5.5,-.8)
    ax.grid(alpha=0.5)
   
    ax.set_ylabel('$\mathrm{log_{10}}\,(\phi \,/\, \mathrm{Mpc^{-3}} \, \mathrm{dex^{-1}})$', size=15)


sm = plt.cm.ScalarMappable(cmap='plasma', norm=plt.Normalize(vmin=0, vmax=1))
sm._A = []  # # fake up the array of the scalar mappable
cbaxes = fig.add_axes([0.3, 0.67, 0.025, 0.1])
cbar = plt.colorbar(sm, ticks=ticks, cax=cbaxes)
cbar.ax.set_yticklabels(dbins)
cbar.ax.set_ylabel('$\mathrm{log_{10}}(1 \,+\,\delta)$', size=15, rotation=90)

ax3.set_xlabel('$\mathrm{log_{10}} \, (\mathsf{SFR} \,/\, \mathrm{M_{\odot}} \, \mathrm{yr^{-1}})$', size=15)

fname='images/sfrf_overdensity.png'
print(fname)
plt.savefig(fname, dpi=150, bbox_inches='tight')
# plt.show()
plt.close()


## overdensity fits ##
fig, ax = plt.subplots(1,1)
ax.plot(dbins,od_norm[0,:],color='C0')
p = np.polyfit(dbins[~np.isinf(od_norm[0,:])],
                       od_norm[0,~np.isinf(od_norm[0,:])],deg=1,full=False)
print(p)
ax.plot(dbins,dbins*p[0] + p[1],color='C0',ls='dashed')

ax.plot(dbins,od_norm[1,:],color='C1')
p = np.polyfit(dbins[~np.isinf(od_norm[1,:])],
                       od_norm[1,~np.isinf(od_norm[1,:])],deg=1,full=False)
print(p)
ax.plot(dbins,dbins*p[0] + p[1],color='C1',ls='dashed')

ax.plot(dbins,od_norm[2,:],color='C2')
p = np.polyfit(dbins[~np.isinf(od_norm[2,:])],
                       od_norm[2,~np.isinf(od_norm[2,:])],deg=1,full=False)
print(p)
ax.plot(dbins,dbins*p[0] + p[1],color='C2',ls='dashed')

# plt.show()
                                       
