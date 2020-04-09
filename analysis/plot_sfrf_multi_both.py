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

from sim_data_sfrf import get_sim_sfrf
from obs_data_sfrf import get_obs_sfrf
simd = get_sim_sfrf()
obs_df = get_obs_sfrf()

import flares
fl = flares.flares('../../flares/data/flares.hdf5')

tags = fl.tags
ticks = np.linspace(0.05, .95, len(tags))
colors = [ cm.viridis(i) for i in ticks ]


zeds = [float(tag[5:].replace('p','.')) for tag in tags]

model = models.DoubleSchechter()

# binLimits = np.linspace(0, 3.2, 17)
# bins = np.logspace(0.1, 3.1, 16)
binLimits = np.linspace(-3.0, 3.0, 21)
bins = np.logspace(-2.85, 2.85, 20)
print(np.log10(bins))
print(binLimits)


## ---- Plot GSMF
fig, ((ax1,ax2,ax3),(ax4,ax5,ax6)) = plt.subplots(2,3,figsize=(11,10))

plt.subplots_adjust(hspace=0.05, wspace=0.05)

axes = [ax1,ax2,ax3,ax4,ax5,ax6]

for c,ax,tag,z in zip(colors,axes,tags,zeds):

    ## ---- Get fit
    sample_ID = 'flares_sfrf_%s'%(tag)
    a = analyse.analyse(ID='samples', model=model, sample_save_ID=sample_ID, verbose=False)
    
    model.update_params(a.median_fit)
    x = np.linspace(-0.3,3,100)
    ax.plot(x, a.model.log10phi(x), color=c, lw=3)

    ax.text(0.74, 0.91, '$z = %.1f$'%z, transform=ax.transAxes, size=14)
    ax.set_xlim(-1,3)
    ax.set_ylim(-8,-1)
    ax.grid(alpha=0.5)
    
    ## ---- Models
    #for author,ls in zip(simd.keys(),['solid','dashed']):
    author = 'yung19'
    ls='dashed'
    zeds = np.array(list(simd[author].keys()))
    #.astype(float)
    mask = (zeds.astype(float) < (z + 0.5)) & (zeds.astype(float) > (z - 0.5))

    if np.sum(mask) > 0:
        zed = zeds[mask][0]
        ax.plot(simd[author][zed]['logSFR'], 
                simd[author][zed]['logphi'], 
                ls=ls,marker='None',color='grey',label='Yung+19')

    
    author = 'bluetides'
    ls='solid'
    zeds = np.array(simd[author]['z']).astype(float)
    mask = (zeds < (z + 0.5)) & (zeds > (z - 0.5))

    if np.sum(mask) > 0:
        z = zeds[mask][0]
        ax.plot(simd[author]['log10SFR'], 
                simd[author][z], 
                ls=ls,marker='None',color='grey',label='Bluetides')


    zeds = np.array(list(simd['lgals15'].keys())).astype(float)
    mask = (zeds < (z + 0.5)) & (zeds > (z - 0.5))

    if np.sum(mask) > 0:
        author = 'lgals15'
        z = zeds[mask][0]
        ax.plot(simd[author][z]['log10SFR'],
                simd[author][z]['logphi'],
                ls='-.',marker='None',color='grey',label='L-Galaxies+15')

        author = 'lgals20'
        z = zeds[mask][0]
        ax.plot(simd[author][z]['log10SFR'],
                simd[author][z]['logphi'],
                ls='dotted',marker='None',color='grey',label='L-Galaxies+20')


    ## ---- Observations
    for author,label,ms in zip(['smit12','katsianis_bouwens'],['Smit+12','Katsianis+17'],['o','d','s','<']):

        mask = (obs_df[author]['z'] < (z + 0.5)) & (obs_df[author]['z'] > (z - 0.5))

        if np.sum(mask) > 0:
            phi = obs_df[author]['phi'][mask]
            lo = np.log10(phi) - np.log10(phi - obs_df[author]['sigma'][mask])
            hi = np.log10(phi + obs_df[author]['sigma'][mask]) - np.log10(phi)

            ax.errorbar(obs_df[author]['log10SFR'][mask],
                        np.log10(obs_df[author]['phi'][mask]),
                        yerr=[lo,hi], ls='none',marker=ms,
                        markeredgecolor='black', markeredgewidth=0.3,
                        color='grey',label=label)



for ax in [ax1,ax4]:
    ax.set_ylabel('$\mathrm{log_{10}}\,(\phi \,/\, \mathrm{Mpc^{-3}} \, \mathrm{dex^{-1}})$', size=14)

for ax in [ax4,ax5,ax6]:
    ax.set_xlabel('$\mathrm{log_{10}} \, (\mathsf{SFR} \,/\, \mathrm{M_{\odot}} \, \mathrm{yr^{-1}})$', size=14)

for ax in [ax5,ax6,ax2,ax3]:
    ax.set_yticklabels([])

for ax in [ax1,ax2,ax3]:
    ax.set_xticklabels([])

from matplotlib.lines import Line2D
line_smit = Line2D([0], [0], color='grey', linestyle='none', marker='o')
line_kats = Line2D([0], [0], color='grey', linestyle='none', marker='d')
line_yung = Line2D([0], [0], color='grey', linestyle='dashed')
line_blue = Line2D([0], [0], color='grey', linestyle='solid')
line_lgals15 = Line2D([0], [0], color='grey', linestyle='-.')
line_lgals20 = Line2D([0], [0], color='grey', linestyle='dotted')
ax6.legend([line_smit,line_kats,line_yung,line_blue,line_lgals15,line_lgals20], 
           ['$\mathrm{Smit\mathsf{+}12}$',
            '$\mathrm{Katsianis\mathsf{+}17}$',
            '$\mathrm{Yung\mathsf{+}19}$',
            '$\mathrm{Wilkins\mathsf{+}17}$',
            '$\mathrm{Henriques\mathsf{+}15}$',
            '$\mathrm{Henriques\mathsf{+}20}$'], frameon=False, loc=3);

# plt.show()
imgf='images/sfrf_multi_both.png'
print(imgf)
fig.savefig(imgf, dpi=150, bbox_inches='tight')

