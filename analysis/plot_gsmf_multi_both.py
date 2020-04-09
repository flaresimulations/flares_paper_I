"""
Plot fitted GSMF and original data
"""
import sys
import numpy as np

import matplotlib as mpl
from matplotlib import cm
import matplotlib.pyplot as plt
# mpl.rcParams['text.usetex'] = True
# plt.style.use('classic')

import fitDF.fitDF as fitDF
import fitDF.models as models
import fitDF.analyse as analyse

from obs_data_gsmf import get_obs_gsmf
from sim_data_gsmf import get_sim_gsmf
obs_df = get_obs_gsmf()
simd = get_sim_gsmf()

import flares
fl = flares.flares('')

tags = fl.tags
ticks = np.linspace(0.05, .95, len(tags))
colors = [ cm.viridis(i) for i in ticks ]

# print(tag)
zeds = [float(tag[5:].replace('p','.')) for tag in tags]

model = models.DoubleSchechter()

massBinLimits = np.linspace(7.45, 13.25, 60)
massBins = np.logspace(7.50, 13.2, 59)


fig, ((ax1,ax2,ax3),(ax4,ax5,ax6)) = plt.subplots(2,3,figsize=(11,10))

plt.subplots_adjust(hspace=0.05, wspace=0.05)

axes = [ax1,ax2,ax3,ax4,ax5,ax6]

for c,ax,tag,z in zip(colors,axes,tags,zeds):

    ## ---- Get fit
    sample_ID = 'flares_gsmf_%s'%(tag)
    a = analyse.analyse(ID='samples', model=model, sample_save_ID=sample_ID, verbose=False)
    
    model.update_params(a.median_fit)
    ax.plot(np.log10(massBins), a.model.log10phi(np.log10(massBins)), color=c, lw=3)

    ax.text(0.74, 0.91, '$z = %.1f$'%z, transform=ax.transAxes, size=14)
    ax.set_xlim(8.0, 11.9)
    ax.set_ylim(-10, -1.1) 
    ax.grid(alpha=0.5)
    
    ## ---- Observations
    for author,ms in zip(obs_df.keys(),['o','d','s','<']):
        
        mask = (obs_df[author]['z'] < (z + 0.5)) & (obs_df[author]['z'] > (z - 0.5))
        
        if np.sum(mask) > 0:
            phi = obs_df[author]['obs_phi'][mask]
            lo = np.log10(phi) - np.log10(phi - obs_df[author]['obs_sigma_lower'][mask])
            hi = np.log10(phi + obs_df[author]['obs_sigma_upper'][mask]) - np.log10(phi)
           
            ax.errorbar(np.log10(obs_df[author]['obs_M'][mask]), 
                        np.log10(obs_df[author]['obs_phi'][mask]), 
                        yerr=[lo,hi], ls='none',marker=ms,color='grey',
                        markeredgecolor='black', markeredgewidth=0.3,
                        label=author.replace('_','-'))


    ## ---- Models
    for author,ls in zip(simd.keys(),['solid','dashed','-.','dotted',(0,(3,1,1,1,1,1))]):
        zeds = np.array(list(simd[author].keys()))
        mask = (zeds.astype(float) < (z + 0.5)) & (zeds.astype(float) > (z - 0.5))

        if np.sum(mask) > 0:
            zed = zeds[mask][0]
            ax.plot(simd[author][zed]['logMstar'],
                    simd[author][zed]['logphi'],
                    ls=ls,marker='None',color='grey',label=author.replace('_','-'))





for ax in [ax1,ax4]:
    ax.set_ylabel('$\mathrm{log_{10}}\,(\phi \,/\, \mathrm{Mpc^{-3}} \, \mathrm{dex^{-1}})$', size=14)

for ax in [ax4,ax5,ax6]:
    ax.set_xlabel('$\mathrm{log_{10}} \, (M_{*} \,/\, M_{\odot})$', size=14)

for ax in [ax5,ax6,ax2,ax3]:
    ax.set_yticklabels([])

for ax in [ax1,ax2,ax3]:
    ax.set_xticklabels([])


#ax6.legend(frameon=False, loc=3);
from matplotlib.lines import Line2D
line_gonz = Line2D([0], [0], color='grey', linestyle='none', marker='o')
line_dunc = Line2D([0], [0], color='grey', linestyle='none', marker='d')
line_stef = Line2D([0], [0], color='grey', linestyle='none', marker='s')
line_song = Line2D([0], [0], color='grey', linestyle='none', marker='<')
ax5.legend([line_gonz,line_dunc,line_song,line_stef],
           ['$\mathrm{Gonzalez\mathsf{+}11}$','$\mathrm{Duncan\mathsf{+}14}$',
            '$\mathrm{Song\mathsf{+}15}$','$\mathrm{Stefanon\mathsf{+}17}$'], frameon=False, loc=3,);

line_yung = Line2D([0], [0], color='grey', linestyle='dashed')
line_fire = Line2D([0], [0], color='grey', linestyle='solid')
line_lgals15 = Line2D([0], [0], color='grey', linestyle='-.')
line_lgals20 = Line2D([0], [0], color='grey', linestyle='dotted')
line_btides = Line2D([0], [0], color='grey', linestyle=(0,(3,1,1,1,1,1)))
ax6.legend([line_btides,line_yung,line_fire,line_lgals15,line_lgals20],
           ['$\mathrm{Wilkins\mathsf{+}17}$','$\mathrm{Yung\mathsf{+}19}$','$\mathrm{Ma\mathsf{+}17}$',
            '$\mathrm{Henriques\mathsf{+}15}$','$\mathrm{Henriques\mathsf{+}20}$'], frameon=False, loc=3,);



# plt.show()
imgf='images/gsmf_multi_both.png'
print(imgf)
fig.savefig(imgf, dpi=150, bbox_inches='tight')
