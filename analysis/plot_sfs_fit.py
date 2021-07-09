
"""
Plot fitted GSMF and original data
"""
import json
import sys
import numpy as np


import matplotlib as mpl
from matplotlib import cm
from matplotlib import gridspec
import matplotlib.pyplot as plt

from methods import piecewise_linear

from methods import mass_bins, binned_weighted_quantile
exec(open("./obs_data_sfs.py").read())

from astropy.cosmology import Planck15 as cosmo

import flares
fl = flares.flares(fname='../../flares/data/flares.hdf5')

tags = fl.tags
zeds = np.array([float(tag[5:].replace('p','.')) for tag in tags])


## ---- Plot
fig = plt.figure(figsize=(5,10), dpi=150)

gs = gridspec.GridSpec(7,1)
gs.update(wspace=0.24, hspace=0.1)

ax1 = plt.subplot(gs[5:7,0])
ax2 = plt.subplot(gs[:3,0])
# ax3 = plt.subplot(gs[2])
ax4 = plt.subplot(gs[3:5,0])

x0,y0,m1,m2 = [],[],[],[]
x0_err,y0_err,m1_err,m2_err = [],[],[],[]

for tag in fl.tags:
    print(tag)
    # try:
    with open('samples/sfs_fit_%s.json'%tag) as f:
        p = json.load(f)
    # except:
    #     print(tag)
    #     continue

    x0.append(p['x0']['median'])
    y0.append(p['y0']['median'])
    m1.append(p['m1']['median'])
    m2.append(p['m2']['median'])

    x0_err.append(p['x0']['std'])
    y0_err.append(p['y0']['std'])
    m1_err.append(p['m1']['std'])
    m2_err.append(p['m2']['std'])


x0 = np.array(x0)
y0 = np.array(y0)
m1 = np.array(m1)
m2 = np.array(m2)

print(zeds)
print("x0",x0 + 9.7)
print("alpha_1",m1)
print("alpha_2",m2)
print("beta",-1 * m1 * x0 + y0)

ax2.errorbar(zeds, m1, yerr=m1_err)
             #label=label, 
             #color=C, marker=None, linestyle=linesty) # m1
ax2.set_ylabel('$\\alpha \;\; (1,2)$', size=15)
ax2.set_ylim(-.3, 1.45)

ax2.errorbar(zeds[1:], m2[1:], yerr=m2_err[1:])#, label=label, 
         #color=C, marker=None, linestyle=linesty) # m2

ax1.errorbar(zeds[1:], x0[1:] + 9.7, yerr=x0_err[1:], color='C1')#, label=label, 
         #color=C, marker=None, linestyle=linesty)
ax1.set_ylabel('$x_{0} + 9.7$', size=15)    

ax4.errorbar(zeds[1:], (-1 * m1 * x0 + y0)[1:], yerr=m1_err[1:],color='C1')
#, label=label, color=C, marker=None) # c1
ax4.set_ylabel('$\\beta$', size=15)

#_norm = (-1 * m1[0] * x0[0] + y0[0])[1:]
ax4.plot(zeds[1:], 1. / (cosmo.age(zeds[1:])), color='C2', label='$1 / t \; (\mathrm{Gyr})$')
ax4.legend(frameon=False,loc=4)

## Observations ##
s = 7

s15_mask = (s15_params['Name'] == 'Salmon15')
salmon15_artist = ax2.errorbar(s15_params['z'][s15_mask], s15_params['a'][s15_mask], 
                s15_params['a_err'][s15_mask], label='$\mathrm{Salmon\mathsf{+}15}$\n$[9.0-\infty]$', fmt='p', color='grey', ms=s)
ax4.errorbar(s15_params['z'][s15_mask], 9.7 * s15_params['a'][s15_mask] + s15_params['b'][s15_mask], 
         color='grey', ms=s, fmt='p')

s15_mask = (s15_params['Name'] == 'Behroozi et al. 2013b ')
behroozi13_artist = ax2.errorbar(s15_params['z'][s15_mask], s15_params['a'][s15_mask], s15_params['a_err'][s15_mask], 
                          label='$\mathrm{Behroozi\mathsf{+}13}$\n$(-\infty)$', fmt='>', color='grey', ms=s)
ax4.errorbar(s15_params['z'][s15_mask], 9.7 * s15_params['a'][s15_mask] + s15_params['b'][s15_mask], 
         color='grey', fmt='>', ms=s)


# wmask = whitaker14_params['z'] > 1.2
# whitaker14_artist = ax2.errorbar(whitaker14_params['z'][wmask], whitaker14_params['alpha'][wmask], 
#              whitaker14_params['alpha_err'][wmask], 
#              label='Whitaker+14 \n[($\sim$9.0-9.4)->$\infty$]', fmt='X', color='forestgreen', ms=s)
# ax2.errorbar(whitaker14_params['z'][wmask], whitaker14_params['alpha2'][wmask], 
#              whitaker14_params['alpha2_err'][wmask], 
#              label='Whitaker+14', fmt='X', color='forestgreen', ms=s, markerfacecolor='none')
# ax1.errorbar(whitaker14_turnover['z'], whitaker14_turnover['turnover'], 
#              whitaker14_turnover['turnover_err'], 
#              label='Whitaker+14', fmt='X', color='forestgreen', ms=s, markerfacecolor='none')


shivaei15_artist = ax2.errorbar(shivaei15['z'], shivaei15['alpha'], 
             shivaei15['alpha_err'], label='$\mathrm{Shivaei\mathsf{+}15}$\n$[9.5-11.5]$', fmt='*', color='grey', ms=s)

# donnari19_artist = ax2.errorbar(donnari19['z'], donnari19['alpha'], 
#              donnari19['alpha_err'], label='Donnari+19 \n[9->10.5]', fmt='D', color='yellowgreen', ms=s)
# ax4.errorbar(donnari19['z'], donnari19['beta'], 
#              donnari19['beta_err'], fmt='D', color='yellowgreen', ms=s)
# ax1.errorbar(donnari19_turnover['z'], donnari19_turnover['turnover'], 
#             fmt='D', color='yellowgreen', ms=s)


schreiber15_artist = ax2.plot([1.5,6],[1.0,1.0],linestyle='dotted', 
                              color='grey', label='$\mathrm{Schreiber\mathsf{+}15}$\n$[(\sim9-11)-\infty]$')
ax2.errorbar([4.25],[0.8],xerr=0.75, color='grey',fmt='*', 
             label='Schreiber+15', capsize=10, ms=s)

# tasca15_artist = ax1.errorbar(tasca15['z'], np.log10(tasca15['x0']), label='Tasca+15', 
#                               color='forestgreen', fmt='o', ms=s, markerfacecolor='none')

santini17_artist = ax2.errorbar(s17_params['z'], s17_params['alpha'], fmt='s', 
                               label='$\mathrm{Santini\mathsf{+}17}$\n$[(8.3-8.7)-\infty]$', color='grey', ms=s)
ax4.errorbar(s17_params['z'], s17_params['beta'], s17_params['beta_err'], 
            fmt='s', color='grey', ms=s)

z = np.linspace(1,6,100) # [float(tag[5:].replace('p','.')) for tag in tags]
t = cosmo.age(z).value
speagle14_artist = ax2.plot(z,0.84-0.026*t, linestyle='-.', 
                            color='grey', label='$\mathrm{Speagle\mathsf{+}14}$\n$[9.7-11.1]$')
ax4.plot(z,(0.84-0.026*t)*9.7 - (6.51 - 0.11 * t), linestyle='-.', color='grey')

# s14_artist = ax2.scatter(s14['zmed'], s14['alpha'], color='forestgreen', s=2, 
#                          label='Speagle+14\n(compilation)')
# ax4.scatter(s14['zmed'], s14['alpha']*9.7 + s14['beta'], 
#                           color='forestgreen', s=2)# s14['beta_err'], , fmt='.')#, label='Speagle+14\n(compilation)')

# s15_mask = (s15_params['Name'] == 'Lu et al. ')
# ax2.errorbar(s15_params['z'][s15_mask], s15_params['a'][s15_mask], s15_params['a_err'][s15_mask], 
#                           label='Lu+13', color='yellowgreen')
# ax4.plot(s15_params['z'][s15_mask], 9.7 * s15_params['a'][s15_mask] + s15_params['b'][s15_mask], 
#          color='yellowgreen')

# s15_mask = (s15_params['Name'] == 'Davé et al. 2013 ')
# ax2.errorbar(s15_params['z'][s15_mask], s15_params['a'][s15_mask], s15_params['a_err'][s15_mask], 
#                           label='Dave+13', color='yellowgreen')
# ax4.scatter(s15_params['z'][s15_mask], 9.7 * s15_params['a'][s15_mask] + s15_params['b'][s15_mask], 
#             color='yellowgreen')

# s15_mask = (s15_params['Name'] == 'Somerville et al. ')
# somerville15_artist = ax2.errorbar(s15_params['z'][s15_mask], s15_params['a'][s15_mask], s15_params['a_err'][s15_mask], 
#                           label='Somerville+12', color='yellowgreen', fmt='<')
# ax4.scatter(s15_params['z'][s15_mask], 9.7 * s15_params['a'][s15_mask] + s15_params['b'][s15_mask], 
#          color='yellowgreen', marker='<')

## END Observations ##

ax1.set_xlabel('$z$', size=15)

for ax in [ax2,ax4]: ax.set_xticklabels([])

for ax in [ax1,ax2,ax4]: # ,ax3
    ax.grid(axis='x', alpha=0.4)
    ax.set_xlim(2.8,10.2)


from matplotlib.lines import Line2D

hi_line = Line2D([0],[0], color='C1', label='$\\alpha_{2}$')
lo_line = Line2D([0],[0], color='C0', label='$\\alpha_{1}$')
#ax2.legend(handles=[lo_line,hi_line], fontsize=14, frameon=False)

handles = [santini17_artist,salmon15_artist,schreiber15_artist[0],shivaei15_artist,speagle14_artist[0],behroozi13_artist,hi_line,lo_line]
labels = [h.get_label() for h in handles]

"""
Remove error bars from legend
https://swdg.io/2015/errorbar-legends/
"""
new_handles = []
for h in handles:
    #only need to edit the errorbar legend entries
    if isinstance(h, mpl.container.ErrorbarContainer): 
        new_handles.append(h[0])
    else: new_handles.append(h)

legend1 = ax2.legend(handles=new_handles, labels=labels, 
                     ncol=3, frameon=False)#, bbox_to_anchor=[1.16,.5,0.5,0.5])
ax2.add_artist(legend1)

plt.show()
fname = 'images/sfs_fit_evolution.pdf'
print(fname)
# fig.savefig(fname, dpi=300, bbox_inches='tight')
