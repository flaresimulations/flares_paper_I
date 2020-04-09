"""
Fit SFS
"""
import json
import numpy as np
import scipy

import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcParams['text.usetex'] = True

from scipy import optimize

binLimits = np.linspace(0, 3.2, 33)
bins = np.logspace(0.05, 3.15, 32)
print(np.log10(bins))
print(binLimits)

from methods import piecewise_linear,bootstrap_resample

import flares
fl = flares.flares('../../flares/data/flares.hdf5')
sfr = fl.load_dataset('SFR_inst_30')
mstar = fl.load_dataset('Mstar_30')

## Overdensity Weights ##
dat = np.loadtxt(fl.weights, skiprows=1, delimiter=',')
weights = dat[:,8]
index = dat[:,0]


def fit(tag, N=1000, n=int(1e2)): 

    mstar_all = []
    sfr_all = []
    w_all = []

    for i,halo in enumerate(fl.halos):
        w = weights[np.where(["%02d"%i == halo for i in index])[0]]

        mstar_temp = mstar[halo][tag]#[centrals[halo][tag]]
        sfr_temp = sfr[halo][tag]#[centrals[halo][tag]]

        mask = (mstar_temp > 0.) & (sfr_temp > 1e-3)

        w_all = np.append(w_all,np.repeat(w,np.sum(mask)))
        mstar_all = np.append(mstar_all,mstar_temp[mask])
        sfr_all = np.append(sfr_all,sfr_temp[mask])


    
    x = np.log10(mstar_all / 10**9.7)
    y = np.log10(sfr_all)
    p0 = [.2,1.5,1.1,.6]

    # do the bootstrap fit
    p_temp = [None] * N
    errs_temp = [None] * N
    for i in np.arange(N):
        bs_1 = bootstrap_resample(x,n)
        bs_2 = bootstrap_resample(y,n)
        p_temp[i], errs_temp[i] = optimize.curve_fit(piecewise_linear,x,y,p0=p0)
    
    print(errs_temp[0])

    x0 = [p[0] for p in p_temp]
    y0 = [p[1] for p in p_temp]
    m1 = [p[2] for p in p_temp]
    m2 = [p[3] for p in p_temp]
    
    #p_temp = np.array(p_temp)
    #errs_temp = np.array(errs_temp)
    print(np.diag(errs_temp[0]))

    p = {}
    for j,(name,param) in enumerate(zip(['x0','y0','m1','m2'],[x0,y0,m1,m2])):
        p[name] = {}
        p[name]['median'] = np.median(param).tolist()
        p[name]['mean'] = np.mean(param).tolist()
        #p[name]['std'] = np.array(param).std().tolist()
        p[name]['std'] = np.sqrt(np.diag(errs_temp[0]))[j]

   
    with open('samples/sfs_fit_%s.json'%tag, 'w') as outfile:
        json.dump(p, outfile)

    return p 



tags = list(fl.tags)

from schwimmbad import MultiPool
from functools import partial

with MultiPool() as pool:
    values = list(pool.map(fit, tags))
    print(values)
