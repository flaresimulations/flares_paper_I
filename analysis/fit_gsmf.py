"""
Fit GSMF for all redshifts and all simulations
"""

import sys
import numpy as np
import scipy

from emcee.autocorr import integrated_time

import fitDF.fitDF as fitDF
import fitDF.models as models
import fitDF.analyse as analyse

from schwimmbad import MultiPool
from functools import partial

from methods import mass_bins
massBins, massBinLimits = mass_bins()
# print(np.log10(massBins))
# print(massBinLimits)

import flares
fl = flares.flares('../../flares/data/flares.hdf5')
mstar = fl.load_dataset('Mstar_30',arr_type='Subhalo')

## Load Periodic ##
fl_ref = flares.flares(fname='../../flares/data/EAGLE_REF_sp_info.hdf5', sim_type='PERIODIC')
# mstar_ref = fl_ref.load_dataset(name='Mstar_30')
# 
# fl_agn = flares.flares(fname='../../flares/data/EAGLE_AGNdT9_sp_info.hdf5', sim_type='PERIODIC')
# mstar_agn = fl_agn.load_dataset(name='Mstar_30')

## Overdensity Weights ##
dat = np.loadtxt(fl.weights, skiprows=1, delimiter=',')
weights = dat[:,8]
index = dat[:,0]

# if len(weights) != len(fl.halos):
#     print("Warning! number of weights not equal to number of halos")


def fitdf(N_up, N, V, mstar_temp, cprior, name):

    obs = [{'bin_edges': massBinLimits, 'N': N_up, 'volume': V, 'sigma': N_up / np.sqrt(N)}] 
    # print("upscaled N", N_up,
    #       "\nN", N,
    #       "\nsigma", obs[0]['sigma'])

    model = models.DoubleSchechter()

    priors = {}

    # scale = 2.0
    # priors['log10phi*_1'] = scipy.stats.norm(loc=cprior['phi1'], scale=scale)
    # priors['log10phi*_2'] = scipy.stats.norm(loc=cprior['phi2'], scale=scale)
    priors['log10phi*_1'] = scipy.stats.uniform(loc=-8, scale=6.0)
    priors['log10phi*_2'] = scipy.stats.uniform(loc=-8, scale=6.0)
    
    # scale = 2.0
    # priors['alpha_1'] = scipy.stats.norm(loc=cprior['a1'], scale=scale)
    # priors['alpha_2'] = scipy.stats.norm(loc=cprior['a2'], scale=scale)
    priors['alpha_1'] = scipy.stats.uniform(loc=-4.5, scale=3.0)
    priors['alpha_2'] = scipy.stats.uniform(loc=-1.001, scale=0.002)
    
    priors['D*'] = scipy.stats.uniform(loc = 8., scale = 4.0)
 
    fitter = fitDF.fitter(obs, model=model, priors=priors, output_directory='samples')
    fitter.lnlikelihood = fitter.gaussian_lnlikelihood
    samples = fitter.fit(nsamples=int(1e4), burn=1000, sample_save_ID=name, 
                             use_autocorr=True, verbose=True)

    # from methods import switch_samples
    # _samples = switch_samples(samples)
 
    observations = [{'volume':V, 'sample': np.log10(mstar_temp), 'bin_edges': massBinLimits, 'N': N_up}]
    a = analyse.analyse(ID='samples', model=model, sample_save_ID=name, observations=observations)#, samples=_samples)

    for ip, p in enumerate(a.parameters):
        acorr = integrated_time(a.samples[p], quiet=True)
        print("Autocorrelation time %s:"%p,acorr)

    #fig = a.triangle(hist2d = True, ccolor='0.5')
    #plt.savefig('images/%s_posteriors.png'%name)
 
    #fig = a.LF(observations=True,xlabel='$\mathrm{log_{10}}(M_{*} \,/\, M_{\odot})$')
    #plt.savefig('images/%s_fit_MF.png'%name)
 


def fit(tags): 
    
    tag, rtag = tags
    print(tag,rtag)

    phi_all = np.zeros(len(massBins))
    phi_sigma = np.zeros((len(fl.halos),len(massBins)))
    hist_all = np.zeros(len(massBins))
    # V_total = 0.
    R = 14./0.6777

    for i,halo in enumerate(fl.halos):

        w = weights[np.where(["%02d"%i == halo for i in index])[0]]
        
        mstar_temp = mstar[halo][tag]
        mstar_temp = mstar_temp[mstar_temp > 0.]

        V = (4./3) * np.pi * R**3
        # V_total += V
        
        #phi, phi_sigma, hist = fl.calc_df(mstar_temp, tag, V, massBinLimits)
        hist, dummy = np.histogram(np.log10(mstar_temp), bins = massBinLimits)
        hist = np.float64(hist)
        
        phi = (hist / V) / (massBinLimits[1] - massBinLimits[0])
        #phi_sigma[i] = (np.sqrt(hist) / V) / (massBinLimits[1] - massBinLimits[0])

        phi_all += np.array(phi) * w
        #phi_sigma[i] *= w
        hist_all += hist
  
   
    
    #phi_sigma = (fl.poisson_confidence_interval(hist_all.astype(int),0.68) / V_total) /\
     #           (massBinLimits[1] - massBinLimits[0])

    # phi_sigma = np.sqrt(np.sum(np.square(phi_sigma), axis=0))

    V = (3200)**3
    N = models.phi_to_N(phi_all,V,massBinLimits)

    fitdf(N, hist_all, V, mstar_temp, cprior=custom_priors[tag], name='flares_gsmf_%s'%tag)
   
    # ## Ref ##
    # # if len(mstar_ref[rtag]) > 0:
    # Phi, phi_sigma, hist = fl.calc_df(mstar_ref[rtag], rtag, 100**3, massBinLimits)
    # fitdf(hist, 100**3, mstar_ref[rtag] * 1e10, cprior=custom_priors[tag], name='ref_gsmf_%s'%rtag)

    ## AGN ##
    # # if len(mstar_agn[rtag]) > 0:
    # Phi, phi_sigma, hist = fl.calc_df(mstar_agn[rtag], rtag, 50**3, massBinLimits)
    # fitdf(hist, 50**3, mstar_agn[rtag] * 1e10, cprior=custom_priors[tag], name='agn_gsmf_%s'%rtag)
    
    return None



# moving priors
custom_priors = {'005_z010p000': {'phi1':-5.0,'phi2':-5.0,'a1':-2.0,'a2':-1.0}, 
                 '006_z009p000': {'phi1':-6.0,'phi2':-5.0,'a1':-2.0,'a2':-1.0},
                 '007_z008p000': {'phi1':-5.0,'phi2':-4.5,'a1':-2.0,'a2':-1.0},
                 '008_z007p000': {'phi1':-5.0,'phi2':-4.5,'a1':-2.0,'a2':-1.0},
                 '009_z006p000': {'phi1':-5.0,'phi2':-4.0,'a1':-2.0,'a2':-1.0},
                 '010_z005p000': {'phi1':-4.5,'phi2':-3.5,'a1':-2.0,'a2':-1.0}} 


tags = list(zip(fl.tags, fl_ref.tags))

_idx = int(sys.argv[1])
fit((fl.tags[_idx],fl_ref.tags[_idx]))

#with MultiPool(processes=1) as pool:
    # values = list(pool.map(fit, tags))
#    values = list(pool.map(fit, [(fl.tags[0],fl_ref.tags[0])]))

