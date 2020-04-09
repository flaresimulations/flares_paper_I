# from pylab import *

import numpy as np
import pandas as pd

import glob
import re

folder = '/cosma7/data/dp004/dc-love2/codes/flares_idf/analysis/obs_data/'

# simulation h
import flares
fl = flares.flares(fname='../../flares/data/flares.hdf5')
#h = fl.h
h = 0.6777

#directory = '/cosma5/data/Eagle/ScienceRuns/Planck1/L0100N1504/PE/REFERENCE/data'
#tag = '010_z003p984'
#h = E.readAttribute('SUBFIND', directory, tag, "/Header/HubbleParam")

def get_obs_gsmf():

    obs_df = {}
    
    name = "gonzalez_11"
    d = {'z', 'obs_M', 'obs_phi', 'obs_sigma_lower', 'obs_sigma_upper'}
    obs_df[name] = pd.DataFrame(columns = d)
    
    
    obs_h = 0.7     # observational h
    
    phi_norm = 0.67 / pow(obs_h,3) * pow(h,3)  # 0.67 factor to correct to a Chabrier IMF
    
    obs_M = 10**(np.array([7.75, 8.25, 8.75, 9.25, 9.75, 10.25, 10.75]) + np.log10(pow(obs_h, 2)) - np.log10(pow(h, 2)))
    
    
    def gonzalez_norm(obs_M, z, phi, sigma_lower, sigma_upper, phi_norm):
        obs_phi_norm = (10**phi) * phi_norm
        obs_sigma_lower = ((10**phi) - (10**(phi - sigma_lower))) * phi_norm
        obs_sigma_upper = ((10**(phi + sigma_upper)) - (10**phi)) * phi_norm
    
        d = {'z': [z] * len(obs_M), 'obs_M': obs_M, 'obs_phi': obs_phi_norm, 'obs_sigma_lower': obs_sigma_lower, 'obs_sigma_upper': obs_sigma_upper}
        return pd.DataFrame(d)
    
    
    obs_df[name] = obs_df[name].append(gonzalez_norm(obs_M, 3.8,
                                         np.array([-1.9, -1.96, -2.30, -2.53, -3.14, -3.80, -4.43]),
                                         np.array([0.12, 0.12, 0.11, 0.11, 0.13, 0.23, 0.46]),
                                         np.array([0.11, 0.12, 0.10, 0.10, 0.12, 0.20, 0.26]),
                                         phi_norm), ignore_index = True)
    
    obs_df[name] = obs_df[name].append(gonzalez_norm(obs_M, 5.0,
                                         np.array([-2.21, -2.27, -2.60, -2.84, -3.46, -4.12, -4.81]),
                                         np.array([0.20, 0.20, 0.19, 0.17, 0.21, 0.27, 0.45]),
                                         np.array([0.18, 0.18, 0.15, 0.17, 0.17, 0.22, 0.33]),
                                         phi_norm), ignore_index = True)
    
    obs_df[name] = obs_df[name].append(gonzalez_norm(obs_M, 5.9,
                                         np.array([-2.09, -2.16, -2.55, -2.81, -3.52, -4.22, -4.97]),
                                         np.array([0.24, 0.23, 0.21 , 0.22, 0.23, 0.31, 0.70]),
                                         np.array([0.23, 0.23, 0.20, 0.19, 0.21, 0.28, 0.38]),
                                         phi_norm), ignore_index = True)
    
    obs_df[name] = obs_df[name].append(gonzalez_norm(obs_M[:-1], 6.8,  # note obs_M array shorter (no high mass data)
                                         np.array([-2.15, -2.23, -2.70, -3.01, -3.80, -4.53]),
                                         np.array([0.39, 0.38, 0.36, 0.34, 0.37, 0.61]),
                                         np.array([0.41, 0.37, 0.34, 0.34, 0.36, 0.39]),
                                         phi_norm), ignore_index = True)
    
    
    
    
    name="duncan_14"
    # Data obtained directly from the author, provided by Scott Clay
    
    obs_h = 0.7     # observational h
    Duncan_files = glob.glob(folder + 'scott_data/*.cat')
    
    d = {'z', 'obs_M', 'obs_phi', 'obs_sigma_lower', 'obs_sigma_upper'}
    obs_df[name] = pd.DataFrame(columns = d)
    
    phi_norm = 1 / pow(obs_h,3) * pow(h,3)
    
    for i in Duncan_files:
    
        z = int(re.findall('[0-9]',i)[-1]) # need to change to accept two digit z
    
        GSMF = pd.read_csv(i, delim_whitespace = True, names = ['obs_M','obs_phi','obs_sigma_lower','obs_sigma_upper'])
    
        obs_M = 10**(GSMF.obs_M + np.log10(pow(obs_h, 2)) - np.log10(pow(h,2)))
    
        obs_phi_norm = GSMF.obs_phi * phi_norm
        obs_sigma_lower = GSMF.obs_sigma_lower * phi_norm * 9.99999e-1
        obs_sigma_upper = GSMF.obs_sigma_upper * phi_norm
    
        d = {'z': [z] * len(obs_M), 'obs_M': obs_M, 'obs_phi': obs_phi_norm, 'obs_sigma_lower': obs_sigma_lower, 'obs_sigma_upper': obs_sigma_upper}
        obs_df[name] = obs_df[name].append(pd.DataFrame(d), ignore_index = True)
   
   
    name='stefanon_17'
    # https://iopscience.iop.org/article/10.3847/1538-4357/aa72d8

    GSMF = pd.read_csv('obs_data/stefanon17.csv', 
                       delim_whitespace=True,
                       skiprows=2,
                       header=None,
                       names = ['z','obs_M','obs_phi','obs_sigma_lower','obs_sigma_upper'])

    GSMF['obs_M'] = 10**GSMF['obs_M']
    GSMF['obs_phi'] = GSMF['obs_phi']*1e-5
    GSMF['obs_sigma_lower'] = -(GSMF['obs_sigma_lower']*1e-5)
    GSMF['obs_sigma_upper'] = (GSMF['obs_sigma_upper']*1e-5)# - GSMF['obs_phi']

    obs_df[name] = GSMF



    name="song_15"
    # http://arxiv.org/pdf/1507.05636.pdf
    
    d = {'z', 'obs_M', 'obs_phi', 'obs_sigma_lower', 'obs_sigma_upper'}
    obs_df[name] = pd.DataFrame(columns = d)
    
    # observational h
    obs_h = 0.7
    phi_norm = 0.67 / pow(obs_h,3) * pow(h,3) # factor of 0.67 to convert to a Chabrier IMF
    
    obs_M = 10**(np.array([7.25,7.75,8.25,8.75,9.25,9.75,10.25,10.75,11.25]) + np.log10(pow(obs_h,2)) - np.log10(pow(h,2)))
    
    log_obs_phi = np.array([-1.57,-1.77,-2.00,-2.22,-2.52,-2.91,-3.41,-4.11,-5.00])
    obs_phi = 10**log_obs_phi * phi_norm
    obs_sigma_upper = 10**(log_obs_phi + np.array([0.21,0.15,0.13,0.09,0.09,0.12,0.13,0.22,0.24])) * phi_norm
    obs_sigma_upper -= obs_phi
    obs_sigma_lower = 10**(log_obs_phi - np.array([0.16,0.14,0.10,0.09,0.09,0.05,0.08,0.21,0.97])) * phi_norm
    obs_sigma_lower = obs_phi - obs_sigma_lower
    
    d = {'z': [4] * len(obs_M), 'obs_M': obs_M, 'obs_phi': obs_phi, 'obs_sigma_lower': obs_sigma_lower, 'obs_sigma_upper': obs_sigma_upper}
    obs_df[name] = obs_df[name].append(pd.DataFrame(d))
    
    log_obs_phi = np.array([-1.47,-1.72,-2.01,-2.33,-2.68,-3.12,-3.63,-4.40,-5.96])
    obs_phi = 10**log_obs_phi * phi_norm
    obs_sigma_upper = 10**(log_obs_phi + np.array([0.24,0.20,0.16,0.15,0.07,0.09,0.13,0.15,0.98])) * phi_norm
    obs_sigma_upper -= obs_phi
    obs_sigma_lower = 10**(log_obs_phi - np.array([0.21,0.20,0.16,0.10,0.14,0.11,0.11,0.35,1.98])) * phi_norm
    obs_sigma_lower = obs_phi - obs_sigma_lower
    
    d = {'z': [5] * len(obs_M), 'obs_M': obs_M, 'obs_phi': obs_phi, 'obs_sigma_lower': obs_sigma_lower, 'obs_sigma_upper': obs_sigma_upper}
    obs_df[name] = obs_df[name].append(pd.DataFrame(d))
    
    log_obs_phi = np.array([-1.47,-1.81,-2.26,-2.65,-3.14,-3.69,-4.55,-5.96,-7.50])
    obs_phi = 10**log_obs_phi * phi_norm
    obs_sigma_upper = 10**(log_obs_phi + np.array([0.35,0.23,0.21,0.15,0.12,0.12,0.19,0.52,1.30])) * phi_norm
    obs_sigma_upper -= obs_phi
    obs_sigma_lower = 10**(log_obs_phi - np.array([0.32,0.28,0.16,0.15,0.11,0.13,0.24,0.32,0.99])) * phi_norm
    obs_sigma_lower = obs_phi - obs_sigma_lower
    
    d = {'z': [6] * len(obs_M), 'obs_M': obs_M, 'obs_phi': obs_phi, 'obs_sigma_lower': obs_sigma_lower, 'obs_sigma_upper': obs_sigma_upper}
    obs_df[name] = obs_df[name].append(pd.DataFrame(d))
    
    log_obs_phi = np.array([-1.63,-2.07,-2.49,-2.96,-3.47,-4.11,-5.12,-6.57,-8.59])
    obs_phi = 10**log_obs_phi * phi_norm
    obs_sigma_upper = 10**(log_obs_phi + np.array([0.54,0.45,0.38,0.32,0.32,0.41,0.68,1.14,2.34])) * phi_norm
    obs_sigma_upper -= obs_phi
    obs_sigma_lower = 10**(log_obs_phi - np.array([0.54,0.41,0.32,0.30,0.35,0.57,1.10,1.37,3.23])) * phi_norm
    obs_sigma_lower = obs_phi - obs_sigma_lower
    
    d = {'z': [7] * len(obs_M), 'obs_M': obs_M, 'obs_phi': obs_phi, 'obs_sigma_lower': obs_sigma_lower, 'obs_sigma_upper': obs_sigma_upper}
    obs_df[name] = obs_df[name].append(pd.DataFrame(d))
    
    log_obs_phi = np.array([-1.73,-2.28,-2.88,-3.45,-4.21,-5.31,-6.93,-8.83,-12.10])
    obs_phi = 10**log_obs_phi * phi_norm
    obs_sigma_upper = 10**(log_obs_phi + np.array([1.01,0.84,0.75,0.57,0.63,1.01,1.61,2.57,4.26])) * phi_norm
    obs_sigma_upper -= obs_phi
    obs_sigma_lower = 10**(log_obs_phi - np.array([0.84,0.64,0.57,0.60,0.78,1.64,2.57,5.75,13.54])) * phi_norm
    obs_sigma_lower = obs_phi - obs_sigma_lower
    
    d = {'z': [8] * len(obs_M), 'obs_M': obs_M, 'obs_phi': obs_phi, 'obs_sigma_lower': obs_sigma_lower, 'obs_sigma_upper': obs_sigma_upper}
    obs_df[name] = obs_df[name].append(pd.DataFrame(d))
   
    return obs_df
