import glob
import os
import re
import numpy as np


def get_sim_sfrf():
    simd = {}

    simd['bluetides'] = {}
    simd['bluetides']['z'] = [13,12,11,10,9,8]
    simd['bluetides']['log10SFR'] = [-0.5 ,-0.3 ,-0.1 ,0.1 ,0.3 ,0.5 ,0.7 ,0.9 ,1.1 ,1.3 ,1.5 ,1.7 ,1.9 ,2.1 ,2.3 ,2.5]
    simd['bluetides'][13] = [-3.93 ,-4.21 ,-4.54 ,-4.92 ,-5.28 ,-5.67 ,-6.04 ,-6.4 ,None ,None ,None ,None ,None ,None ,None ,None]
    simd['bluetides'][12] = [-3.5 ,-3.77 ,-4.07 ,-4.42 ,-4.72 ,-5.13 ,-5.49 ,-5.82 ,-6.23 ,None ,None ,None ,None ,None ,None ,None]
    simd['bluetides'][11] = [-3.11 ,-3.36 ,-3.62 ,-3.93 ,-4.22 ,-4.54 ,-4.89 ,-5.21 ,-5.63 ,-5.98 ,-6.32 ,None ,None ,None ,None ,None]
    simd['bluetides'][10] = [-2.8 ,-3.02 ,-3.26 ,-3.52 ,-3.78 ,-4.06 ,-4.39 ,-4.69 ,-5.05 ,-5.41 ,-5.77 ,-6.21 ,None ,None ,None ,None]
    simd['bluetides'][9] = [-2.51 ,-2.71 ,-2.93 ,-3.16 ,-3.39 ,-3.64 ,-3.9 ,-4.14 ,-4.44 ,-4.79 ,-5.12 ,-5.44 ,-5.91 ,-6.34 ,None ,None]
    simd['bluetides'][8] = [-2.22 ,-2.41 ,-2.61 ,-2.8 ,-3.01 ,-3.22 ,-3.46 ,-3.7 ,-3.97 ,-4.25 ,-4.57 ,-4.91 ,-5.28 ,-5.65 ,-6.12 ,-6.57]


    ## Yung+19
    files = glob.glob('sim_data/yung19/SFRF*.dat')
    simd['yung19'] = {}
    
    for f in files:
        z = re.findall(r"[-+]?\d*\.\d+|[-+]?\d+", os.path.basename(f))[0]
        dat = np.loadtxt(f,skiprows=3)
        simd['yung19'][z] = {'logSFR': dat[:,0], 'logphi': np.log10(dat[:,1])}


    simd['lgals15'] = {}
    simd['lgals20'] = {}
    lgals_dir = '/cosma7/data/dp004/dc-irod1/FLARES/'
    snapnums = [11, 12, 13, 15, 17, 19]
    zeds = [9.72, 8.93, 8.22, 6.97, 5.92, 5.03]
    for z,snapnum in zip(zeds,snapnums):
        dat = np.load(lgals_dir+'SFRF/data_' + str(snapnum) + '_2015_SFRF.npy')
        simd['lgals15'][z] = {'log10SFR': dat[0], 'logphi': np.log10(dat[1])}
        dat = np.load(lgals_dir+'SFRF/data_' + str(snapnum) + '_2020_SFRF.npy')
        simd['lgals20'][z] = {'log10SFR': dat[0], 'logphi': np.log10(dat[1])}


    return simd
