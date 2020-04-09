import glob
import os
import re
import numpy as np


def get_sim_gsmf():
    simd = {}
    simd['fire'] = {}
    
    simd['fire'][5] = {'logMstar': [3.83,4.50,5.17,5.84,6.51,7.18,7.85,8.52,9.19,9.85],
               'logphi':   [1.10,0.95,0.67,0.29,-0.15,-0.61,-1.34,-2.09,-2.32,-3.87]}
    
    simd['fire'][6] = {'logMstar': [3.83, 4.48, 5.13, 5.78, 6.43, 7.08, 7.74, 8.39, 9.04, 9.69],
               'logphi':   [1.20,0.95,0.68,0.22,-0.25,-0.83,-1.39,-2.08,-2.68,-3.87]}
    
    simd['fire'][7] = {'logMstar': [3.81,4.42,5.04,5.65,6.26,6.88,7.49,8.11,8.72,9.33],
               'logphi':   [1.29,0.93,0.63,0.15,-0.34,-0.93,-1.43,-2.10,-2.82,-3.91]}
    
    simd['fire'][8] = {'logMstar': [3.86, 4.57, 5.29, 6.00, 6.72, 7.43, 8.15, 8.86],
               'logphi':   [1.19 ,  0.79 , 0.29 , -0.32, -1.05, -1.65, -2.45, -3.76]} 
    
    simd['fire'][9] = {'logMstar': [3.83, 4.50, 5.17, 5.84, 6.50, 7.17, 7.84, 8.51],
               'logphi':   [1.08 ,  0.69 , 0.18 , -0.42, -1.14, -1.78, -2.52, -3.56]} 
    
    simd['fire'][10] = {'logMstar': [3.80, 4.41, 5.01, 5.62, 6.22, 6.83, 7.43, 8.04],
               'logphi':   [0.98 ,  0.61 , 0.15 , -0.49, -1.16, -1.77, -2.36, -3.34]} 
    
    simd['fire'][11] = {'logMstar': [3.82, 4.46, 5.10, 5.74, 6.38, 7.03, 7.67],
               'logphi':   [0.81 ,  0.48 , -0.18, -0.86, -1.67, -2.63, -2.75]} 
    
    simd['fire'][12] = {'logMstar': [3.81, 4.42, 5.04,  5.65,  6.26,  6.88,  7.49],
               'logphi':   [0.65 ,  0.33 , -0.27, -1.11, -1.78, -2.66, -3.13]} 
    
    ## Yung+19
    files = glob.glob('sim_data/yung19/SMF*.dat')
    simd['yung19'] = {}
    
    for f in files:
        z = re.findall(r"[-+]?\d*\.\d+|[-+]?\d+", os.path.basename(f))[0]
        dat = np.loadtxt(f,skiprows=3)
        simd['yung19'][z] = {'logMstar': dat[:,0], 'logphi': np.log10(dat[:,1])}


    simd['Henriques+15'] = {}
    simd['Henriques+20'] = {}
    lgals_dir = '/cosma7/data/dp004/dc-irod1/FLARES/'
    snapnums = [11, 12, 13, 15, 17, 19]
    zeds = [9.72, 8.93, 8.22, 6.97, 5.92, 5.03]
    for z,snapnum in zip(zeds,snapnums):
        dat = np.load(lgals_dir+'SMF/data_' + str(snapnum) + '_2015.npy')
        simd['Henriques+15'][z] = {'logMstar': dat[0], 'logphi': np.log10(dat[1])}
        dat = np.load(lgals_dir+'SMF/data_' + str(snapnum) + '_2020.npy')
        simd['Henriques+20'][z] = {'logMstar': dat[0], 'logphi': np.log10(dat[1])}


    simd['Wilkins+17'] = {}
    mstar = [8.1,8.3,8.5,8.7,8.9,9.1,9.3,9.5,9.7,9.9,10.1,10.3]
    simd['Wilkins+17'][10] = {'logMstar': mstar[:8], 
            'logphi': [-3.72,-4.01,-4.30,-4.61,-4.96,-5.34,-5.65,-6.11]}
    simd['Wilkins+17'][9] = {'logMstar': mstar[:10], 
            'logphi': [-3.22,-3.45,-3.71,-3.98,-4.27,-4.58,-4.93,-5.28,-5.65,-6.07]}
    simd['Wilkins+17'][8] = {'logMstar': mstar, 
            'logphi': [-2.76,-2.97,-3.19,-3.42,-3.66,-3.92,-4.20,-4.50,-4.87,-5.17,-5.59,-5.96]}


    return simd
