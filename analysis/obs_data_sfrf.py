# from pylab import *

import numpy as np
import pandas as pd

import glob
import re

# folder = '/cosma7/data/dp004/dc-love2/codes/flares_idf/analysis/obs_data/'

# simulation h
import flares
fl = flares.flares(fname='../../flares/data/flares.hdf5')
h = 0.6777


def get_obs_sfrf():
    
    obs_df = {}  
    
    obs_df['katsianis_bouwens'] = {} 
    obs_df['katsianis_bouwens']['z'] = np.array([8,8,8,8,8,8,8,8,7,7,7,7,7,7,7,7,7,7,6,6,6,6,6,6,6,6,6,6,5,5,5,5,5,5,5,5,5,5,5,5])
    obs_df['katsianis_bouwens']['log10SFR'] = np.log10([43.269 ,21.704 ,10.891 ,5.469  ,2.850  ,1.803  ,0.902  ,0.359  ,73.533 ,41.186 ,
                                          23.070 ,12.921 ,7.239  ,4.235  ,2.674  ,1.687  ,0.534  ,0.171  ,141.748 ,77.951 ,
                                          42.862 ,23.585 ,12.974 ,7.132  ,3.921  ,1.860  ,0.742  ,0.309  ,382.081 ,208.215 ,
                                          113.462 ,61.828 ,33.695 ,18.369 ,10.001 ,5.452  ,2.974  ,1.280  ,0.512  ,0.203])

    obs_df['katsianis_bouwens']['phi'] = 1e-2 * np.array([0.0010 ,0.0026 ,0.0116 ,0.0120 ,0.0662 ,0.1066 ,0.2120 ,0.5480 ,
                                                         0.0002 ,0.0062 ,0.0090 ,0.0362 ,0.0578 ,0.1224 ,0.1697 ,0.3212 ,
                                                         1.0925 ,1.5901 ,0.0004 ,0.0028 ,0.0100 ,0.0330 ,0.0598 ,0.1305 ,
                                                         0.2330 ,0.3554 ,1.2496 ,2.5517 ,0.0004 ,0.0012 ,0.0063 ,0.0189 ,
                                                         0.0495 ,0.1270 ,0.1925 ,0.2486 ,0.3900 ,0.8343 ,1.6080 ,4.5640]) 

    obs_df['katsianis_bouwens']['sigma'] = 1e-2 * np.array([0.0006 ,0.0010 ,0.0030 ,0.0050 ,0.0208 ,0.0452 ,0.0680 ,0.2080 ,
                                                           0.0004 ,0.0017 ,0.0028 ,0.0064 ,0.0114 ,0.0187 ,0.0331 ,0.0894 ,
                                                           0.2731 ,0.5499 ,0.0004 ,0.0012 ,0.0024 ,0.0047 ,0.0077 ,0.0015 ,
                                                           0.0026 ,0.0598 ,0.2581 ,0.7857 ,0.0004 ,0.0006 ,0.0015 ,0.0026 ,
                                                           0.0047 ,0.0086 ,0.0125 ,0.0175 ,0.0319 ,0.0101 ,0.0331 ,0.0133])


    # correct Salpeter->Chabrier IMF
    obs_df['katsianis_bouwens']['log10SFR'] = np.log10(10**obs_df['katsianis_bouwens']['log10SFR'] * 0.63)


    name = 'mashian'
    ## NOTE: SFRs require recalibrating by a factor of 0.63 
    ## due to the updated Kennicutt & Evans+12 calibrations
    out = {'z': [4.9,5.9,6.8,7.9],
           'log10phi*': [-3.25,-3.45,-3.67,-3.79],
           'log10SFR*': [1.75,1.62,1.54,1.31],
           'alpha': [-1.59,-1.62,-1.76,-1.79]
           }

    # correct Salpeter->Chabrier IMF
    # out['log10SFR*'] = np.log10(10**np.array(out['log10SFR*']) * 0.63)

    obs_df[name] = pd.DataFrame(out) 


    name = 'smit12'
    
    obs_df[name] = pd.read_csv('obs_data/smit12.csv',
                              delim_whitespace=True,
                              skiprows=3,
                              header=None,
                              names = ['z','log10SFR','phi','sigma'])

    # correct Salpeter->Chabrier IMF
    obs_df[name]['log10SFR'] = np.log10(10**obs_df[name]['log10SFR'] * 0.63)


    return obs_df

