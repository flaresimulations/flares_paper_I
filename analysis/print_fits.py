
import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcParams['text.usetex'] = True

import fitDF.fitDF as fitDF
import fitDF.models as models
import fitDF.analyse as analyse

model = models.DoubleSchechter()

import flares
fl = flares.flares(fname='')

_format = 'latex'

print("\n###################\nGSMF:\n####################\n")
for tag in fl.tags:
    z = int(float(tag[5:].replace('p','.')))

    a = analyse.analyse(ID='samples', model=model, 
                        sample_save_ID='flares_gsmf_%s'%tag,verbose=False)

    _params = {}
    for ip, p in enumerate(a.parameters):
        _params[p] = np.zeros(3)
        _params[p][1] = np.round(np.percentile(a.samples[p], 50),3)
        _params[p][0] = np.round(_params[p][1] - np.percentile(a.samples[p], 16),3)
        _params[p][2] = np.round(np.percentile(a.samples[p], 84) - _params[p][1],3)

    # print(_params)

    print(f"{ z } &\
            ${_params['D*'][1]}_{{ - {_params['D*'][0]} }}^{{ +{_params['D*'][2]} }}$ &\
            ${_params['log10phi*_1'][1]}_{{-{_params['log10phi*_1'][0]}}}^{{+{_params['log10phi*_1'][2]} }}$ &\
            ${_params['log10phi*_2'][1]}_{{-{_params['log10phi*_2'][0]}}}^{{+{_params['log10phi*_2'][2]} }}$ &\
            ${_params['alpha_1'][1]}_{{-{_params['alpha_1'][0]}}}^{{+{_params['alpha_1'][2]}}}$ \\\\\
            ")



print("\n###################\nSFRF:\n###################\n")
for tag in fl.tags:#[[5,6,7,8,9,10]]: 
    a = analyse.analyse(ID='samples', model=model, 
                        sample_save_ID='flares_sfrf_%s'%tag,verbose=False)
   

    _params = {}
    for ip, p in enumerate(a.parameters):
        _params[p] = np.zeros(3)
        _params[p][1] = np.round(np.percentile(a.samples[p], 50),3)
        _params[p][0] = np.round(_params[p][1] - np.percentile(a.samples[p], 16),3)
        _params[p][2] = np.round(np.percentile(a.samples[p], 84) - _params[p][1],3)

    # print(_params)

    print(f"{ z } &\
            ${_params['D*'][1]}_{{ - {_params['D*'][0]} }}^{{ +{_params['D*'][2]} }}$ &\
            ${_params['log10phi*_1'][1]}_{{-{_params['log10phi*_1'][0]}}}^{{+{_params['log10phi*_1'][2]} }}$ &\
            ${_params['log10phi*_2'][1]}_{{-{_params['log10phi*_2'][0]}}}^{{+{_params['log10phi*_2'][2]} }}$ &\
            ${_params['alpha_1'][1]}_{{-{_params['alpha_1'][0]}}}^{{+{_params['alpha_1'][2]}}}$ \\\\\
            ")


