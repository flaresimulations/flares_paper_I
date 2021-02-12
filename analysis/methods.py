
import numpy as np

def phi_flares(property_dict, bins, binLimits):

    phi_all = np.zeros(len(bins))
    phi_sigma = np.zeros((len(fl.halos),len(bins)))
    hist_all = np.zeros(len(bins))
    V_total = 0.    
    
    for i,halo in enumerate(fl.halos):
    
        w = weights[np.where(["%02d"%i == halo for i in index])[0]]
        # print(i,halo,w)
    
        _temp = property_dict[halo][tag]
        #mstar_temp = mstar_temp[mstar_temp > 0.]
    
        V = (4./3) * np.pi * R**3
        V_total += V
    
        hist, dummy = np.histogram(np.log10(_temp), bins = binLimits)
        hist = np.float64(hist)
        phi = (hist / V) / (binLimits[1] - binLimits[0])
        phi_sigma[i] = (np.sqrt(hist) / V) / (binLimits[1] - binLimits[0])
    
        phi_all += np.array(phi) * w #* V_factor
        phi_sigma[i] *= w
        hist_all += hist
    
    
    phi_sigma = np.sqrt(np.sum(np.square(phi_sigma), axis=0))

    return phi, phi_sigma, hist


def mass_bins():
    massBinLimits = np.linspace(7.95, 13.35, 28)
    massBins = np.logspace(8.05, 13.25, 27)
    # massBinLimits = np.linspace(6.85, 13.45, 23)
    # massBins = np.logspace(7.0, 13.3, 22) 
    return massBins, massBinLimits

def switch_samples(samples):
    # switch samples from bimodal distribution
    mask = samples['alpha_1'] > samples['alpha_2']

    a1_temp = samples['alpha_1'][mask]
    a2_temp = samples['alpha_2'][mask]

    samples['alpha_1'][mask] = a2_temp
    samples['alpha_2'][mask] = a1_temp

    p1_temp = samples['log10phi*_1'][mask]
    p2_temp = samples['log10phi*_2'][mask]

    samples['log10phi*_1'][mask] = p2_temp
    samples['log10phi*_2'][mask] = p1_temp

    return samples



def weighted_quantile(values, quantiles, sample_weight=None, 
                      values_sorted=False, old_style=False):
    """ Very close to numpy.percentile, but supports weights.
    NOTE: quantiles should be in [0, 1]!
    :param values: numpy.array with data
    :param quantiles: array-like with many quantiles needed
    :param sample_weight: array-like of the same length as `array`
    :param values_sorted: bool, if True, then will avoid sorting of
        initial array
    :param old_style: if True, will correct output to be consistent
        with numpy.percentile.
    :return: numpy.array with computed quantiles.
    """
    
    # do some housekeeping
    values = np.array(values)
    quantiles = np.array(quantiles)
    if sample_weight is None:
        sample_weight = np.ones(len(values))
    sample_weight = np.array(sample_weight)
    assert np.all(quantiles >= 0) and np.all(quantiles <= 1), \
        'quantiles should be in [0, 1]'

    # if not sorted, sort values array
    if not values_sorted:
        sorter = np.argsort(values)
        values = values[sorter]
        sample_weight = sample_weight[sorter]

    
    weighted_quantiles = np.cumsum(sample_weight) - 0.5 * sample_weight
    if old_style:
        # To be convenient with numpy.percentile
        weighted_quantiles -= weighted_quantiles[0]
        weighted_quantiles /= weighted_quantiles[-1]
    else:
        weighted_quantiles /= np.sum(sample_weight)
    return np.interp(quantiles, weighted_quantiles, values)


def binned_weighted_quantile(x,y,weights,bins,quantiles=[0.5]):

    if not isinstance(quantiles,list):
        quantiles = [quantiles]

    out = np.full((len(bins)-1,len(quantiles)),np.nan)

    for i,(b1,b2) in enumerate(zip(bins[:-1],bins[1:])):
        mask = (x >= b1) & (x < b2)
        if np.sum(mask) > 0:
            _temp = weighted_quantile(y[mask],quantiles,sample_weight=weights[mask])
            out[i,:] = _temp 

    return np.squeeze(np.array(out))



def piecewise_linear(x, x0, y0, m1, m2):
    """
    Fit a piecewise linear regression with two parts:
    
        y1 = m1*x + c1  (x <= x0)
        y2 = m2*x + c2  (x > x0)
    
    where,
    
        c1 = -m1*x0 + y0
        c2 = -m2*x0 + y0
        
        y0 = c2 + m2*x0 = c1 + m1*x0
        
    """
    return np.piecewise(x, [x < x0], [lambda x: m1*x + y0 - m1*x0, lambda x: m2*x + y0 - m2*x0])


def bootstrap_resample(X, n=None):
    """ Bootstrap resample an array_like
    Parameters
    ----------
    X : array_like
      data to resample
    n : int, optional
      length of resampled array, equal to len(X) if n==None
    Results
    -------
    returns resample index
    """
    if n == None: n = len(X)
    
    resample_i = np.floor(np.random.rand(n)*len(X)).astype(int)
    return resample_i


