import numpy as np
import pandas as pd

# s14 = pd.read_csv('obs_data/speagl14_highz.csv')

santini17 = pd.read_csv('obs_data/santini17.csv', delimiter=',')
santini17['mstar'] = 10**santini17['mstar'] * (0.62 / 1.06)
santini17['sfr'] = 10**santini17['sfr'] * (0.62 / 1.06)


# santini17 rebinned (see Paola private comm.)
s17_rebin = pd.DataFrame()
s17_rebin['z_low'] = [1.3,1.3,1.3,2,2,2,3,3,3,4,4,4]
s17_rebin['z_high'] = [2,2,2,3,3,3,4,4,4,5,5,5]
s17_rebin['mass'] = [8.400,9.200,9.950,8.400,9.200,9.950,8.400,9.200,9.950,8.400,9.200,9.950]
# s17_rebin['mass'] = 10**s17_rebin['mass'] * (0.62 / 1.06) # imf correction
s17_rebin['sigma'] = [0.505975,0.463054,0.253123,0.437175,0.372676,0.220388,0.351086,
                      0.367125,0.458413,0.513641,0.359223,0.458002]


salmon15 = pd.read_csv('obs_data/salmon15.csv', delimiter=',')
salmon15.columns = ['logM','logSFR','sigma_MAD','sigma_MC','z']
salmon15['logM'] -= 0.25 # Salpeter -> Chabrier
salmon15['logSFR'] -= 0.25 # Salpeter -> Chabrier
salmon15['logM'] += np.log10(pow(0.7, 2)) - np.log10(pow(0.6777, 2))
salmon15['sigma_intr'] = salmon15['sigma_MAD']**2 - salmon15['sigma_MC']**2
salmon15.loc[salmon15['sigma_intr'] <= 0.,'sigma_intr'] = \
    salmon15.loc[salmon15['sigma_intr'] <= 0.,'sigma_MC']**2
salmon15['sigma_intr'] = np.sqrt(salmon15['sigma_intr'])


s17_params = pd.DataFrame()
s17_params['z'] = [1.65, 2.5, 3.5, 4.5, 5.5]
s17_params['alpha'] = [1.04,1.16,1.02,0.94,0.92]
s17_params['alpha_err'] = [0.03,0.03,0.04,0.06,0.15]
s17_params['beta'] =[1.01,1.22,1.37,1.37,1.99]
s17_params['beta_err'] =[0.04,0.03,0.03,0.05,0.13]


s15_params = pd.read_csv('obs_data/salmon15_params.csv')
s15_params.columns = ['Name','a','a_err','b','b_err','b when a = 1','sigma_mad','chi2','z']


tasca15 = pd.DataFrame()
tasca15['z'] = [2, 3] # 0.4, 1.1
tasca15['x0'] = [1.5e10, 2.5e10] # 8e8, 1e10

shivaei15 = pd.DataFrame()
shivaei15['z'] = [2.0,2.4]
shivaei15['alpha'] = [0.65,0.58]
shivaei15['alpha_err'] = [0.08,0.10]

# IllustrisTNG
donnari19 = pd.DataFrame()
donnari19['z'] = [1.75,2]
donnari19['alpha'] = [0.77,0.77]
donnari19['alpha_err'] = [0.02,0.03]
donnari19['beta'] = donnari19['alpha']*9.7 + np.array([-6.97,-6.83])
donnari19['beta_err'] = [0.24,0.25]

donnari19_turnover = {}
donnari19_turnover['z'] = [1.75]
donnari19_turnover['turnover'] = [10.7]

whitaker14 = pd.read_csv('obs_data/whitaker14.txt', sep='\t')

whitaker14_params = pd.DataFrame()
whitaker14_params['z'] = [0.75,1.25,1.75,2.25]
whitaker14_params['alpha'] = [0.94,0.99,1.04,0.91]
whitaker14_params['alpha_err'] = [0.03,0.04,0.05,0.06]
whitaker14_params['alpha2'] = [0.14,0.51,0.62,0.67]
whitaker14_params['alpha2_err'] = [0.08,0.07,0.06,0.06]
whitaker14_params['beta'] = [1.11,1.31,1.49,1.62]
whitaker14_params['beta_err'] = [0.03,0.02,0.02,0.02]

whitaker14_turnover = pd.DataFrame()
whitaker14_turnover['z'] = [0.75,1.5,2.25]
whitaker14_turnover['turnover'] = [10.0,10.2,10.5] 
whitaker14_turnover['turnover_err'] = [0.1,0.1,0.3] 


def sargent14(m,z):
    A = [1.85,2.05,2.38]
    B = 0.16 #[0.31,0.16,0.09]
    C = 1.54 # [1.22,1.54,1.86][::-1]
    nu = -0.2
    N = [0.092,0.095,0.097]
#     N = (0.095/1e9) * 10**(nu * (m - np.log10(5e10)))
#     return N * np.exp((A*z)/(1+B*z**C)) * 10**m

    N_factor = 1e-9 * 10**(nu * (m - np.log10(5e10)))
    return [np.log10(N[i] * N_factor * np.exp((A[i]*z)/(1+B*z**C)) * 10**m) \
        for i in np.arange(3)]
    
    

def Schreiber15(m, z, imf=True):
    r = np.log10(1.+z)
    m -= 9
    m0 = 0.5 # +- 0.07
    m1 = 0.36 # +- 0.3
    a0 = 1.5 # +- 0.15
    a1 = 0.3 # +- 0.08
    a2 = 2.5 # +- 0.6
    
    if imf is True: # Salpeter -> Chabrier
        m -= 0.25
        m0 -= 0.25
    
    return m - m0 + (a0 * r) - \
            a1 * np.max([0,m - m1 - (a2 * r)])**2

    

def whitaker_2012_fit(mstar, z):
    """    
    Args:
    mstar - log_10(stellar mass)
    z - redshift
    
    Returns:
    log_10(SFR)
    """
    
    alpha = 0.7 - 0.13*z
    beta = 0.38 + 1.14*z - 0.19*z**2
    
    return alpha * (mstar - 10.5) + beta


def speagle14(Mstar, t):
    """
    Args:
    Mstar - array, log10 stellar mass
    t - age of the universe, Gyr
    """
    phi = (0.84 - 0.026 * t) * Mstar - (6.51 - 0.11 * t)
    return phi

