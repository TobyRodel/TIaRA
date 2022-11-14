# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import fits

def LC_read(path):
    '''Reads a lightcurve and returns useful metadata in a pandas dataframe
    and timeseries data'''
    lc = fits.open(path)
    head1 = lc[0].header
    head2 = lc[1].header
    good = (lc[1].data['QUALITY'] == 0)
    f = lc[1].data['PDCSAP_FLUX']
    e = lc[1].data['PDCSAP_FLUX_ERR']
    t = lc[1].data['TIME']
    df = {'TIC-ID':[head1['OBJECT']], 'Sector':[head1['SECTOR']], 
      'Camera':[head1['CAMERA']], 'Chip':[head1['CCD']], 
      'RA':[head1['RA_OBJ']], 'DEC':[head1['DEC_OBJ']],
      'T-MAG':[head1['TESSMAG']], 'T_Eff':[head1['TEFF']],
      'Log-g':[head1['LOGG']], 'Metalicity':[head1['MH']],
      'R_star':[head1['RADIUS']], 'M_star':[np.power(head1['Radius'], 1.25)], 
      'Livetime':[head2['LIVETIME']], 'deadtime-cor':[head2['DEADC']],
      'RMS-Noise':[head2['CDPP2_0']], 'Crowd':[head2['CROWDSAP']],
      'Target-Frac':[head2['FLFRCSAP']], 'Variability':[head2['PDCVAR']]}
    log = pd.DataFrame(df)
    flux = f[good]
    time = t[good]
    fluxerr = e[good]
    return log, flux, time, fluxerr

Path = 'lightcurves\hlsp_tess-spoc_tess_phot_0000000033910247-s0001_tess_v1_lc.fits'

LOG, FLUX, TIME, FLUXERR = LC_read(Path)

print(LOG)
print(LOG['R_star'])
print(LOG['M_star'])

plt.errorbar(x=TIME, y=FLUX, yerr=FLUXERR)
plt.xlabel('Time (TESS BJD)')
plt.ylabel('Flux ($e^- s^{-1}$)')
plt.show()   
    