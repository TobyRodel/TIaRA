# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import pandas as pd
import os
import textwrap
from astropy.io import fits
    
def spoc_lc_path(TIC_ID, Sector):
    if len(str(TIC_ID)) == 16:
        tid = str(TIC_ID)
    elif len(str(TIC_ID)) < 16:
        leading = '0'*(16-len(str(TIC_ID)))
        tid = leading + str(TIC_ID)
    tid1, tid2, tid3, tid4 = textwrap.wrap(tid, 4)
    sector = 's' + '0'*(4-len(str(Sector))) + str(Sector)
    SEC = 'S' + '0'*(2-len(str(Sector))) + str(Sector)
    lc_file = 'hlsp_tess-spoc_tess_phot_'+str(tid)+'-'+str(sector)+'_tess_v1_lc.fits'
    path_SCRTP = os.path.join('/','storage', 'astro2', 'phsqzm', 'TESS', 'SPOC_30min', str(SEC))
    path_SPOC = os.path.join('target', tid1, tid2, tid3, tid4)
    path_final = os.path.join(path_SCRTP, path_SPOC, lc_file)
    return path_final


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

year1 = np.arange(1,14)
#print(y1targets)
 
    