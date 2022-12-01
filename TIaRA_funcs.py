#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  1 11:37:05 2022

@author: phrpkq
"""

import os
import numpy as np
import pandas as pd
import textwrap
from astropy.io import fits

def tarlist_compile(sectors):
    '''Converts TESS-SPOC target lists for a given list of sectors into a 
    dictionary'''
    data = np.array(list())
    secdict = {}
    for s in sectors:
        leading = '0'*(4-len(str(s)))
        sid ='s'+leading+str(s)
        dat = pd.read_csv(os.path.join('..', 'target-lists',sid+'.csv'),
                          usecols=['#TIC_ID'])
        secdict[sid] = dat.to_numpy(dtype=str)
        data = np.append(data, dat.to_numpy(dtype=str))
    ticlist = np.unique(data)
    return ticlist, secdict

def sector_search(tics, secs):
    '''Takes a unique list of TIC IDs and a dictionary of targets for different
    sectors and outputs a pandas dataframe with a True/False value for whether 
    each target is observed in each sector'''
    df = pd.DataFrame({'OBJECT':tics})
    for sec in secs.keys():
        secindex = np.isin(tics, secs[sec])
        df[sec] = secindex
    return df

def tic_choose(targets, row):
    '''Outputs a target id and a list of sectors from a given row of a pandas
    dataframe (when loading the csv to a dataframe set index_col=0)'''
    tic_id = str(targets.iloc[row]['OBJECT'])
    sector_bool = np.array(targets.iloc[row], dtype=bool)[1:]
    sectors = np.arange(1, targets.shape[1])
    return tic_id, sectors[sector_bool]

def spoc_lc_path(TIC_ID, Sector):
    '''Generates the path to the TESS-SPOC lightcurve for a given TIC-ID and 
    sector on the warwick SCRTP system'''
    if len(str(TIC_ID)) == 16:
        tid = str(TIC_ID)
    elif len(str(TIC_ID)) < 16:
        leading = '0'*(16-len(str(TIC_ID)))
        tid = leading + str(TIC_ID)
    tid1, tid2, tid3, tid4 = textwrap.wrap(tid, 4)
    sector = 's' + '0'*(4-len(str(Sector))) + str(Sector)
    SEC = 'S' + '0'*(2-len(str(Sector))) + str(Sector)
    lc_file = ('hlsp_tess-spoc_tess_phot_'+str(tid)+'-'+str(sector)
               +'_tess_v1_lc.fits')
    path_SCRTP = os.path.join('/','storage', 'astro2', 'phsqzm', 'TESS', 
                              'SPOC_30min', str(SEC))
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

def planetmaker(number, rates, radius_low, radius_up, period_low, period_up):
    '''Generates parameters for a given number of planets in a planetary system
    from occurence rates'''
    rowdex = np.arange(len(rates))
    anglew = np.random.uniform(low=-90., high=90., size=number)
    e = np.random.beta(a=1.03, b=13.6, size=number)
    rows = np.random.choice(rowdex, p=rates, size=number)
    radp = np.random.uniform(low=radius_low[rows], high=radius_up[rows], 
                             size=number)
    per = np.random.uniform(low=period_low[rows], high=period_up[rows], 
                            size=number)
    return radp, per, anglew, e