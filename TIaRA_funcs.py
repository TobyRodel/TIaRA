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

def LC_read(paths):
    '''Takes an array of file paths and returns stellar and sector parameters in pandas dataframes and timeseries data as numpy arrays'''
    lc = fits.open(paths[0])
    h1 = lc[0].header
    
    #Read in stellar parameters
    ID = h1['object']
    RA = h1['RA_OBJ']
    DEC = h1['DEC_OBJ']
    MAG = h1['TESSMAG']
    TEMP = h1['TEFF']
    LOGG = h1['LOGG']
    MH = h1['MH']
    RSTAR = h1['RADIUS']
    #Estimate Mass of star, needs updating to a better model
    MSTAR = np.power(RSTAR, 1.25)

    #Place parameters into dataframe
    Star_params = pd.DataFrame({'TICID':ID, 'RA':RA, 'DEC':DEC,
        'MAG':MAG, 'TEMP':TEMP, 'LOGG':LOGG, 'MH':MH, 'RSTAR':RSTAR, 'MSTAR':MSTAR})
    
    #Create empty arrays to populate with sector data in loop
    cams = np.array([])
    chips = np.array([])
    ltimes = np.array([])
    dtimes = np.array([])
    noises = np.array([])
    crwds = np.array([])
    tfracs = np.array([])
    vars = np.array([])
    #Create empty arrays to populate with timeseries data
    flux = np.array([])
    fluxerr = np.array([])
    time = np.array([])
    #Loop through lightcurves for sector data
    for path in paths:
        lc = fits.open(path)
        h1 = lc[0].header
        h2 = lc[1].header

        f = lc[1].data['PDCSAP_FLUX']
        e = lc[1].data['PDCSAP_FLUX_ERR']
        t = lc[1].data['TIME']

        #Create mask for time series data, filtering out marked systematics and NaNs
        good = ((lc[1].data['QUALITY'] == 0)&(np.isnan(f)==False)&(np.isnan(t)==False)&(np.isnan(e)==False))
        #Apply mask
        F = f[good]
        T = t[good]
        E = e[good]

        #Read in sector data
        cam = h1['CAMERA']
        chip = h1['CCD']
        ltime = h2['LIVETIME']
        dtime = h2['DEADC']
        noise2hr = h2['CDPP2']
        crwd = h2['CROWDSAP']
        tfrac = h2['FLFRCSAP']
        var = h2['PDCVAR']

        #Update numpy arrays
        cams = np.append(cams, cam)
        chips = np.append(chips, chip)
        ltimes = np.append(ltimes, ltime)
        dtimes = np.append(dtimes, dtime)
        noises = np.append(noises, noise2hr)
        cwrds = np.append(crwds, crwd)
        tfracs = np.append(tfracs, tfrac)
        vars = np.append(vars, var)
        #Update timeseries arrays
        flux = np.append(flux, F, axis=0)
        fluxerr = np.append(fluxerr, E, axis=0)
        time = np.append(time, T, axis=0)
    #Create dataframe of sector data
    sectors = pd.DataFrame({'CAM': cams, 'CCD': chips, 'LIVETIME':ltimes, 'DEADTIMECORR': dtimes, 'NOISE':noises, 'CROWD': crwds, 'TFRAC': tfracs, 'VAR': vars})
    return Star_params, sectors, time, flux, fluxerr

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