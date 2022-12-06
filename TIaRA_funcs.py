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
import scipy.constants as const

M_sun = 1.989e+30
R_sun = 695700000
R_Earth = 6378100
day = 86400
AU = 1.496e+11

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

def spoc_lc_path(TIC_ID, Sectors):
    '''Generates the path to the TESS-SPOC lightcurve for a given TIC-ID and 
    sector on the warwick SCRTP system'''
    if len(str(TIC_ID)) == 16:
        tid = str(TIC_ID)
    elif len(str(TIC_ID)) < 16:
        leading = '0'*(16-len(str(TIC_ID)))
        tid = leading + str(TIC_ID)
    tid1, tid2, tid3, tid4 = textwrap.wrap(tid, 4)
    paths = np.array([], dtype=str)
    for Sector in Sectors:
        sector = 's' + '0'*(4-len(str(Sector))) + str(Sector)
        SEC = 'S' + '0'*(2-len(str(Sector))) + str(Sector)
        lc_file = ('hlsp_tess-spoc_tess_phot_'+str(tid)+'-'+str(sector)
                +'_tess_v1_lc.fits')
        path_SCRTP = os.path.join('/','storage', 'astro2', 'phsqzm', 'TESS', 
                              'SPOC_30min', str(SEC))
        path_SPOC = os.path.join('target', tid1, tid2, tid3, tid4)
        path_final = os.path.join(path_SCRTP, path_SPOC, lc_file)
        paths = np.append(paths, path_final)
    return paths

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

class star:
    def __init__(self, paths):
        lc_file = fits.open(paths[-1])
        hdul0 = lc_file[0].header
        self.id =  hdul0['OBJECT'] #TIC ID
        self.ra = hdul0['RA_OBJ'] #Right ascension
        self.dec = hdul0['DEC_OBJ'] #Declination
        self.mag = hdul0['TESSMAG'] #TESS magnitude
        self.temp = hdul0['TEFF'] #Effective temperature of star
        self.log_g = hdul0['LOGG'] #Log gravity of star
        self.mh = hdul0['MH'] #Metalicity of star
        self.rad = hdul0['RADIUS'] #Radius in solar units
        if self.rad != None and float(self.rad)>0:
            self.mass = np.power(self.rad, 1.25) #Mass of star estimated from radius using power law
        else:
            self.mass = None
        #Assign spectral type based on temperature
        if self.temp != None:
            if 7500. <= self.temp < 10000.:
                self.spectral_type = 'A'
            if 6000. <= self.temp < 7500.:
                self.spectral_type = 'F'
            if 5200. <= self.temp < 6000.:
                self.spectral_type = 'G'
            if 3700. <= self.temp < 5200.:
                self.spectral_type = 'K'
            if 2400. <= self.temp < 3700.:
                self.spectral_type = 'M'
        else:
            self.spectral_type = None
        self.lightcurves = []
        for p in paths:
            self.lightcurves.append(self.lightcurve(p))
    class lightcurve:
        def __init__(self, path):
            lc_file = fits.open(path)
            data = lc_file[1].data
            rawflux = data['PDCSAP_FLUX']
            rawfluxerr = data['PDCSAP_FLUX_ERR']
            rawtime = data['TIME']
            qual = data['QUALITY']
            good = ((qual == 0)&(np.isnan(rawflux)==False)&(np.isnan(rawtime)==False)&(np.isnan(rawfluxerr)==False))
            self.flux = rawflux[good] #PDSCAP flux with mask applied to remove points with bad quality flags and NaNs
            self.flux_err = rawfluxerr[good] #PDCSAP flux error with mask applied to remove points with bad quality flags and NaNs
            self.time = rawtime[good] #Time series in TESS BJD
            hdul0 = lc_file[0].header
            hdul1 = lc_file[1].header
            self.sector = hdul0['SECTOR'] #Sector of lightcurve
            self.cam = hdul0['CAMERA'] #Camera of lightcurve
            self.ccd = hdul0['CCD'] #CCD detector of lightcurve
            self.livetime = hdul1['LIVETIME'] #Livetime of lightcurve
            self.deadtime = hdul1['DEADC'] #Deadtime correction
            self.crowd = hdul1['CROWDSAP'] #Ratio of background flux to target flux
            self.target_frac = hdul1['FLFRCSAP'] #Fraction of target flux in aperture
            self.var = hdul1['PDCVAR'] #Variability measure
            self.cadence = hdul1['TIMEDEL']
            if hdul1['CDPP2_0'] > 0:
                self.noise2hr = hdul1['CDPP2_0'] #Two hour RMS noise measurement from fits header
                self.headernoise = True
            else:
                self.noise2hr = np.std(self.flux)*np.sqrt((1/12)*self.cadence)
                self.headernoise = False
    class planet:
        def __init__(self, params):
            self.radius = params['radius'] #Radius in Earth Radii
            self.period = params['period'] #Period in days
            self.e = params['e'] #Orbital eccentricity
            self.periastron = params['periastron'] #Periastron angle in radians
            self.probability = params['probability'] # Geometric probability of transit
            self.a = params['a'] #Orbital semi-major axis in AU
    def inject_planets(cls, number, rate_tables):
        '''Injects a given number of planets into a star based from occurence rates stored as a dictionary of pandas dataframes'''
        anglew = np.random.uniform(low=-0.5*np.pi, high=0.5*np.pi, size=number)
        e = np.random.beta(a=1.03, b=13.6, size=number)
        rate_table = rate_tables[cls.spectral_type]
        rowdex = rate_table.shape[0]
        rates = rate_table['f'].to_numpy()
        rows = np.random.choice(rowdex, p=rates, size=number)
        #This currently assumes the period bin edges are the 1st and 2nd columns of the occurence rate dataframe
        per = np.random.uniform(low=rate_table.iloc[rows, 0].to_numpy(), high=rate_table.iloc[rows, 1].to_numpy(), 
                            size=number)
        #This currently assumes the radius bin edges are the 3rd and 4th columns of the occurence rate dataframe 
        radp = np.random.uniform(low=rate_table.iloc[rows, 2].to_numpy(), high=rate_table.iloc[rows, 3].to_numpy(), 
                             size=number)
        a = np.cbrt((np.square(per*day)*const.G*(cls.mass*M_sun))/(4*np.pi))/AU
        prob = ((cls.rad*R_sun+radp*R_Earth)/a)*((1+e*np.sin(anglew))/(1-np.square(e)))
        Params = pd.DataFrame({'period':per, 'radius': radp, 'e':e, 'periastron':anglew, 'a':a, 'probability':prob})
        cls.Planets = []
        for i in range(number):
            planeti = cls.planet(Params.iloc[i])
            cls.Planets.append(planeti)
        
# fgk = pd.read_csv(os.path.join('Occurrence_rates','norm_rates_FGK.csv'))
# tables = {'F': fgk, 'G':fgk, 'K':fgk}
# target_list_path = os.path.join('target-lists', 'year1targets.csv')
# target_list = pd.read_csv(target_list_path, index_col=None)
# tic, secs = tic_choose(target_list, 127)
# Paths = spoc_lc_path(tic, secs)
# Star = star(Paths)
# Star.inject_planets(2, tables)
# print(vars(Star.Planets[0]))
# print(vars(Star.Planets[1]))

