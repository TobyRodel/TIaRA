import os
import numpy as np
import pandas as pd
import textwrap
from astropy.io import fits

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

class star:
    def __init__(self, filename):
        lc_file = fits.open(filename)
        hdul0 = lc_file[0].header
        self.id =  hdul0['OBJECT'] #TIC ID
        self.ra = hdul0['RA_OBJ'] #Right ascension
        self.dec = hdul0['DEC_OBJ'] #Declination
        self.mag = hdul0['TESSMAG'] #TESS magnitude
        self.temp = hdul0['TEFF'] #Effective temperature of star
        self.log_g = hdul0['LOGG'] #Log gravity of star
        self.mh = hdul0['MH'] #Metalicity of star
        self.rad = hdul0['RADIUS'] #Radius in solar units
        self.mass = np.power(self.rad, 1.25) #Mass of star estimated from radius using power law
        #Assign spectral type based on temperature
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
    class lightcurve:
        def __init__(self, filename):
            lc_file = fits.open(filename)
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
            else:
                self.noise2hr = np.std(self.flux)*np.sqrt((1/12)*self.cadence)
