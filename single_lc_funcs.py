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

class star:
    def __init__(self, path):
        self.issues = 0
        try:
            hdul0 = fits.getheader(filename=path, ext=0)
            self.st_id =  hdul0['OBJECT'] #TIC ID
            self.ra = hdul0['RA_OBJ'] #Right ascension
            self.dec = hdul0['DEC_OBJ'] #Declination
            self.mag = hdul0['TESSMAG'] #TESS magnitude
            self.t_eff = hdul0['TEFF'] #Effective temperature of star
            self.log_g = hdul0['LOGG'] #Log gravity of star
            self.mh = hdul0['MH'] #Metalicity of star
            self.st_rad = hdul0['RADIUS'] #Radius in solar units
            if self.st_rad != None and float(self.st_rad)>0:
                self.mass = np.power(self.st_rad, 1.25) #Mass of star estimated from radius using power law
            else:
                self.mass = None
            #Assign spectral type based on temperature
            if self.t_eff != None:
                if 7500. <= self.t_eff < 10000.:
                    self.spectral_type = 'A'
                if 6000. <= self.t_eff < 7500.:
                    self.spectral_type = 'F'
                if 5200. <= self.t_eff < 6000.:
                    self.spectral_type = 'G'
                if 3700. <= self.t_eff < 5200.:
                    self.spectral_type = 'K'
                if 2400. <= self.t_eff < 3700.:
                    self.spectral_type = 'M'
            else:
                self.spectral_type = None
            self.lc = self.lightcurve(path)
        except FileNotFoundError:
            self.st_id = None
            self.ra = None
            self.dec = None
            self.mag = None
            self.t_eff = None
            self.log_g = None
            self.mh = None
            self.st_rad = None
            self.mass = None
            self.lc = None
            self.issues += 1
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
            self.cadence2hr = (np.mean(np.diff(self.time)))*12.
            if hdul1['CDPP2_0'] > 0:
                self.noise2hr = hdul1['CDPP2_0'] #Two hour RMS noise measurement from fits header
                self.headernoise = True
            else:
                self.noise2hr = None
                self.headernoise = False

class planets(star):
    def __init__(self, path, number, rate_tables):
            super().__init__(path=path)
            self.pl_anglew = np.random.uniform(low=-0.5*np.pi, high=0.5*np.pi, size=number) #Periastron angle
            self.pl_e = np.random.beta(a=1.03, b=13.6, size=number) #Orbital eccentricity
            rate_table = rate_tables[self.spectral_type]
            rowdex = rate_table.shape[0]
            rates = rate_table['f'].to_numpy()
            rows = np.random.choice(rowdex, p=rates, size=number)
            #This currently assumes the period bin edges are the 1st and 2nd columns of the occurence rate dataframe
            self.pl_period = np.random.uniform(low=rate_table.iloc[rows, 0].to_numpy(), high=rate_table.iloc[rows, 1].to_numpy(), size=number) #Orbital period in days
            #This currently assumes the radius bin edges are the 3rd and 4th columns of the occurence rate dataframe 
            self.pl_rad = np.random.uniform(low=rate_table.iloc[rows, 2].to_numpy(), high=rate_table.iloc[rows, 3].to_numpy(), size=number) #Planet radius in Earth Radii
            self.a = np.cbrt((np.square(self.pl_period*day)*const.G*(self.mass*M_sun))/(4*np.pi)) #Semi major axis in m
            self.p_transit = ((self.st_rad*R_sun+self.pl_rad*R_Earth)/self.a)*((1+self.pl_e*np.sin(self.pl_anglew))/(1-np.square(self.pl_e))) # Probability of transit
            self.k = ((self.pl_rad*R_Earth)/(self.st_rad*R_sun)) #Calcutae ratio of R_pl to R_star
            self.p_obs = (self.lc.time[1]-self.lc.time[0])/self.pl_period #Calculate probability of planet being observed
            #Begin calculating TSM
            self.pl_T_est = self.t_eff*np.sqrt(self.st_rad*R_sun/self.a)*(0.25**0.25) # Estimated T_eq of planet with zero Albedo and perfect heat distribution
            #Estimate mass using a power law, seperate for lower and higher mass planets
            Mrange = (self.pl_rad > 1.23)
            self.pl_M_est = np.empty_like(self.pl_rad)
            self.pl_M_est[Mrange] = 1.436*np.power(self.pl_rad[Mrange],1.7)
            self.pl_M_est[np.invert(Mrange)] = 0.9718*np.power(self.pl_rad[np.invert(Mrange)],3.58)
            #Calculate scale factor for the 4 different radius bins
            sf1 = (self.pl_rad <= 1.5)
            sf2 = ((self.pl_rad > 1.5)&(self.pl_rad<=2.75))
            sf3 = ((self.pl_rad>2.75)&(self.pl_rad<=4.))
            sf4 = (self.pl_rad>4.)
            scale = np.empty_like(self.pl_rad)
            scale[sf1] = 0.190
            scale[sf2] = 1.26
            scale[sf3] = 1.28
            scale[sf4] = 1.15
            self.TSM = scale*(np.power(self.pl_rad,3)*self.pl_T_est)/(self.pl_M_est*self.st_rad**2)*np.power(10., self.mag*-0.2)
class signals(planets):
    def __init__(self,path, number, rate_tables, N_b, N_phase):
        super().__init__(path=path, number=number, rate_tables=rate_tables)
        self.N_pl = number
        self.N_b = N_b
        self.N_ph = N_phase
        self.N_sig = N_b*N_phase
        self.N_tot = self.N_sig*self.N_pl
        #self.b = np.random.uniform(low=0, high=(1.+ self.k), size=(number, N_b)) #Calculate impact parameter b
        self.b = np.random.uniform(low=0., high=1., size=(number, N_b)) #Don't generate grazing transits
        self.cosi = np.empty_like(self.b)
        self.T_dur = np.empty_like(self.b)
        for i in range(number):
            self.cosi[i] = (self.st_rad*R_sun*self.b[i]*(1.+self.pl_e[i]*np.sin(self.pl_anglew[i])))/(self.a[i]*(1.-np.square(self.pl_e[i]))) # Calculate cosi from b
            self.T_dur[i] = (self.pl_period[i]*day/np.pi)*np.arcsin((self.st_rad*R_sun/self.a[i])*np.sqrt(1+self.k[i]-np.square(self.b[i]))/np.sqrt(1-np.square(self.cosi[i]))) # Duration of transit in seconds
        #self.cosi = (self.st_rad*R_sun*self.b*(1.+self.pl_e*np.sin(self.pl_anglew)))/(self.a*(1.-np.square(self.pl_e))) 
        T_min = self.lc.time[0]-(self.T_dur/day)
        T_max = self.lc.time[-1]+(self.T_dur/day)
        self.T_0 = np.empty(shape=(number, N_b, N_phase))
        self.ph_offset = np.empty_like(self.T_0)
        for i in range(number):
            self.T_0[i] = np.random.uniform(low=T_min[i], high=T_max[i], size=(N_b, N_phase)) # Create T_0 to ensure transit
            self.ph_offset[i] = ((self.T_0[i]-T_min[i]) / self.pl_period[i])%1
    def Signal_to_noise(self, per, a, b, cosi, k, T_0):
        T_dur = (per*day/np.pi)*np.arcsin((self.st_rad*R_sun/a)*np.sqrt(1+k-np.square(b))/np.sqrt(1-np.square(cosi)))
        phi = ((self.lc.time - T_0)%per)/per
        In_transit = (phi >= 0) & (phi <= ((T_dur/day)/per))
        T_trans = In_transit.sum()*(self.lc.cadence2hr)
        SN = (np.square(k)*1000000)*(np.sqrt(T_trans*0.5)/self.lc.noise2hr)
        return SN
    def get_results(self):
        self.SN = np.empty_like(self.T_0)
        iters = np.shape(self.T_0)
        for i in range(iters[0]):
            for j in range(iters[1]):
                for k in range(iters[2]):
                    self.SN[i,j,k] = self.Signal_to_noise(per=self.pl_period[i], a=self.a[i] ,b=self.b[i,j], cosi=self.cosi[i,j], k=self.k[i], T_0=self.T_0[i,j,k])
    def to_pandas(self):
        '''Places relvant object parameters in a pandas dataframe'''
        #Save values to np arrays to prepare for saving to pd dataframe:
        tic_id = np.repeat(self.st_id, self.N_tot)
        ra = np.repeat(self.ra, self.N_tot)
        dec = np.repeat(self.dec, self.N_tot)
        mag = np.repeat(self.mag, self.N_tot)
        t_eff = np.repeat(self.t_eff, self.N_tot)
        spec_type = np.repeat(self.spectral_type, self.N_tot)
        p_obs = np.repeat(self.p_obs, self.N_sig)
        t_eq = np.repeat(self.pl_T_est, self.N_sig)
        TSM = np.repeat(self.TSM, self.N_sig)
        p_transit = np.repeat(self.p_transit, self.N_sig)
        radius = np.repeat(self.pl_rad, self.N_sig)
        period = np.repeat(self.pl_period, self.N_sig)
        b = np.repeat(self.b, self.N_ph)
        offset = np.ndarray.flatten(self.ph_offset)
        SN = np.ndarray.flatten(self.SN)
        df = pd.DataFrame({'TIC_ID':tic_id, 'RA':ra, 'DEC':dec, 'TESSMAG':mag, 'ST_T_EFF':t_eff, 'TYPE':spec_type, 'PL_PERIOD':period, 'PL_RADIUS':radius, 'P_TRANS':p_transit, 'P_OBS':p_obs, 'IMPACT':b, 'OFFSET':offset,'PL_T_EST':t_eq, 'TSM':TSM,  'S/N':SN})
        return df


