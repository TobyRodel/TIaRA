# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import os
from astropy.io import ascii
from scipy import constants

def planetmaker(number, rates, radius_low, radius_up, period_low, period_up):
    '''Generates parameters for a given number of planets in a planetary system from occurence rates'''
    rowdex = np.arange(len(rates))
    cosi = np.random.random_sample()
    anglew = np.random.uniform(low=-180., high=180., size=number)
    e = np.random.beta(a=1.03, b=13.6, size=number)
    rows = np.random.choice(rowdex, p=rates, size=number)
    radp = np.random.uniform(low=radius_low[rows], high=radius_up[rows], size=number)
    per = np.random.uniform(low=period_low[rows], high=period_up[rows], size=number)
    return radp, per, cosi, anglew, e

#Define constants
M_sun = 1.989e+30
R_sun = 695700000
R_Earth = 6378100
day = 86400
AU = 1.496e+11

#Generate a random temperature. In the final pipeline Temperature will be chosen from TICv8
T_eff = np.random.uniform(3700.0, 7500.0)
#Assign spectral type based on temperature
if 7500<= T_eff <10000:
    spectral_type = 'A'
if 6000<= T_eff <7500:
    spectral_type = 'F'
if 5200<= T_eff <6000:
    spectral_type = 'G'
if 3700<= T_eff < 5200:
    spectral_type = 'K'
if 2400<= T_eff < 3700:
    spectral_type = 'M'

print(spectral_type, 'type star')

#Choose mass and radius dependent on spectral type. In the final pipeline these will come from TICv8
if spectral_type == 'A':
    M_star = np.random.uniform(2.1, 1.4)
    R_star = np.random.uniform(1.4, 1.8)
if spectral_type == 'F':
    M_star = np.random.uniform(1.04, 1.4)
    R_star = np.random.uniform(1.15, 1.4)
if spectral_type == 'G':
    M_star = np.random.uniform(0.8, 1.04)
    R_star = np.random.uniform(0.96, 1.15)
if spectral_type == 'K':
    M_star = np.random.uniform(0.45, 0.8)
    R_star = np.random.uniform(0.7, 0.96)
if spectral_type == 'M':
    M_star = np.random.uniform(0.08, 0.45)
    R_star = np.random.uniform(0.1, 0.8)

print('Mass = ', M_star, 'Solar masses')
print('Radius =', R_star, 'Solar radii')
#Generate the number of planets in the system based on spectral type. These values of mu are temporary
if spectral_type in ['A', 'F', 'G', 'K']:
    n_planet = np.random.poisson(0.689)
if spectral_type == 'M':
    n_planet = np.random.poisson(2.5)

print('Hosts',n_planet, 'planets')

Path = os.path.join('..', 'Occurrence_rates', 'norm_rates_FGK.csv')
table = ascii.read(Path, format='csv', header_start=0, data_start=1)
occ = np.array(table['f'], dtype=float)
radbinlower = np.array(table['rmin'], dtype=float)
radbinupper = np.array(table['rmax'], dtype=float)
pbinlower = np.array(table['pmin'], dtype=float)
pbinupper = np.array(table['pmax'], dtype=float)

if spectral_type in ['F', 'G', 'K']:
    radii, periods, cosincs, peris, eccs = planetmaker(n_planet, occ, radbinlower, radbinupper, pbinlower, pbinupper)
    a = np.cbrt((constants.G*(M_star*M_sun)*np.square(periods*day))/(4*np.square(np.pi)))
    b = ((a*np.square(cosincs))/(R_star*R_sun))*((1-np.square(eccs))/(1+eccs*np.sin(peris)))
    k = (radii*R_Earth)/(R_star*R_sun)
    Transit = (0<= b) & (b< 1)
    Grazing = (b >= 1) & (b < (1+k))
    incs = np.degrees(np.arccos(cosincs))
    for i in range(n_planet):
        print('Planet', i)
        print('Radius = ',radii[i], 'Earth Radii')
        print('Period =', periods[i], 'days')
        print('Inclination angle =', incs, 'deg')
        print('Eccentricity =', eccs[i])
        print('Semi major axis =',a[i]/AU, 'AU')
        print('Impact parameter, b =', b[i])
        print('Transit?', Transit[i])
        print('Grazing?', Grazing[i])
        
    