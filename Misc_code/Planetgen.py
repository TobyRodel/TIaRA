# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii

def planetmaker(number, rates, radius_low, radius_up, period_low, period_up):
    '''Generates parameters for a given number of planets in a planetary system'''
    sample = np.arange(number)
    cosi = np.random.random_sample()
    radp = np.empty_like(sample)
    per = np.empty_like(sample)
    anglew = np.empty_like(sample)
    e = np.empty_like(sample)
    for i in sample:
        rows = np.arange(len(rates))
        r = np.random.choice(rows, p=rates)
        radp[i] = np.random.uniform(low=radius_low[r]*1000, high=radius_up[r]*1000)
        per[i] = np.random.uniform(low=period_low[r]*1000, high=period_up[r]*1000)
        anglew[i] = np.random.uniform(-90., 90.)
        e[i] = np.random.beta(1.03, 13.6)
        print('Planet', i)
        print('Radius =', radp[i]/1000, 'chosen from bin:', radius_low[r], radius_up[r])
        print('Period =', per[i]/1000, 'chosen from bin:', period_low[r], period_up[r])
        print('Periastron angle = ', anglew[i])
        print('Inclination angle = ', np.degrees(np.arccos(cosi)))
        print('eccentricity = ', e[i])
    return radp/1000, per/1000, cosi, anglew, e

table = ascii.read('occurrence_rates - norm_rates_FGK.csv', format='csv', header_start=0, data_start=1)
occ = np.array(table['f'], dtype=float)
radbinlower = np.array(table['rmin'], dtype=float)
radbinupper = np.array(table['rmax'], dtype=float)
pbinlower = np.array(table['pmin'], dtype=float)
pbinupper = np.array(table['pmax'],dtype=float)

radii, periods, incs, peris, eccs = planetmaker(10, occ, radbinlower, radbinupper, pbinlower, pbinupper)

fig = plt.figure(figsize=[12.8, 9.6])
plt.xscale('symlog', base=2)
plt.yscale('symlog', base=2)
plt.scatter(periods, radii, s=1, c='#FF0000')
plt.xlabel('Period (days)')
plt.xlim(0.5, 500)
plt.ylim(0.5, 16)
plt.ylabel('Radius of planet ($R_E$)')
plt.xticks([0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256, 500], [0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256, 500])
plt.yticks([0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.5, 3, 4, 6, 8, 12, 16], [0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.5, 3, 4, 6, 8, 12, 16])
plt.grid()
#plt.savefig('Simulated_planet_pop.pdf', format='pdf')
plt.show()