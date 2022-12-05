import os
import numpy as np
import pandas as pd
import textwrap
from astropy.io import fits
import matplotlib.pyplot as plt

import TIaRA_funcs

tex_fonts = {
    # Use LaTeX to write all text
    "text.usetex": True,
    "font.family": "serif",
    # Use 10pt font in plots, to match 10pt font in document
    "axes.labelsize": 13,
    "font.size": 12,
    # Make the legend/label fonts a little smaller
    "legend.fontsize": 13,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12
}

plt.rcParams.update(tex_fonts)

target_lists1 = pd.read_csv(os.path.join('target-lists', 's0001.csv'))
tics = target_lists1['#TIC_ID'].to_numpy(dtype=str)
rows = np.arange(low=0, high=target_lists1.shape[0])
print(rows)
noise = []
headnoise = []
mag = []

for i in rows:
    path = TIaRA_funcs.spoc_lc_path(tics[i], [1])
    star = TIaRA_funcs.star(path[0])
    lc = star.lightcurve(path[0])
    noise.append(lc.noise2hr)
    headnoise.append(lc.headernoise)
    mag.append(star.mag)

header_noise = np.array(noise)[headnoise]
head_mag = np.array(mag)[headnoise]
calc_noise = np.array(noise)[np.invert(headnoise)]
calc_mag = np.array(mag)[np.invert(headnoise)]
anomalies = (np.array(noise)>5000)

fig, ax = plt.subplots(1, 1)
ax.scatter(x=head_mag, y=header_noise, s=1, label='header noise')
ax.scatter(x=calc_mag, y=calc_noise, s=1, label='Manually calculated noise')
ax.set_xlabel('TESS magnitude')
ax.set_ylim(0, 5000)
ax.set_ylabel('Noise')
ax.set_title('Sector 1 noise')
plt.text(4,2500,('Number of anomalies=',anomalies.sum))
plt.text(4,2000, ('number of missing heads = ', len(calc_noise)))
plt.gca().invert_xaxis()
plt.legend()
plt.savefig('noiseVmag.pdf', format='pdf')