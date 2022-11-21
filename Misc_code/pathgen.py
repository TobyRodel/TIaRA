#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 11:05:26 2022

@author: phrpkq
"""

import os
import textwrap
import numpy as np
import pandas as pd
import astropy.constants as const
from astropy.io import fits

def spoc_lc_path(TIC_ID, Sector):
    if len(str(TIC_ID)) == 16:
        tid = str(TIC_ID)
    if len(str(TIC_ID)) < 16:
        leading = '0'*(16-len(str(TIC_ID)))
        tid = leading + str(TIC_ID)
    tid1, tid2, tid3, tid4 = textwrap.wrap(tid, 4)
    sector = 's' + '0'*(4-len(str(Sector))) + str(Sector)
    SEC = 'S' + '0'*(2-len(str(Sector))) + str(Sector)
    lc_file = 'hlsp_tess-spoc_tess_phot_'+str(tid)+'-'+str(sector)+'_tess_v1_lc.fits'
    #print(lc_file)
    path_SCRTP = os.path.join('/','storage', 'astro2', 'phsqzm', 'TESS', 'SPOC_30min', str(SEC))
    path_SPOC = os.path.join('target', tid1, tid2, tid3, tid4)
    path_final = os.path.join(path_SCRTP, path_SPOC, lc_file)
    return path_final
