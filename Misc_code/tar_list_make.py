# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 15:45:18 2022

@author: Student
"""

import os
import numpy as np
import pandas as pd

def tarlist_compile(sectors):
    data = np.array(list())
    secdict = {}
    for s in sectors:
        leading = '0'*(4-len(str(s)))
        sid ='s'+leading+str(s)
        dat = pd.read_csv(os.path.join('..', 'target-lists',sid+'.csv'),usecols=['#TIC_ID'])
        secdict[sid] = dat.to_numpy(dtype=str)
        data = np.append(data, dat.to_numpy(dtype=str))
    ticlist = np.unique(data)
    return ticlist, secdict

def sector_search(tics, secs):
    df = pd.DataFrame({'OBJECT':tics})
    for sec in secs.keys():
        secindex = np.zeros_like(tics)
        for i in range(len(tics)):
            if tics[i] in secs[sec]:
                secindex[i] = 1
        df[sec] = secindex
    return df

year1 = np.arange(1,14)
tarlist, tarsecs = tarlist_compile(year1)
obslist = sector_search(tarlist, tarsecs)