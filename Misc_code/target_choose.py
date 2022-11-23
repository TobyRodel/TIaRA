# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 15:32:47 2022

@author: Student
"""

import pandas as pd
import numpy as np
import os

def tic_choose(targets, row):
    '''Outputs a target id and a list of sectors from a given row of a pandas
    dataframe (when loading the csv to a dataframe set index_col=0'''
    tic_id = str(targets.iloc[row]['OBJECT'])
    sector_bool = np.array(targets.iloc[row], dtype=bool)[1:]
    sectors = np.arange(1, targets.shape[1])
    return tic_id, sectors[sector_bool]

spoc_list = pd.read_csv(os.path.join('..', 'target-lists', 'year1targets.csv')
                        , index_col=0)
length = spoc_list.shape[0]
r = np.random.randint(0, length)
target, sectors = tic_choose(spoc_list, r)
print(target)
print(sectors)