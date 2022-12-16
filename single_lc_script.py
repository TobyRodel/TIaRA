import os
import numpy as np
import pandas as pd
import textwrap
from astropy.io import fits
import single_lc_funcs as fn

#Load in occurence rate tables
fgk = pd.read_csv(os.path.join('Occurrence_rates','norm_rates_FGK.csv'))
tables = {'F': fgk, 'G':fgk, 'K':fgk}

#Filter target list to include only single sector targets
target_list_path = os.path.join('target-lists', 'year1targets.csv')
target_list = pd.read_csv(target_list_path, index_col=None)
target_list['N_sec'] = target_list.drop('OBJECT', axis=1).sum(axis=1)
mask = target_list['N_sec'] == 1
target_list = target_list[mask]
del target_list['N_sec']

#Choose 1000 random rows from target lists
choices = np.random.randint(low=0, high=target_list.shape[0], size=1000)

#Loop over chosen targets, check if files exist and skip over errors
for i in choices:
    tic, sec = fn.tic_choose(target_list, i)
    path = fn.spoc_lc_path(tic, sec)[0]
    if os.path.exists(path) == True:
        try:
            data = fn.signals(path, number=10, rate_tables=tables, N_b=10, N_phase=10)
        except:
            print('Something went wrong')
            continue
        data.get_results()
        output = data.to_pandas()
        output_path = 'test_output_1000.csv'
        output.to_csv(output_path, mode='a', header=not os.path.exists(output_path))