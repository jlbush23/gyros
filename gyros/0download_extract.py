#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 14 11:12:37 2022

@author: bush
"""

#%% 0) imports

import os 
import pandas as pd
import numpy as np



import ffi_cpm as fc

pdir = r'/Users/bush/Documents/TOI1224/second_cut'
bench_fn = os.path.join(pdir, 'membership_banyan_dr3.csv')

bench = pd.read_csv(bench_fn, dtype = {'tic':np.str_})


#%%

down_dir = os.path.join(pdir,'lightcurves')
if os.path.exists(pdir) == False: os.mkdir(pdir)

test_fn = os.path.join(pdir,'TOI1224_query_success.csv')

group_name = 'TOI1224-Friends'

# xtract_param = pd.DataFrame(data = {'bkg_subtract':[True, False, False, True],
#                                     'bkg_n':[300, 0 , 0 , 300],
#                                     'k':[100, 50, 50, 50], 
#                                     'n':[100, 100, 100, 100], 
#                                     'l2_reg_cpm':[0.1, 0.1, 0.1, 1], 
#                                     'l2_reg_poly':[0.01, 0.1, 0.1, 0.1], 
#                                     'apt_size':[1, 2, 4, 9],
#                                     'exclusion_size':[5, 5, 5, 5], 
#                                     'pred_pix_method':["similar_brightness","similar_brightness","similar_brightness","similar_brightness"],
#                                     'add_poly':[False, False, False, False], 
#                                     'poly_scale':[2,2,2,2], 
#                                     'poly_num_terms':[4,4,4,4]},
#                             index = ['cpm_flux1','cpm_flux2','cpm_flux4','cpm_flux9'])

xtract_param = pd.DataFrame(data = {'bkg_subtract':[True, True],
                                    'bkg_n':[300, 300],
                                    'k':[100, 50], 
                                    'n':[100, 100], 
                                    'l2_reg_cpm':[0.1, 1], 
                                    'l2_reg_poly':[0.01, 0.1], 
                                    'apt_size':[1, 9],
                                    'exclusion_size':[5,  5], 
                                    'pred_pix_method':["similar_brightness","similar_brightness"],
                                    'add_poly':[False,  False], 
                                    'poly_scale':[2,2], 
                                    'poly_num_terms':[4,4]},
                            index = ['cpm_flux1','cpm_flux9'])

#%%
fc.ffi_cpm(query_df = bench, 
           query_fn = test_fn, 
           download_dir = down_dir, 
           group_name = group_name,
           xtract_param = xtract_param,
           keep_tesscut = False)

