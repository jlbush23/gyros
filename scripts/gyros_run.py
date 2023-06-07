#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 14 11:12:37 2022

@author: bush
"""

#%% 0) Import modules and sample list for query.

import os 
import pandas as pd
import numpy as np

import gyros.ffi_cpm as fc

### Set your project directory path with the 'pdir' variable
# Example is set to a new directory that will be created called "test"

# The code first checks if this directory exists and creates the 
# director if it doesn't. If it does exist, that's cool, too.

pdir = r'test/'

if os.path.exists(pdir) == False: os.mkdir(pdir)

### Set the path name for the sample of stars you want lightcurves for
# Example is set to 10 stars from the MELANGE-4 Association in samples/

sample_fn = '/samples/MELANGE-4_sample.csv'

### Read in the sample as a pandas dataframe
# Code only needs the 'csv' to have RA/Dec columns
# Code first queries the TIC via MAST to find TIC ID,
# then searches for FFI data using the lightkurve package

sample_query = pd.read_csv(sample_fn)

# If you have TIC IDs for your sample already, make sure the 
# column is called 'tic' and read in that column as a string 
# data type (example below commented out)

# sample_query = pd.read_csv(sample_fn, dtype = {'tic':np.str_})


#%% 1) Inputs for the query function.

### Name a subdirectory in the project directory to store
#   the lightcurves
download_dir = os.path.join(pdir,'lightcurves')

### Name the sample you are querying 
sample_name = 'TOI1224-Friends'

### Name an output success file that reports on FFI data 
#   availability and contains all info available in the 
#   TIC for that star if a TIC match was found
output_fn = os.path.join(pdir,'TOI1224_query_success.csv')


### Set xtract parameters for the CPM extraction through the 
#   unpopular package by creating a pandas dataframe

# Option 1: Extract one lightcurve for each star with parameters
# tuned for one aperture pixel
xtract_param = pd.DataFrame(data = {'bkg_subtract':[True],
                                    'bkg_n':[300],
                                    'k':[100], 
                                    'n':[100], 
                                    'l2_reg_cpm':[0.1], 
                                    'l2_reg_poly':[0.01], 
                                    'apt_size':[1],
                                    'exclusion_size':[5], 
                                    'pred_pix_method':["similar_brightness"],
                                    'add_poly':[False], 
                                    'poly_scale':[2], 
                                    'poly_num_terms':[4]},
                            index = ['cpm_flux1'])


# Option 2: Extract 2 different lightcurves for each star with 
# two different sized apertures
# xtract_param = pd.DataFrame(data = {'bkg_subtract':[True, True],
#                                     'bkg_n':[300, 300],
#                                     'k':[100, 50], 
#                                     'n':[100, 100], 
#                                     'l2_reg_cpm':[0.1, 1], 
#                                     'l2_reg_poly':[0.01, 0.1], 
#                                     'apt_size':[1, 9],
#                                     'exclusion_size':[5,  5], 
#                                     'pred_pix_method':["similar_brightness","similar_brightness"],
#                                     'add_poly':[False,  False], 
#                                     'poly_scale':[2,2], 
#                                     'poly_num_terms':[4,4]},
#                             index = ['cpm_flux1','cpm_flux9'])

# Option 3: Extract 4 different lightcurves for each star with 
# four different sized apertures
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



#%%
fc.ffi_cpm(query_df = sample_query, 
           query_fn = sample_fn, 
           download_dir = download_dir, 
           group_name = group_name,
           xtract_param = xtract_param,
           keep_tesscut = False)

