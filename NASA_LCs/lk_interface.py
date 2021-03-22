# -*- coding: utf-8 -*-
"""
Created on Thu Feb 25 10:35:45 2021

@author: jlbus
"""
import os 
import pandas as pd
import numpy as np

import lightkurve as lk
from lightkurve import search_lightcurve

import pickle as pkl

from astropy.table import Table
#from astropy.table import vstack


def get_lk_LCs(tic):
    #search MAST archive with tic. eventually needs to change to search_lightcurve,
    #but issues with .download_all() function from search_lightcurve
    lcf_search = search_lightcurve('tic ' + str(tic))
    
    #save search table and remove all HLSP from results 
    #(so far, only QLP HLSP seem compatitable with lightkurve)
    lcf_search_table = lcf_search.table
    lcf_search_table =lcf_search_table[lcf_search_table['obs_collection'] == 'TESS']
    
    # lcf_df = lcf_search_table.to_pandas()
    # all_lcs = []
    # for i,row in lcf_df.iterrows():
    #     temp_search_table = lcf_search_table[lcf_search_table['obs_id'] == row['obs_id']]
    #     temp_search = lk.SearchResult(table = temp_search_table)
    #     try:
    #         all_LCs.append(temp_search.download())
    #     except:
    #         print("Error lk download.")
    
    #download all lightcurves from search result
    lcf_search = lk.SearchResult(table = lcf_search_table)
    lcf = lcf_search.download_all()
    all_lcs = lcf.data #save all lightcurves to all_lcs list
    
    spoc120_lc = spoc120(lcf_search, all_lcs) #extract spoc120 LC if available
    
    return(all_lcs,spoc120_lc,lcf_search_table)
    
    # with open(lc_fn,'wb') as outfile:
    #     pkl.dump((pdc_lc_df,sap_lc_df),outfile)
        
def spoc120(lcf_search, all_lcs):
    '''
    Extracts 2-minute cadence SPOC LCs from a lightkurve query.

    Parameters
    ----------
    lcf_search : lk.search_lightcurvefile object
        DESCRIPTION.
    all_lcs : lk.search_lightcurvefile.download_all().data object
        DESCRIPTION.

    Returns
    -------
    spoc_120_lc : pandas dataframe
        SPOC 2 minute cadence for all sectors with sector label column

    '''
    search_df = lcf_search.table.to_pandas()
    spoc120_index = list(search_df[(search_df['author'] == 'SPOC') & (search_df['exptime'] == 120)].index.to_numpy(dtype = 'int'))
    if (len(spoc120_index) > 0) & (len(all_lcs) > 0):
        spoc120_LCs = [all_lcs[i] for i in spoc120_index]
    
        for i,lc in enumerate(spoc120_LCs):
            sector = lc.sector
            
            temp_lc = lc.to_pandas().reset_index()[['time','pdcsap_flux','pdcsap_flux_err','sap_flux','sap_flux_err','quality']]
            
            sector_repeats = np.repeat(a = sector, repeats = len(temp_lc))
            temp_lc['sector'] = sector_repeats
            
            if i == 0:
                spoc120_full_LC = temp_lc
            else:
                spoc120_full_LC = pd.concat([spoc120_full_LC,temp_lc])
                
        return(spoc120_full_LC)
    else:
        return(pd.DataFrame())
        