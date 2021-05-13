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

from tess_cpm.interface import cpm_interface as cpm_int

def tpf_sap(tic):

    #search targetpixelfiles for given TIC ID
    search_res = lk.search_targetpixelfile('TIC ' + str(tic))
    #initialize SPOC found,first found
    spoc_found = False
    spoc_first = False
    
    
    lc_holder = []
    for i in range(len(search_res)):
        #select search result object
        search_i = search_res[i]
        #skip if not SPOC, 120 s exposure
        if (search_i.author.data[0] == 'SPOC') & (search_i.exptime.data[0] == 120):
            print("Found SPOC " + str(search_i.mission[0]) + " data for TIC " + str(tic) + "!")
            spoc_found = True
            if (spoc_first == False) & (spoc_first is not None):
                spoc_first = True
        else:
            continue
        tpf = search_res[i].download()
        if spoc_first == True:
            # get apertures, median image, and  WCS data for contamination plot
            pipeline_mask = tpf.pipeline_mask
            threshold_mask = tpf.create_threshold_mask(threshold = 10, reference_pixel = 'center')
            median_im = np.nanmedian(tpf.hdu[1].data['FLUX'],axis = 0)
            im_header = tpf.hdu[2].header #used for later WCS projection
            #reset spoc_first to None
            spoc_first = None
            
        # get lightcurve with pipeline mask
        lc = tpf.to_lightcurve(aperture_mask = tpf.pipeline_mask)
        lc = lc.remove_outliers(sigma = 5.0)
        lc['sector'] = np.repeat(a = lc.sector, repeats = len(lc)) #add sector label for my plotting functions
        lc_holder.append(lc.to_pandas().reset_index(drop = False)) #store in lc_holder
        
    if spoc_found == False:
        print("No SPOC data found for TIC " + str(tic) + ".")
        spoc_lc = pd.DataFrame()
        pipeline_mask = np.array([])
    else:
        spoc_lc = pd.concat(lc_holder) #combine lc into 1 pandas dataframe
        
    return(spoc_lc,pipeline_mask,threshold_mask,median_im,im_header)
        
def spoc120(tic):
    # '''
    # Extracts 2-minute cadence SPOC LCs from a lightkurve query.

    # Parameters
    # ----------
    # lcf_search : lk.search_lightcurvefile object
    #     DESCRIPTION.
    # all_lcs : lk.search_lightcurvefile.download_all().data object
    #     DESCRIPTION.

    # Returns
    # -------
    # spoc_120_lc : pandas dataframe
    #     SPOC 2 minute cadence for all sectors with sector label column

    # '''
    # search_df = lcf_search.table.to_pandas()
    # spoc120_index = list(search_df[(search_df['author'] == 'SPOC') & (search_df['exptime'] == 120)].index.to_numpy(dtype = 'int'))
    # if (len(spoc120_index) > 0) & (len(all_lcs) > 0):
    #     spoc120_LCs = [all_lcs[i] for i in spoc120_index]
    
    #     for i,lc in enumerate(spoc120_LCs):
    #         sector = lc.sector
            
    #         temp_lc = lc.to_pandas().reset_index()[['time','pdcsap_flux','pdcsap_flux_err','sap_flux','sap_flux_err','quality']]
            
    #         sector_repeats = np.repeat(a = sector, repeats = len(temp_lc))
    #         temp_lc['sector'] = sector_repeats
            
    #         if i == 0:
    #             spoc120_full_LC = temp_lc
    #         else:
    #             spoc120_full_LC = pd.concat([spoc120_full_LC,temp_lc])
                
    #     return(spoc120_full_LC)
    # else:
    #     return(pd.DataFrame())
    
    #search light curve for given TIC ID
    search_res = lk.search_lightcurve('TIC ' + str(tic))
    #initialize SPOC found,first found
    spoc_found = False
    spoc_first = False
    
    lc_holder = []
    for i in range(len(search_res)):
        #select search result object
        search_i = search_res[i]
        #skip if not SPOC, 120 s exposure
        if (search_i.author.data[0] == 'SPOC') & (search_i.exptime.data[0] == 120):
            print("Found SPOC " + str(search_i.mission[0]) + " data for TIC " + str(tic) + "!")
            spoc_found = True
            if (spoc_first == False) & (spoc_first is not None):
                spoc_first = True
        else:
            continue
        lk_lc = search_res[i].download()
        lk_lc = lk_lc.remove_outliers(sigma = 5.0)
        lk_lc_df = lk_lc.to_pandas().reset_index(drop=False)
        
        lk_lc_df['sector'] = np.repeat(a = lk_lc.sector, repeats = len(lk_lc)) #add sector label for my plotting functions
        lc_holder.append(lk_lc_df) #store in lc_holder
        
    if spoc_found == False:
        print("No SPOC data found for TIC " + str(tic) + ".")
        spoc_lc = pd.DataFrame()
    else:
        spoc_lc = pd.concat(lc_holder) #combine lc into 1 pandas dataframe
        
    return(spoc_lc)

def lk_tesscut(tic,ra = None,dec = None,size = 32):
    #search light curve for given TIC ID
    search_res = lk.search_tesscut('TIC ' + str(tic))
    #initialize SPOC found,first found
    tesscut_found = False
    tesscut_first = False
    
    lc_holder = []
    for i in range(len(search_res)):
        #select search result object
        search_i = search_res[i]
        #skip if not TESScut
        if (search_i.author[0] == 'TESScut'): #& (search_i.exptime.data[0] == 120):
            print("Found " + str(search_i.mission[0]) + " TESScut data for TIC " + str(tic) + "!")
            tesscut_found = True
            if (tesscut_first == False) & (tesscut_first is not None):
                tesscut_first = True
        else:
            continue
        #download this sector's tesscut
        lk_tesscut_obj = search_res[i].download(cutout_size = size)
        #instantiate cpm_obj
        cpm_obj = cpm_int(tic = tic,ra = ra,dec = dec)
        #get cpm_lc for this sector by passing lk_tess_obj to cpm_obj
        if i == 0:
            med_im_header = True
            lc_df,median_im,im_header = cpm_obj.lk_cpm_lc(lk_tesscut_obj = lk_tesscut_obj,
                                                          med_im_header = med_im_header)
        else:
            med_im_header = False
            lc_df = cpm_obj.lk_cpm_lc(lk_tesscut_obj = lk_tesscut_obj,
                                      med_im_header = med_im_header)
        #append to lc_holder for later concatenation
        lc_holder.append(lc_df) #store in lc_holder
        
        #save median_im and im_header if i == 0
        if i == 0:
            median_im = cpm_obj.median_im
            im_header = cpm_obj.im_header
        
        del cpm_obj
        del lk_tesscut_obj
        
    if tesscut_found == False:
        print("No TESScut data found for TIC " + str(tic) + ".")
        tesscut_lc = median_im = im_header = pd.DataFrame()
    else:
        tesscut_lc = pd.concat(lc_holder) #combine lc into 1 pandas dataframe
 
    return(tesscut_lc, median_im, im_header)
 
def get_lk_LCs(tic):
    #search MAST archive with tic. eventually needs to change to search_lightcurve,
    #but issues with .download_all() function from search_lightcurve
    lcf_search = search_lightcurve('TIC ' + str(tic))
    
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