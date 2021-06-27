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
from astropy.io import fits

import math
import time

#from tess_cpm.interface import cpm_interface as cpm_int
from scipy.optimize import curve_fit

import tess_cpm

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
        
        #delete stuff
        fn = lk_lc.FILENAME
        del lk_lc
        os.remove(path = fn)
        
    if spoc_found == False:
        print("No SPOC data found for TIC " + str(tic) + ".")
        spoc_lc = pd.DataFrame()
    else:
        spoc_lc = pd.concat(lc_holder) #combine lc into 1 pandas dataframe
        
    return(spoc_lc)

def k2sff_LC(epic):
    #search light curve for given KIC ID
    search_res = lk.search_lightcurve('EPIC ' + str(epic))
    #initialize SPOC found,first found
    lc_found = False
    lc_first = False
    
    try:
        authors = search_res.table['author']
    except:
        k2sff_avail = False
        k2sff_lc = pd.DataFrame()
        return(k2sff_lc, k2sff_avail)
        
    if 'K2SFF' in authors:
        k2sff_avail = True
    else:
        k2sff_avail = False
        k2sff_lc = pd.DataFrame()
        return(k2sff_lc, k2sff_avail)
    
    lc_holder = []
    for i in range(len(search_res)):
        #select search result object
        search_i = search_res[i]
        #skip if not SPOC, 120 s exposure
        if (search_i.author.data[0] == 'K2SFF'):# & (search_i.exptime.data[0] == 120):
            print("Found " + str(search_i.mission[0]) + " data for EPIC " + str(epic) + "!")
            lc_found = True
            if (lc_first == False) & (lc_first is not None):
                lc_first = True
        else:
            continue
        lk_lc = search_res[i].download()
        lk_lc = lk_lc.remove_outliers(sigma = 5.0)
        lk_lc_df = lk_lc.to_pandas().reset_index(drop=False)
        
        lk_lc_df['campaign'] = np.repeat(a = lk_lc.CAMPAIGN, repeats = len(lk_lc)) #add sector label for my plotting functions
        lc_holder.append(lk_lc_df) #store in lc_holder
        
        #delete stuff
        fn = lk_lc.FILENAME
        del lk_lc
        os.remove(path = fn)
        
    if lc_found == False:
        print("No K2SFF data found for EPIC " + str(epic) + ".")
        k2sff_lc = pd.DataFrame()
    else:
        k2sff_lc = pd.concat(lc_holder) #combine lc into 1 pandas dataframe
        
    return(k2sff_lc,k2sff_avail)

def k2sc_LC(epic):
    #search light curve for given KIC ID
    search_res = lk.search_lightcurve('EPIC ' + str(epic))
    #initialize SPOC found,first found
    lc_found = False
    lc_first = False
    
    try:
        authors = search_res.table['author']
    except:
        k2sc_avail = False
        k2sc_lc = pd.DataFrame()
        return(k2sc_lc, k2sc_avail)
        
    if 'K2SC' in authors:
        k2sc_avail = True
    else:
        k2sc_avail = False
        k2sc_lc = pd.DataFrame()
        return(k2sc_lc, k2sc_avail)
    
    lc_holder = []
    for i in range(len(search_res)):
        #select search result object
        search_i = search_res[i]
        #skip if not SPOC, 120 s exposure
        if (search_i.author.data[0] == 'K2SC'):# & (search_i.exptime.data[0] == 120):
            print("Found " + str(search_i.mission[0]) + " data for EPIC " + str(epic) + "!")
            lc_found = True
            if (lc_first == False) & (lc_first is not None):
                lc_first = True
        else:
            continue
        lk_lc = search_res[i].download()
        lk_lc = lk_lc.remove_outliers(sigma = 5.0)
        lk_lc_df = lk_lc.to_pandas().reset_index(drop=False)
        
        lk_lc_df['campaign'] = np.repeat(a = lk_lc.CAMPAIGN, repeats = len(lk_lc)) #add sector label for my plotting functions
        lc_holder.append(lk_lc_df) #store in lc_holder
        
        #delete stuff
        fn = lk_lc.FILENAME
        del lk_lc
        os.remove(path = fn)
        
    if lc_found == False:
        print("No K2SC data found for EPIC " + str(epic) + ".")
        k2sc_lc = pd.DataFrame()
    else:
        k2sc_lc = pd.concat(lc_holder) #combine lc into 1 pandas dataframe
        
    return(k2sc_lc,k2sc_avail)

def kepler_prime_LC(kic):
    #search light curve for given KIC ID
    search_res = lk.search_lightcurve('KIC ' + str(kic))
    #initialize SPOC found,first found
    lc_found = False
    lc_first = False
    
    try:
        authors = search_res.table['author']
    except:
        kepler_avail = False
        kepler_lc = pd.DataFrame()
        return(kepler_lc, kepler_avail)
        
    if 'Kepler' in authors:
        kepler_avail = True
    else:
        kepler_avail = False
        kepler_lc = pd.DataFrame()
        return(kepler_lc, kepler_avail)
    
    lc_holder = []
    for i in range(len(search_res)):
        #select search result object
        search_i = search_res[i]
        #skip if not SPOC, 120 s exposure
        if (search_i.author.data[0] == 'Kepler'):# & (search_i.exptime.data[0] == 120):
            print("Found " + str(search_i.mission[0]) + " data for KIC " + str(kic) + "!")
            lc_found = True
            if (lc_first == False) & (lc_first is not None):
                lc_first = True
        else:
            continue
        lk_lc = search_res[i].download()
        lk_lc = lk_lc.remove_outliers(sigma = 5.0)
        lk_lc_df = lk_lc.to_pandas().reset_index(drop=False)
        
        lk_lc_df['quarter'] = np.repeat(a = lk_lc.QUARTER, repeats = len(lk_lc)) #add sector label for my plotting functions
        lc_holder.append(lk_lc_df) #store in lc_holder
        
        #delete stuff
        fn = lk_lc.FILENAME
        del lk_lc
        os.remove(path = fn)
        
    if lc_found == False:
        print("No Kepler data found for KIC " + str(kic) + ".")
        kepler_lc = pd.DataFrame()
    else:
        kepler_lc = pd.concat(lc_holder) #combine lc into 1 pandas dataframe
        
    return(kepler_lc,kepler_avail)

def lk_tesscut(tic,ra = None,dec = None,size = 32,
               bkg_subtract = False, bkg_n=40, k = 5, n = 35, l2_reg = [0.1], exclusion_size = 5, apt_size = 1, 
               pred_pix_method = "similar_brightness", add_poly = False, poly_scale = 2, poly_num_terms = 4):
    #search light curve for given TIC ID
    search_res = lk.search_tesscut('TIC ' + str(tic))
    #initialize SPOC found,first found
    tesscut_found = False
    tesscut_first = False
    
    try:
        authors = search_res.table['author']
    except:
        tc_avail = False
        tesscut_lc = median_im = im_header = pd.DataFrame()
        return(tesscut_lc, median_im, im_header,tc_avail)
        
    if 'TESScut' in authors:
        tc_avail = True
    else:
        tc_avail = False
        tesscut_lc = median_im = im_header = pd.DataFrame()
        return(tesscut_lc, median_im, im_header,tc_avail)
 
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
        
        # #instantiate cpm_obj
        # cpm_obj = cpm_int(tic = tic,ra = ra,dec = dec)
        
        # #get cpm_lc for this sector by passing lk_tess_obj to cpm_obj
        # if i == 0:
        #     med_im_header = True
        #     lc_df,median_im,im_header = cpm_obj.lk_cpm_lc(lk_tesscut_obj = lk_tesscut_obj,
        #                                                   med_im_header = med_im_header)
        # else:
        #     med_im_header = False
        #     lc_df = cpm_obj.lk_cpm_lc(lk_tesscut_obj = lk_tesscut_obj,
        #                               med_im_header = med_im_header)
        
        #get cpm_lc for this sector by passing lk_tess_obj to lk_cpm_lc function 
        if i == 0:
            med_im_header = True
            lc_df,median_im,im_header = lk_cpm_lc(lk_tesscut_obj = lk_tesscut_obj,
                                                  med_im_header = med_im_header)
        else:
            med_im_header = False
            lc_df = lk_cpm_lc(lk_tesscut_obj = lk_tesscut_obj,
                                      med_im_header = med_im_header)        
        
        #append to lc_holder for later concatenation
        lc_holder.append(lc_df) #store in lc_holder
        
        # #save median_im and im_header if i == 0
        # if i == 0:
        #     median_im = cpm_obj.median_im
        #     im_header = cpm_obj.im_header
        
        #delete stuff
        path = lk_tesscut_obj.path
        #del cpm_obj
        del lk_tesscut_obj
        os.remove(path = path)
        
    if tesscut_found == False:
        print("No TESScut data found for TIC " + str(tic) + ".")
        tesscut_lc = median_im = im_header = pd.DataFrame()
    else:
        tesscut_lc = pd.concat(lc_holder) #combine lc into 1 pandas dataframe
 
    return(tesscut_lc, median_im, im_header,tc_avail)

def cpm_multi_lk(tic, sectors = None, size = [32], bkg_subtract = [False], 
                 bkg_n = [40] ,k=[5], n=[35], exclusion_size = [5], apt_size = [1],
                 l2_reg = [[0.1]], pred_pix_method = ["cosine_similarity"], 
                 add_poly = False, poly_scale = 2, poly_num_terms = 4, med_im_header = False):
    
    #search light curve for given TIC ID
    search_res = lk.search_tesscut('TIC ' + str(tic))
    #initialize SPOC found,first found
    tesscut_found = False
    tesscut_first = False
    
    try:
        authors = search_res.table['author']
    except:
        tc_avail = False
        tesscut_lc = median_im = im_header = pd.DataFrame()
        return(multi_lc_df, tc_avail)
        
    if 'TESScut' in authors:
        tc_avail = True
    else:
        tc_avail = False
        tesscut_lc = median_im = im_header = pd.DataFrame()
        return(multi_lc_df, tc_avail)
 
    lc_holder = []
    for i in range(len(search_res)):
        #select search result object
        search_i = search_res[i]
        #skip if not TESScut
        if (search_i.author[0] == 'TESScut'): #& (search_i.exptime.data[0] == 120):
            if sectors is None:
                print("Found " + str(search_i.mission[0]) + " TESScut data for TIC " + str(tic) + "!")
                tesscut_found = True
                if (tesscut_first == False) & (tesscut_first is not None):
                    tesscut_first = True
            else:
                table_sector = search_i.table['sequence_number'][0]
                if (table_sector in np.array(sectors, dtype = 'int')) == False:
                    print("Found TESScut, but not in desired sector.")
                    continue
                else:
                    print("Found " + str(search_i.mission[0]) + " TESScut data for TIC " + str(tic) + "!")
                    tesscut_found = True
                    if (tesscut_first == False) & (tesscut_first is not None):
                        tesscut_first = True
        else:
            continue
        #download this sector's tesscut, run multi extract
        j = 0
        for sz in size:
            lk_tesscut_obj = search_res[i].download(cutout_size = sz)
            
            
            for bkg_sub in bkg_subtract:
                for bkg_N in bkg_n:
                    for N in n:
                        for exclusion in exclusion_size:
                            for apt_s in apt_size:
                                for choose_pix in pred_pix_method:
                                    for K in k:
                                        for reg in l2_reg:
                                            flux_type = choose_pix[0] + '_bkg=' + str(bkg_sub)[0]
                                            if bkg_sub: 
                                                flux_type = flux_type + '_bkgN' + str(bkg_N)
                                            else:
                                                flux_type = flux_type + '_bkgNa'
                                            flux_type = flux_type + '_s' + str(sz) + '_n' + str(N)  + '_ex' + str(exclusion) + '_apt' + str(apt_s) + '_k' + str(K) + '_l2-' + str(reg[0])
                                            print("Flux type: " + str(flux_type))
                                            if med_im_header == False:
                                                temp_lc = lk_cpm_lc(lk_tesscut_obj = lk_tesscut_obj, med_im_header = med_im_header,
                                                                      bkg_subtract = bkg_sub, bkg_n = bkg_N,
                                                                      k=K ,n=N, exclusion_size=exclusion,apt_size = apt_s,
                                                                      l2_reg = reg, pred_pix_method=choose_pix)
                                            if med_im_header == True:
                                                temp_lc,_,_ = lk_cpm_lc(lk_tesscut_obj = lk_tesscut_obj, med_im_header = med_im_header,
                                                                      bkg_subtract = bkg_sub, bkg_n = bkg_N,
                                                                      k=K ,n=N, exclusion_size=exclusion,apt_size = apt_s,
                                                                      l2_reg = reg, pred_pix_method=choose_pix)
                                    
                                            ## add cutout size to flux title!
                                            
                                            if j == 0:
                                                multi_sector_df = temp_lc.rename(columns = {'cpm':flux_type})
                                            if j > 0:
                                                multi_sector_df[flux_type] = temp_lc['cpm']
                                            j=j+1
            #delete stuff
            path = lk_tesscut_obj.path
            #del cpm_obj
            del lk_tesscut_obj
            time.sleep(5)
            os.remove(path = path)                                            
        #add multi sector df to lc holder at end of each sector extraction
        lc_holder.append(multi_sector_df)

    if tesscut_found == False:
        print("No TESScut data found for TIC " + str(tic) + "in desired sectors.")
        multi_lc_df = pd.DataFrame()
    else:
        multi_lc_df = pd.concat(lc_holder) #combine lc into 1 pandas dataframe
    
    return(multi_lc_df, tc_avail)

def lk_cpm_lc(lk_tesscut_obj, med_im_header = False, bkg_subtract = False, 
              bkg_n=40, k = 5, n = 35, l2_reg = [0.1], apt_size = 1,
              exclusion_size = 5, pred_pix_method = "similar_brightness", 
              add_poly = False, poly_scale = 2, poly_num_terms = 4):
    # if self.use_tic == True:
    #     self.cpm_lc_df_fn = "tic" + str(self.tic) + "_cpm_LC.pkl"
        
    # if self.use_tic == False:
    #     self.cpm_lc_df_fn = "ra" + str(self.ra) + "_dec" + str(self.dec) + "_cpm_LC.pkl"
        
    
    # TESS_cuts = os.listdir(path_to_tesscuts)
    # tesscut_fits = []
    # for cut in TESS_cuts: 
    #     if cut[-5:] == '.fits': tesscut_fits.append(cut)
    # cpm_lc_df_list = []
    
    # for i,cut in enumerate(tesscut_fits):
    #print(cut)
    path = lk_tesscut_obj.path
    sector = os.path.split(path)[1].split('-')[1][-2:]
    
    #temp_cut_fn = os.path.join(path_to_tesscuts,cut)
    
    with fits.open(path, mode="readonly") as hdu:
        x_cen = int(math.floor(hdu[1].header["1CRPX4"]))
        y_cen = int(math.floor(hdu[1].header["2CRPX4"]))
        if med_im_header == True:
            #store median image
            median_im = np.nanmedian(hdu[1].data['FLUX'],axis = 0)
            im_header = hdu[2].header #used for later WCS projection
    
    temp_source = tess_cpm.Source(path, remove_bad=True, bkg_subtract = bkg_subtract, bkg_n = bkg_n)  
    if apt_size == 1:          
        temp_source.set_aperture(rowlims=[y_cen,y_cen], collims=[x_cen, x_cen])  
    if apt_size == 2:
        rowlims,collims = cpm_apt_two_by_one(median_im,x_cen,y_cen)
        temp_source.set_aperture(rowlims = rowlims, collims = collims)
    if apt_size == 4:
        rowlims,collims = cpm_apt_two_by_two(median_im,x_cen,y_cen)
        temp_source.set_aperture(rowlims = rowlims, collims = collims)
    if apt_size == 9:
        temp_source.set_aperture(rowlims=[y_cen-1,y_cen+1], collims=[x_cen-1, x_cen+1])   
             
    temp_source.add_cpm_model(exclusion_size = exclusion_size, n=n, predictor_method = pred_pix_method); 
    # temp_source.plot_cutout(l = 1, h = 99,show_aperture = True)
    # _ = temp_source.models[0][0].plot_model() #plot selected pixels 
    
    if add_poly == True:
        temp_source.add_poly_model(scale = poly_scale, num_terms = poly_num_terms); 
        temp_source.set_regs(l2_reg)# needs to be list, like [0.01,0.1] - first in list is for cpm model, second in list is for poly model          
    else:
        temp_source.set_regs(l2_reg) #In the current implementation the value is the reciprocal of the prior variance (i.e., precision) on the coefficients
    temp_source.holdout_fit_predict(k=k)   
    #temp_source.models[0][0].summary_plot()         
    time = temp_source.time            
    flux = temp_source.get_aperture_lc(data_type="cpm_subtracted_flux")            
    sector = np.repeat(a=sector, repeats = len(time))            
    lc_table = Table([time,flux,sector], names = ['time','cpm','sector'])            
    lc_df = lc_table.to_pandas()  
    
    del temp_source
    
    if med_im_header == True:
        return(lc_df,median_im,im_header)
    else:
        return(lc_df)
    # if len(cpm_lc_df_list) > 0: 
    #     cpm_lc_df = pd.concat(cpm_lc_df_list)
    #     self.lc_df = cpm_lc_df
    #     return(cpm_lc_df)
    # else:
    #     return(pd.DataFrame())
    
def cpm_apt_two_by_two(median_im,x_cen,y_cen):
    center_flux = median_im[y_cen,x_cen]
    
    apt_dict = {'apt1':median_im[y_cen:y_cen+2,x_cen:x_cen+2],
                'apt2':median_im[y_cen-1:y_cen+1,x_cen:x_cen+2],
                'apt3':median_im[y_cen-1:y_cen+1,x_cen-1:x_cen+1],
                'apt4':median_im[y_cen:y_cen+2,x_cen-1:x_cen+1]}
    
    dev_holder = []
    for apt in apt_dict.keys():
        temp_dev = np.nanmean([np.abs(flux - center_flux) for flux in apt_dict[apt]])
        temp_df = pd.DataFrame(data = {'apt':[apt],'dev':[temp_dev]})
        dev_holder.append(temp_df)
    dev_df = pd.concat(dev_holder)
    best_apt = dev_df['apt'].iloc[np.where(dev_df['dev'] == np.nanmin(dev_df['dev']))][0]
    
    rows_dict = {'apt1':[y_cen,y_cen+1], 'apt2':[y_cen-1,y_cen], 'apt3':[y_cen-1,y_cen], 'apt4':[y_cen,y_cen+1]}
    cols_dict = {'apt1':[x_cen,x_cen+1], 'apt2':[x_cen,x_cen+1], 'apt3':[x_cen-1,x_cen], 'apt4':[x_cen-1,x_cen]}
    
    rowlims = rows_dict[best_apt]
    collims = cols_dict[best_apt]
    
    return(rowlims,collims)

def cpm_apt_two_by_one(median_im,x_cen,y_cen):
    center_flux = median_im[y_cen,x_cen]
    
    apt_dict = {'apt1':median_im[y_cen:y_cen+2,x_cen],
                'apt2':median_im[y_cen-1:y_cen+1,x_cen],
                'apt3':median_im[y_cen,x_cen:x_cen+2],
                'apt4':median_im[y_cen,x_cen-1:x_cen+1]}
    
    dev_holder = []
    for apt in apt_dict.keys():
        temp_dev = np.nanmean([np.abs(flux - center_flux) for flux in apt_dict[apt]])
        temp_df = pd.DataFrame(data = {'apt':[apt],'dev':[temp_dev]})
        dev_holder.append(temp_df)
    dev_df = pd.concat(dev_holder)
    best_apt = dev_df['apt'].iloc[np.where(dev_df['dev'] == np.nanmin(dev_df['dev']))][0]
    
    rows_dict = {'apt1':[y_cen,y_cen+1], 'apt2':[y_cen-1,y_cen], 'apt3':[y_cen,y_cen], 'apt4':[y_cen,y_cen]}
    cols_dict = {'apt1':[x_cen,x_cen], 'apt2':[x_cen,x_cen], 'apt3':[x_cen,x_cen+1], 'apt4':[x_cen-1,x_cen]}
    
    rowlims = rows_dict[best_apt]
    collims = cols_dict[best_apt]
    
    return(rowlims,collims)
    
    
    

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


def poly_lc_fit(x,y,y_err = None,poly_order = 3):
    if poly_order == 3:
        popt, pcov = curve_fit(poly3, x, y)
    if poly_order == 6:
        popt, pcov = curve_fit(poly6,x,y)
    
    return(popt,pcov)
    
def poly3(x,a,b,c,d):
    return(a + b*x + c*x**2 + d*x**3) 

def poly6(x,a,b,c,d,e,f,g):
    return(a + b*x + c*x**2 + d*x**3 + e*x**4 + f*x**5 + g*x**6)
        
















     