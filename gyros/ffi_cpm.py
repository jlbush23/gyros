# -*- coding: utf-8 -*-
"""
Created on Mon Feb 14 11:32:36 2022

@author: jlbus
"""

#%% 0) imports

import os 
import pandas as pd
import numpy as np

import sys

from astropy.io.ascii import read

import lightkurve as lk

from astropy.io import fits
import math
import tess_cpm

import astropy.units as u
from astropy.table import Table
from astropy.table import QTable
from astropy.time import Time

import shutil

import gyros.catalog_queries as catQ
import gyros.rotation_tools as rt

from tqdm import tqdm


#%% 2) function for downloading all available FFIs 

def lk_tesscut(tic,ra = None,dec = None,size = 50,
               # bkg_subtract = True, bkg_n=300, k = 100, n = 100, 
               # l2_reg = [0.1], exclusion_size = 5, apt_size = 1, 
               # pred_pix_method = "similar_brightness", 
               # add_poly = False, poly_scale = 2, poly_num_terms = 4
               download_dir = None):
    
    if os.path.exists(download_dir) == False: os.mkdir(download_dir)
    
    #search light curve for given TIC ID
    search_res = lk.search_tesscut('TIC ' + str(tic))
    #initialize SPOC found,first found
    tesscut_found = False
    tesscut_first = False
    
    try:
        authors = search_res.table['author']
    except:
        tc_avail = False
        tesscut_lc = pd.DataFrame() # median_im = im_header = pd.DataFrame()
        return() #, median_im, im_header,tc_avail)
        
    if 'TESScut' in authors:
        tc_avail = True
    else:
        tc_avail = False
        tesscut_lc = pd.DataFrame() #median_im = im_header = pd.DataFrame()
        # return(tesscut_lc)#, median_im, im_header,tc_avail)
        return()
    
    sectors = ''
    if len(search_res.table) > 1:
        for i,sector in enumerate([sect[-2:] for sect in search_res.table['mission']]):
            if i < len(search_res.table) -1 : 
                sectors = sectors + sector + ', '
            else:
                sectors = sectors + sector
    else:
        sectors = search_res.table['mission'][0][-2:]
        
    print("Downloading TESS FFI Cutouts for Sectors " + sectors + "!")
    tesscut_collection = search_res.download_all(cutout_size = size,
                                                 download_dir = download_dir)
    
    del tesscut_collection
    
    return()
    # # lc_holder = []
    # for i in range(len(search_res)):
    #     #select search result object
    #     search_i = search_res[i]
    #     #skip if not TESScut
    #     if (search_i.author[0] == 'TESScut'): #& (search_i.exptime.data[0] == 120):
    #         print("Found " + str(search_i.mission[0]) + " TESScut data for TIC " + str(tic) + "!")
    #         tesscut_found = True
    #         if (tesscut_first == False) & (tesscut_first is not None):
    #             tesscut_first = True
    #     else:
    #         continue
    #     #download this sector's tesscut
    #     lk_tesscut_obj = search_res[i].download(cutout_size = size,
    #                                             download_dir = download_dir)

#%% 3)  single cpm lc aperture

def cpm_extraction(tesscut_path, med_im_header = False, 
              bkg_subtract = True, 
              bkg_n=300, k = 100, n = 100, l2_reg = [0.1], apt_size = 1,
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
    
    # sector = os.path.split(tesscut_path)[1].split('-')[1][-2:]
    
    #temp_cut_fn = os.path.join(path_to_tesscuts,cut)
    
    with fits.open(tesscut_path, mode="readonly") as hdu:
        x_cen = int(math.floor(hdu[1].header["1CRPX4"]))
        y_cen = int(math.floor(hdu[1].header["2CRPX4"]))
        # btjd_ref_int = hdu[0].header['BJDREFI']
        # btjd_ref_frac = hdu[0].header['BJDREFF']
        # btjd_ref = btjd_ref_int + btjd_ref_frac
        if med_im_header == True:
            #store median image
            median_im = np.nanmedian(hdu[1].data['FLUX'],axis = 0)
            im_header = hdu[2].header #used for later WCS projection
    
    temp_source = tess_cpm.Source(tesscut_path, remove_bad=True, bkg_subtract = bkg_subtract, bkg_n = bkg_n)  
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
    # print("Add CPM model")
    temp_source.add_cpm_model(exclusion_size = exclusion_size, n=n, predictor_method = pred_pix_method); 
    # temp_source.plot_cutout(l = 1, h = 99,show_aperture = True)
    # _ = temp_source.models[0][0].plot_model() #plot selected pixels 
    
    if add_poly == True:
        temp_source.add_poly_model(scale = poly_scale, num_terms = poly_num_terms); 
        temp_source.set_regs(l2_reg)# needs to be list, like [0.01,0.1] - first in list is for cpm model, second in list is for poly model          
    else:
        # print("Add L2")
        temp_source.set_regs(l2_reg) #In the current implementation the value is the reciprocal of the prior variance (i.e., precision) on the coefficients
    temp_source.holdout_fit_predict(k=k);   
    #temp_source.models[0][0].summary_plot()         
    time = Time(temp_source.time, format = 'jd')           
    flux = temp_source.get_aperture_lc(data_type="cpm_subtracted_flux")            
    # sector = np.repeat(a=sector, repeats = len(time))    
    # print("lc table")        
    lc_table = QTable([time,flux], names = ['time','cpm_flux'])            
    # lc_df = lc_table.to_pandas()  
    
    del temp_source
    
    if med_im_header == True:
        return(lc_table,median_im,im_header)
    else:
        return(lc_table)
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

#%% x) cpm to fits writer
def cpm_to_fits(tesscut_path, xtract_param, save_path = None):
    i = 0
    lc_dict = {}
    with fits.open(tesscut_path, mode="readonly") as hdu:
        og_primary_header = hdu[0].header
    for apt,row in xtract_param.iterrows():
        if row['add_poly'] == False:
            l2_reg = [row['l2_reg_cpm']]
        else:
            l2_reg = [row['l2_reg_cpm'],row['l2_reg_poly']]
        if i == 0:
            lc_tab,med_im,im_header = cpm_extraction(tesscut_path, med_im_header = True,
                                                     bkg_subtract = row['bkg_subtract'], 
                                                     bkg_n=row['bkg_n'], k = row['k'], n = row['n'], l2_reg = l2_reg, 
                                                     apt_size = row['apt_size'],
                                                     exclusion_size = row['exclusion_size'], pred_pix_method = row['pred_pix_method'], 
                                                     add_poly = row['add_poly'], poly_scale = row['poly_scale'], poly_num_terms = row['poly_num_terms'])
        else:
            lc_tab = cpm_extraction(tesscut_path, med_im_header = False,
                                                     bkg_subtract = row['bkg_subtract'], 
                                                     bkg_n=row['bkg_n'], k = row['k'], n = row['n'], l2_reg = l2_reg, 
                                                     apt_size = row['apt_size'],
                                                     exclusion_size = row['exclusion_size'], pred_pix_method = row['pred_pix_method'], 
                                                     add_poly = row['add_poly'], poly_scale = row['poly_scale'], poly_num_terms = row['poly_num_terms'])
            
        lc_dict[apt] = lc_tab
        i += 1
        del lc_tab
    
    ### LC BinTableHDU
    lc_bin_table_cols = [fits.Column(name = apt, array = lc_dict[apt]['cpm_flux'], format = 'f4') for apt,row in xtract_param.iterrows()]
    lc_bin_table_cols.insert(0, fits.Column(name = 'time', array = lc_dict[xtract_param.index[0]]['time'].value, format = 'D'))
    
    hdu_lc_table = fits.BinTableHDU.from_columns(lc_bin_table_cols, name = 'CPM_LC')
    
    ### median ImageHDU
    image_hdu = fits.ImageHDU(data = med_im, header = im_header, name = 'TESSFFI_Cutout')
    
    ### original primary HDU header
    primary_hdu = fits.PrimaryHDU(data = None, header = og_primary_header)
    
    
    ### extraction BinTableHDU
    # [fits.Column(name = col, array = xtract_param[col]) for col in xtract_param.columns]
    xtract_table = Table.from_pandas(xtract_param)
    hdu_xtract_table = fits.BinTableHDU(data = xtract_table, name = 'XTRACT_PARAMS')
    
    ### create HDU List
    hdul = fits.HDUList([primary_hdu, hdu_lc_table, hdu_xtract_table, image_hdu])
    
    if save_path is not None:
        hdul.writeto(save_path, overwrite = True)
        hdul.close()
        return()
    else:
        return(hdul)
    
#%% x) extract cpm from all available sectors, write fits, delete ffis

def cpm_all_sectors(tic_dir,
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
                                                index = ['cpm_flux1']),
                    tesscut_dir = None,
                    keep_tesscut = False):
    if tesscut_dir is None:
        tesscut_dir = os.path.join(tic_dir,'tesscut')
     
    if os.path.exists(tesscut_dir) == True:
        fits_paths = os.listdir(tesscut_dir)
    
        if len(fits_paths) == 0:
            print('No fits paths found for ' + tic_dir)
            return()
    
        else:
        ##save dir
            save_dir = os.path.join(tic_dir,'cpm')
            if os.path.exists(save_dir) == False: os.mkdir(save_dir)
            
            for path in fits_paths:
                tesscut_path = os.path.join(tesscut_dir,path)
                
                ### add in n x n cut in save_fn
                
                save_fn = os.path.split(tic_dir)[-1] + '-' + path[:10] + '_CPM.fits'
                save_path = os.path.join(save_dir,save_fn)
                
                try:
                    _ = cpm_to_fits(tesscut_path = tesscut_path, 
                                    xtract_param = xtract_param, 
                                    save_path = save_path)
                except:
                    print("Issue extracting CPM, skipping")
            
            
            if keep_tesscut == False: 
                if os.path.exists(tesscut_dir) == True: shutil.rmtree(tesscut_dir)
            
            return()
    else:
        return()
    
#%% x) update CPM fits files with rotations

def LS_fits_update(tic_dir,lc_type = 'cpm', max_per = 30):
    
    #lc_dir
    lc_dir = os.path.join(tic_dir,lc_type)
    if os.path.exists(lc_dir) == False:
        print("No light curves available for type '" + lc_type + "'.")
        return()
    ## get fits paths
    fits_paths = os.listdir(lc_dir)
    if len(fits_paths) == 0:
        print("No available FITS files.")
        return()
    
    
    for path in fits_paths:
        fits_fn = os.path.join(lc_dir,path)
        
        hdul = fits.open(fits_fn, mode = 'update')
        
        hdu_names = [hdu.name for hdu in hdul]
        
        cpm_lc = hdul[1].data
        flux_cols = [col.name for col in cpm_lc.columns if col.name != 'time']
        
        res_df_holder = []
        p_gram_dict = {}
        for col in flux_cols:
            res,p_gram = rt.my_LS(time = cpm_lc['time'],
                                  flux = cpm_lc[col],
                                  max_per = max_per)
            res.insert(loc = 0, column = 'flux ', value = col)
            res_df_holder.append(res)
            p_gram_dict[col] = p_gram
        
        ## save res, add to HDU
        res_df = pd.concat(res_df_holder)  
        ls_res_bin_table_hdu = fits.BinTableHDU(data = Table.from_pandas(res_df), name = 'LS_RES')
            
        ## update primary header with LS_Per1, LS_Per2, and powers            
        # for key in res_dict:
        #     res_df = res_dict[key]
        #     for col in res_df.columns:
        #         name = key + "_" + col
        #         hdul[0].header[name] = res_df.squeeze()[col];
        
        if 'LS_RES' in hdu_names:
            hdul['LS_RES'] = ls_res_bin_table_hdu
        else:
            hdul.append(ls_res_bin_table_hdu)
        
        ## save periodgrams to BinTable, add to HDU
            
        ls_p_gram_bin_table_cols = [fits.Column(name = key + '-power', array = p_gram_dict[key]['power'], format = 'f4') for key in p_gram_dict]
        ls_p_gram_bin_table_cols.insert(0, fits.Column(name = 'period', array = list(p_gram_dict.values())[0]['period'], format = 'f4'))
        
        ## save new fits with update option
        ls_p_gram_bintable_hdu = fits.BinTableHDU.from_columns(ls_p_gram_bin_table_cols, name = 'LS_PERIODOGRAM')
        
        if 'LS_PERIODOGRAM' in hdu_names:
            hdul['LS_PERIODOGRAM'] = ls_p_gram_bintable_hdu
        else:
            hdul.append(ls_p_gram_bintable_hdu)
    
        hdul.close()
        
        
        del p_gram_dict
        del res_df
    
    print("Successfully added rotations :)")
    return()
    
#%% x) downlaod/extract function

def ffi_cpm(query_df, query_fn, download_dir, group_name,
            size = 50,
            extract_cpm = True,
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
                                                index = ['cpm_flux1']),
            tesscut_dir = None,
            keep_tesscut = False,
            verbose = True):
    
    if os.path.exists(download_dir) == False: os.mkdir(download_dir)

    ### check for TIC IDs
    if 'tic' not in query_df.columns:
        print("Querying each object for TIC ID.")
        _ = catQ.get_tic_bulk(query_df = query_df)
        query_df.to_csv(query_fn, index = False)
    
    ### loop through query_df, download ffi, extract cpm, delete ffi
    success = []
    for i,row in tqdm(query_df.iterrows(), total = len(query_df)):
        ### print info to terminal
        print(str(i+1) + "/" + str(len(query_df)))
        
        tic = row['tic']
        
        if str(tic) == 'nan':
            if verbose == True: print("No TIC value for this target. Moving to next.")
            success.append("no_tic_match")
            continue
        
        if verbose == True: print("TIC " + str(row['tic']))
        
        tic_dir = os.path.join(download_dir,'tic' + str(tic))
        
        #### REMOVE BEFORE REUSING IN THE FUTURE
        
        ##check if already successfully extracted CPM downloaded
        if ('tic' + tic) in os.listdir(download_dir):
            tesscut_dir = os.path.join(tic_dir, 'tesscut')
            cpm_dir = os.path.join(tic_dir, 'cpm')
            
            if os.path.exists(cpm_dir) == True:
                success.append(True)
                print("Already downloaded and extracted CPM. Moving to next target.")
                continue
            elif os.path.exists(tesscut_dir) == True:
                success.append(False)
                print("TESS Cut downloaded, but some issue caused it to be unavailable.")
                ## check if TESS cut files still available. If so, run extraction. Else, redownload.
                if len(os.listdir(tesscut_dir)) > 0:
                    if extract_cpm == True:
                        _ = cpm_all_sectors(tic_dir = tic_dir,
                                            tesscut_dir = tesscut_dir,
                                            keep_tesscut = keep_tesscut,
                                            xtract_param = xtract_param)

                        _ = LS_fits_update(tic_dir = tic_dir)
                        print("Successfully extracted CPM from pre-downloaded .fits files.")
                        print("Moving to next target :)")
                        continue
                print("Deleting TESS Cut folder and trying again.")
                shutil.rmtree(tesscut_dir)
            else:
                success.append(False)
                print("TIC directory exists, but unsuccessful. Trying again.")
        
        tesscut_dir = None
        #### REMOVE ABOVE CODE BEFORE REUSING FUNCITON IN THE FUTURE
        
        ### try to download ffi
        try:
            lk_tesscut(tic = tic, download_dir = tic_dir, size = size)
            download_success = True
            success.append(download_success)
            if verbose == True: print("Success :)")
        except:
            download_success = False
            success.append(download_success)
            print("Where da hell is the TESS data??")
            if 'TESS_Sectors' in query_df.columns.to_numpy(dtype = 'str'):
                if verbose == True: print("Should be in sectors " + str(row['TESS_Sectors']))
        
        if download_success == False: 
            if verbose == True: print("Download unsuccessful, moving to next target.")
            tesscut_dir = os.path.join(tic_dir, 'tesscut')
            if os.path.exists(tesscut_dir) == True: shutil.rmtree(tesscut_dir)
            continue
        
        ### extract cpm, run rotations, delete ffi
        if extract_cpm == True:
            _ = cpm_all_sectors(tic_dir = tic_dir,
                                tesscut_dir = tesscut_dir,
                                keep_tesscut = keep_tesscut,
                                xtract_param = xtract_param)
            
            _ = LS_fits_update(tic_dir = tic_dir)
                    
        
            
        
    query_df['TESScut_avail'] = download_success
    
    query_df_dir = os.path.split(query_fn)[0]
    
    query_fn2 = os.path.join(query_df_dir,str(group_name) + '_ffi_download_success.csv')
    query_df.to_csv(query_fn2,index = False)