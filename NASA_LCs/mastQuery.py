# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 15:08:18 2020

@author: jlbus
"""

import pandas as pd
import numpy as np

import os
import shutil
import pickle as pkl

from astroquery.mast import Observations
from astroquery.mast import Catalogs

from astropy import units as u
#from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
#from astropy.table import Column

import rotation_tools

class mast_target(object):
    """
    Used for downloading TESS data products from MAST.
    """
    
    def __init__(self, tic=None, ra=None, dec=None):
        """
        Takes in TIC and/or RA/Dec, download directory, and product list.
        Updates: 
            - make tic,ra,dec flexible for float/str input
            - make sure download dir is proper format
            - make sure products is "all" or a list
            - specify ResolveError and No Data Products error exceptions
        """
        
        self.tic  = tic
        self.ra = ra
        self.dec = dec
         
        if tic == None:   
            radii = np.linspace(start = 0.0001, stop = 0.001, num = 19)            
            for rad in radii:
                if self.tic == None:
                    query_string = str(self.ra) + " " + str(self.dec) # make sure to have a space between the strings!
                    obs_table = Catalogs.query_object(query_string, radius = rad*u.deg, catalog = "TIC")
                    obs_df = obs_table.to_pandas()
                    if len(obs_table['ID']) == 1:
                        self.tic = obs_table['ID'][0]
                        self.bp_rp = (obs_table['gaiabp'] - obs_table['gaiarp'])[0]
                        break
                    
                    if len(obs_df[obs_df['GAIA'].to_numpy(dtype = 'str') != '']) == 1:
                        temp_obs_df = obs_df[obs_df['GAIA'].to_numpy(dtype = 'str') != '']
                        self.tic = temp_obs_df['ID'].iloc[0]
                        self.bp_rp = (temp_obs_df['gaiabp'] - temp_obs_df['gaiarp']).iloc[0]              
                        break
                    
                    if len(np.unique(obs_df[obs_df['HIP'].to_numpy(dtype = 'str') != '']['HIP'])) == 1:
                        self.tic = obs_table['ID'][0]
                        self.bp_rp = (obs_table['gaiabp'] - obs_table['gaiarp'])[0]    
                        break
#                    
#                    if len(obs_table[obs_table['typeSrc'] == "tmgaia2"]) == 1:
#                        self.tic = obs_table['ID'][0]
#                        self.bp_rp = (obs_table['gaiabp'] - obs_table['gaiarp'])[0]
#                        break
            
            if self.tic == None:
                self.tic = "tic issue"
                #self.bp_rp = 9999
        
        if ra == None:
            query_string = "tic " + self.tic # make sure to have a space between the strings!
            obs_table = Catalogs.query_object(query_string, radius = 0.001*u.deg, catalog = "TIC")
            #obs_df = obs_table.to_pandas()
            self.ra = obs_table['ra'][0]
            self.dec = obs_table['dec'][0]
        
    def download(self,products = ['LC'], download_dir=None):
        """
        Dowload, reorganize, and simple rename desired data products.
        Updates:
            - get SPOC_df and product list ONLY option (maybe separate query() function)
            - simple rename function
        """
        self.target_path = None
        self.parent_folder = download_dir
        if self.tic != "tic issue":
            if self.tic != None:
                if download_dir != None:
                    self.folder_name = os.path.join(download_dir,self.tic) 
        
        if products == "all":
            self.products = ['LC','DVM','TP','DVT','DVS','DVR']
        else:
            self.products = products
        
        def reorganize(): #bring all files to top level
            print("Reorganizing...")
            ### EXTRACT ALL FILES
            #get down to where sector folders are, get list of folder names
            get_down = os.path.join(self.folder_name,'mastDownload','TESS')
            sub_folders = os.listdir(get_down)
            
            for sub in sub_folders:
                #get file list of subfolder
                sub_path = os.path.join(get_down,sub)
                file_list = os.listdir(sub_path)
                
                #copy files into parent self.folder_name
                for file in file_list:
                    file_path = os.path.join(sub_path,file)
                    new_loc = os.path.join(self.folder_name,file)
                    shutil.copyfile(src = file_path, dst = new_loc)
             
            ### DELETE 'mastDownload' FOLDER 
            delete_dir = os.path.join(self.folder_name,'mastDownload')
            shutil.rmtree(delete_dir)
        
        #def rename(): simple rename for all data products: tic#######_sector####_'type'.'type'
        
        def LC_extract(): #extract LC data from all sectors, label by sector
            print("Extracting LC...")
            
            # IDENTIFY LC files
            
            file_list = os.listdir(self.folder_name)            
            LC_list = []            
            for file in file_list:
                file_type = os.path.splitext(file)[-2].split("_")[-1]
                if file_type == 'lc':
                    LC_list.append(file)            
            self.LC_list = LC_list            
            # COMBINE LC'S INTO SINGLE DATAFRAME            
            temp_df_list = []                    
            for lc in LC_list:                
                # append LC data
                lc_path = os.path.join(self.folder_name,lc)                
                temp_df = self.get_LC_data(lc_path)                
                # add sector label column
                sect = lc.split("-")[1].split("s")[-1][-2:]
                temp_sector = np.repeat(a=sect,repeats = len(temp_df['time']))
                temp_df['sector'] = temp_sector
                
                temp_df_list.append(temp_df)
                
            lc_df = pd.concat(temp_df_list, ignore_index = True)
            self.spoc_lc = lc_df       
            
        ### DO THE DOWNLOADING        
        query_string = "tic " + self.tic     
        try:
            obs_table = Observations.query_object(query_string, radius = 0.0005*u.deg)        
            self.SPOC_table = obs_table[obs_table['provenance_name'] == 'SPOC']            
            self.SPOC_df = self.SPOC_table.to_pandas()     
            self.FFI_ids = self.SPOC_df[self.SPOC_df['dataproduct_type'] == 'image']['obsid'].to_numpy(dtype='str')            
            self.timeseries_ids = self.SPOC_df[self.SPOC_df['dataproduct_type'] == 'timeseries']['obsid'].to_numpy(dtype='str')            
            Observations.download_products(self.timeseries_ids, download_dir = self.folder_name, productSubGroupDescription = self.products)            
            reorganize()            
            LC_extract()           
                     
            self.query_success = "success"            
        except:
            self.query_success = "fail"
            
            
    def run_rotation_tools(self,flux_type = ['SAP_flux','PDCSAP_flux'], flux_err_avail = True, min_freq = 1/30):
        
        try:
            for flux in flux_type:
                LS_res,LS_periodogram_df = rotation_tools.my_LS_multi_sector(lc_df = self.spoc_lc,
                                                                          flux_type = flux,
                                                                          flux_err_avail=flux_err_avail,
                                                                          min_freq=min_freq)
                AC_res,AC_periodogram = rotation_tools.exo_acf_multi_sector(lc_df = self.spoc_lc,
                                                                            flux_type = flux,
                                                                            flux_err_avail=flux_err_avail,
                                                                            max_per = 1/min_freq)
                amp_df = rotation_tools.amp_multi_sector(lc_df = self.spoc_lc, flux_type = flux)
                rot_fig = rotation_tools.period_graph(target_name = 'TIC ' + str(self.tic),
                                                      lc_df = self.spoc_lc, flux_type = flux,
                                                      LS_res = LS_res, LS_periodogram = LS_periodogram_df,
                                                      AC_res = AC_res, AC_periodogram = AC_periodogram)
                
                if flux == 'SAP_flux':
                    self.LS_res_sap = LS_res
                    self.LS_periodogram_sap = LS_periodogram_df
                    self.AC_res_sap = AC_res
                    self.AC_periodogram_sap = AC_periodogram
                    self.amp_df_sap = amp_df
                    self.rot_fig_sap = rot_fig
                if flux == 'PDCSAP_flux':
                    self.LS_res_pdc = LS_res
                    self.LS_periodogram_pdc = LS_periodogram_df
                    self.AC_res_pdc = AC_res
                    self.AC_periodogram_pdc = AC_periodogram
                    self.amp_df_pdc = amp_df
                    self.rot_fig_pdc = rot_fig
            print("Rotations added!")
        except:
            print("Need to run 'download' function first!")

    ### EVENTUALLY
    
#    def info_summary(self): #create (a list ?) of query summary based on returned observation_table
#        """
#        Takes observation_table and create summary info attributes of the mastObj
#        to be fed into query_reports.
#        """

#    def check_local(self): #check to see what products have already been downloaded
#         """
#         Takes in self and checks download_dir for already downloaded products.
#         """
    
    def get_LC_data(self, lc_file, remove_nan = True):
        
        with fits.open(lc_file) as fits_file:
            #get light curve data time and flux arrays
            data = fits_file[1].data        
            time = data.TIME
            SAP_flux = data.SAP_FLUX
            SAP_flux_err = data.SAP_FLUX_ERR
            PDCSAP_flux = data.PDCSAP_FLUX
            PDCSAP_flux_err = data.PDCSAP_FLUX_ERR
            quality_flag = data.QUALITY
            
        if remove_nan == True:
            #get rid of all data points that don't have time values or don't have flux values
            time_nan = ~np.isnan(time)        
            SAP_flux_nan = ~np.isnan(SAP_flux)        
            SAP_flux_err_nan = ~np.isnan(SAP_flux)            
            PDCSAP_flux_nan = ~np.isnan(PDCSAP_flux)        
            PDCSAP_flux_err_nan = ~np.isnan(PDCSAP_flux)                
            no_nan = np.logical_and(np.logical_and(time_nan,SAP_flux_nan),SAP_flux_err_nan)            
            no_nan = np.logical_and(np.logical_and(no_nan,PDCSAP_flux_nan),PDCSAP_flux_err_nan)
        
            time = time[no_nan]
            SAP_flux = SAP_flux[no_nan]
            SAP_flux_err = SAP_flux_err[no_nan]
            PDCSAP_flux = PDCSAP_flux[no_nan]
            PDCSAP_flux_err = PDCSAP_flux_err[no_nan]
            quality_flag = quality_flag[no_nan]
            
        cols = ['time', 'SAP_flux', 'SAP_flux_err', 'PDCSAP_flux', 'PDCSAP_flux_err','quality_flag']
        dat = Table([time,SAP_flux,SAP_flux_err,PDCSAP_flux,PDCSAP_flux_err,quality_flag], names = cols)        
        df = dat.to_pandas()
        df = pd.DataFrame(df, dtype = 'float')        
        return(df)
    
    def tab2df(self,table):
        df = pd.DataFrame()        
        columns = table.colnames
        for col in columns:
            df[col] = table[col]
        return(df)  
        
        