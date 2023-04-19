# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 10:12:22 2021

@author: jlbus
"""
import pandas as pd
import numpy as np

import os

import matplotlib.pyplot as plt

from astropy import units as u
from astropy.coordinates import SkyCoord

import gyros.catalog_queries as catQ
import gyros.lk_interface as lk_int
import gyros.rotation_tools as rot_tools
import gyros.target_tools as tt

from tess_cpm.interface import cpm_interface as cpm_int

class Target:
    def __init__(self,tic = None, ra = None, dec = None, 
                 query_info = True, kic = None, epic = None):
        '''
        

        Parameters
        ----------
        tic : STR or INT, optional
            Known TESS Input Catalog ID of the object. The default is None.
        ra : FLOAT, optional
            The right ascension, in degrees, of the object. The default is None.
        dec : FLOAT, optional
            The declination, in degrees, of the object. The default is None.
        query_info : BOOL, optional
            If True, query the TIC and Gaia for this object's information. The default is True.
        kic : STR or INT, optional
            Known Kepler Input Catalog ID of the object. The default is None.
        epic : STR or INT, optional
            Known Ecliptc Plane Input Catlaog ID of the object. The default is None.
        Returns
        -------
        Returns initalized Target object with TIC and Gaia Information in the 
        'target_df' attribute, with a SkyCoord attribute called 'coord'.

        '''
        self.tic = tic
        self.ra = ra
        self.dec = dec
        self.kic = kic
        self.epic = epic
        self.available_attributes = []
        
        if self.tic is None:
            if (self.ra is not None) & (self.dec is not None):
                self.tic = str(catQ.get_tic(ra = self.ra, dec = self.dec))
            else:
                print("Please re-initialize Target object with TIC or RA/Dec.")
        if self.ra is None:
            if self.tic is not None:
                self.ra,self.dec = catQ.get_coord_from_ID(id_type = 'TIC',ID = self.tic)
                self.tic = str(self.tic)
            else:
                print("Please re-initialize Target object with TIC or RA/Dec.")
                
        if query_info == True:
            self.target_info_query()
            
        self.target_name = 'TIC ' + str(self.tic)
        
                
    def target_info_query(self):
        ##add target TIC and Gaia info
        if self.tic != 'nan':
            TIC_query,_ = catQ.get_TIC_data(ra = self.ra, dec = self.dec)
            TIC_query = TIC_query.rename(columns = {'ID':'tic','ra':'RA','dec':'DEC'})
            gaia_query = catQ.get_gaia_data(ra = self.ra, dec = self.dec, gaia_kwrgs = 'all')
            gaia_query['tic'] = self.tic
            #create target df
            self.target_df = TIC_query.merge(right = gaia_query, on = 'tic', how = 'left').drop_duplicates(subset = ['tic']).reset_index(drop=True)
            self.available_attributes.append('target_df')
            
            self.sc = SkyCoord(ra = self.target_df.squeeze().ra*u.deg, 
                               dec = self.target_df.squeeze().dec*u.deg, 
                               pm_ra_cosdec=self.target_df.squeeze().pmra*u.mas/u.yr, 
                               pm_dec=self.target_df.squeeze().pmdec*u.mas/u.yr,
                               distance=u.pc*(1000./self.target_df.squeeze().parallax),
                               frame = 'icrs')
            
            self.target_df['x'] = self.sc.distance*np.cos(self.sc.galactic.b)*np.cos(self.sc.galactic.l)
            self.target_df['y'] = self.sc.distance*np.cos(self.sc.galactic.b)*np.sin(self.sc.galactic.l)
            self.target_df['z'] = self.sc.distance*np.sin(self.sc.galactic.b)
            
            ### add option for uvw if radial velocity is available from Gaia
        
        
    def add_spoc_LCs(self,tpf = False,lk_lc=True):
        # self.all_LCs, self.spoc120_lc, self.lk_search_table = lk_int.get_lk_LCs(tic = self.tic)
        # if len(self.spoc120_lc) > 0: self.available_attributes.append('spoc120_lc')
        # if len(self.all_LCs) > 0: self.available_attributes.append('all_LCs')
        if tpf == True:
            tpf_lc,pipeline_mask,threshold_mask,median_im,im_header = lk_int.tpf_sap(tic = self.tic)
            self.tpf_lc = tpf_lc
            self.tpf_lc_dict = {'pipeline_mask':pipeline_mask,
                                 'threshold_mask':threshold_mask,
                                 'median_im':median_im,
                                 'im_header':im_header}
            self.available_attributes.append('tpf_lc')
            self.available_attributes.append('tpf_lc_dict')
        if lk_lc == True:
            #note- the 'flux' in spoc120 LCs is pdcsap_flux
            spoc_lc = lk_int.spoc120(tic = self.tic)
            self.spoc_lc = spoc_lc
            self.available_attributes.append('spoc_lc')
            
    def add_kepler_LC(self):
        ### check if near kepler FOV
        kep_fov_center = SkyCoord("19h22m40s 44d30m00s", frame = 'icrs')
        sep = kep_fov_center.separation(self.sc)
        #query lightkurve for LC if on sky distance is within 12 degrees
        #Kepler FOV is 115 deg^2, so conservative check within radius of 12 degrees
        in_fov = sep < 12. *u.deg 
        
        if in_fov == True:
            #query Kepler input catalog for KIC ID
            if self.kic is None:
                self.kic = catQ.get_kic(ra = self.ra, dec = self.dec)
            
            #only query LK if KIC is available
            if str(self.kic) != 'nan':
                self.kepler_lc, self.kepler_avail = lk_int.kepler_prime_LC(kic = self.kic)
                self.available_attributes.append('kepler_lc')
                return(self.kepler_avail)
            else:
                self.kepler_lc = pd.DataFrame()
                self.kepler_avail = False
                return(self.kepler_avail)
        
        
    def add_cpm_LC(self,bkg_subtract = True, bkg_n = 300, k=100, n=100, 
                   size = 50, l2_reg = [0.1], exclusion_size = 5, apt_size = 1,
                   pred_pix_method = "similar_brightness", 
                   save_lc = False, keep_tesscut = False, 
                   add_poly = False, poly_scale = 2, poly_num_terms = 4):
        # cpm_obj = cpm_int(tic = self.tic)
        # cpm_obj.download_extract(bkg_subtract = bkg_subtract,
        #                          bkg_n = bkg_n,
        #                          k = k,
        #                          n = n,
        #                          size = size,
        #                          l2_reg = l2_reg,
        #                          exclusion_size = exclusion_size,
        #                          pred_pix_method = pred_pix_method,
        #                          save_lc = save_lc, 
        #                          keep_tesscut = keep_tesscut,
        #                          add_poly = add_poly,
        #                          poly_scale = poly_scale, 
        #                          poly_num_terms = poly_num_terms)
        # self.cpm_lc = cpm_obj.lc_df
        # self.median_cpm_im = cpm_obj.median_im
        # self.cpm_im_header = cpm_obj.im_header
        
        self.cpm_lc, self.median_cpm_im, self.cpm_im_header, self.tc_avail = lk_int.lk_tesscut(tic = self.tic,
                                                                                ra = self.ra,
                                                                                dec = self.dec,
                                                                                size = size,
                                                                                apt_size = apt_size,
                                                                                l2_reg = l2_reg,
                                                                                exclusion_size = exclusion_size,
                                                                                pred_pix_method  = pred_pix_method,
                                                                                n = n, k = k,
                                                                                bkg_subtract = bkg_subtract,
                                                                                bkg_n = bkg_n,
                                                                                add_poly = add_poly,
                                                                                poly_scale = poly_scale,
                                                                                poly_num_terms = poly_num_terms)
        if len(self.cpm_lc) > 0: 
            self.available_attributes.append('cpm_lc')
            
        return(self.tc_avail)
    
    def add_k2sff_LC(self):
        
        if self.epic is None:
            self.epic = catQ.get_epic(ra = self.ra, dec = self.dec)
            
        if str(self.epic) != 'nan':
            self.k2sff_lc, self.k2sff_avail = lk_int.k2sff_LC(epic = self.epic)
            self.available_attributes.append('k2sff_lc')
            return(self.k2sff_avail)
        else:
            self.k2sff_lc = pd.DataFrame()
            self.k2sff_avail = False
            return(self.k2sff_avail)
        
    
    def add_cpm_multi_LC(self,sectors = None, size = [32], bkg_subtract = [False], 
                         bkg_n = [40] ,k=[5], n=[35], exclusion_size = [5], apt_size = [1],
                         l2_reg = [[0.1]], pred_pix_method = ["cosine_similarity"], 
                         add_poly = False, poly_scale = 2, poly_num_terms = 4, median_im_header = False,
                         obj_num = None):
        
        self.cpm_multi_lc, self.tc_avail = lk_int.cpm_multi_lk(tic = self.tic, sectors = sectors,med_im_header = median_im_header,
                                                               size = size, bkg_subtract = bkg_subtract, 
                                                               bkg_n = bkg_n, k= k, n= n, exclusion_size = exclusion_size, 
                                                               apt_size = apt_size, l2_reg = l2_reg, 
                                                               pred_pix_method = pred_pix_method, 
                                                               add_poly = add_poly, poly_scale = poly_scale, poly_num_terms = poly_num_terms,
                                                               obj_num = obj_num)
    
    def run_cpm_multi_rots(self):
        flux_list = self.cpm_multi_lc.drop(columns = ['time','sector']).columns.to_numpy(dtype = 'str')
    
        res_df = pd.DataFrame(columns = ['extract_method','cpm_per','cpm_power',
                                           'cpm_ac_per','cpm_amp'])
    
        for flux in flux_list:
            LS_res,_ = rot_tools.my_LS(time = self.cpm_multi_lc['time'].to_numpy(),flux = self.cpm_multi_lc[flux].to_numpy())
            ac_per,_ = rot_tools.exo_acf(time = self.cpm_multi_lc['time'].to_numpy(),flux = self.cpm_multi_lc[flux].to_numpy(),flux_type = 'cpm')
            amp = rot_tools.measure_amp(self.cpm_multi_lc[flux].to_numpy())
            temp_df = pd.DataFrame(data = {'extract_method':[flux],'cpm_per':[LS_res['LS_Per1'][0]],
                                           'cpm_power':[LS_res['LS_Power1'][0]],'cpm_ac_per':[ac_per],
                                           'cpm_amp':[amp]})
            res_df = pd.concat([res_df,temp_df])
        
        
        sap_rot = self.sap_rot
        sap_amp = self.sap_amp
        
        res_df['ls_ac_div'] = np.divide(res_df['cpm_per'],res_df['cpm_ac_per'])
        res_df['cpm-ls_sap_rot_err'] = np.abs(res_df['cpm_per'] - sap_rot)/sap_rot
        res_df['cpm-ac_sap_rot_err'] = np.abs(res_df['cpm_ac_per'] - sap_rot)/sap_rot
        res_df['cpm-ls_div_sap'] = res_df['cpm_per']/sap_rot
        res_df['cpm-ac_div_sap'] = res_df['cpm_ac_per']/sap_rot
        res_df['cpm_div_sap_amp'] = res_df['cpm_amp']/sap_amp
        res_df['cpm_sap_amp_err'] = np.abs(res_df['cpm_amp'] - sap_amp)/sap_amp
        
        self.cpm_multi_rot_res = res_df
    
    def check_ffi_contamination(self,srad = 15.5*20):
        self.gaia_contam = tt.target_contam_gaia(ra = self.ra,dec = self.dec, srad = srad)
         
    
    def contamination_plot(self):
        fig = plt.figure(figsize = (8, 8))
        self.contam_fig = tt.contamination_plot(fig = fig,
                                                gaia_contam = self.gaia_contam,
                                                median_im = self.median_cpm_im,
                                                im_header = self.cpm_im_header,
                                                target_df = self.target_df)
                
    def run_cpm_rots(self,min_freq = 1/30):
        if 'cpm_lc' in self.available_attributes:
            flux_type = ['cpm']
            flux_err_avail = False
            try:
                for flux in flux_type:
                    LS_res,LS_periodogram_df = rot_tools.my_LS_multi_sector(lc_df = self.cpm_lc,
                                                                              flux_type = flux,
                                                                              flux_err_avail=flux_err_avail,
                                                                              min_freq=min_freq)
                    AC_res,AC_periodogram = rot_tools.exo_acf_multi_sector(lc_df = self.cpm_lc,
                                                                                flux_type = flux,
                                                                                flux_err_avail=flux_err_avail,
                                                                                max_per = 1/min_freq)
                    amp_df = rot_tools.amp_multi_sector(lc_df = self.cpm_lc, flux_type = flux)
                
                    if flux == 'cpm':
                        self.cpm_rot_dict = {'LS_res':LS_res,'LS_periodogram':LS_periodogram_df,
                                             'AC_res':AC_res,'AC_periodogram':AC_periodogram,
                                             'amp_df':amp_df}
                        self.available_attributes.append('cpm_rot_dict')
                print("Rotations added!")
            except:
                self.cpm_rot_dict = {}
                print("Need to run 'download' function first!")
        else: 
            self.cpm_rot_dict = {}
            print("Need to run 'download' function first!")
            
    def cpm_rot_fig(self):
        if 'cpm_rot_dict' in self.available_attributes:
            self.target_info_query()
            LS_res = self.cpm_rot_dict['LS_res']
            LS_periodogram = self.cpm_rot_dict['LS_periodogram']
            AC_res = self.cpm_rot_dict['AC_res']
            AC_periodogram = self.cpm_rot_dict['AC_periodogram']
            self.cpm_rot_fig = rot_tools.period_graph(target_name = 'TIC ' + str(self.tic),
                                          lc_df = self.cpm_lc, flux_type = 'cpm',
                                          LS_res = LS_res, LS_periodogram = LS_periodogram,
                                          AC_res = AC_res, AC_periodogram = AC_periodogram,
                                          target_df = self.target_df)
        else:
            print("Need to run_cpm_rots first!")
            
    def run_k2sff_rots(self,min_freq = 1/50):
        if 'k2sff_lc' in self.available_attributes:
            try:
                flux_type = ['flux']
                keep_cols = ['time','flux','sector']
                k2sff_lc = self.k2sff_lc.rename(columns = {'campaign':'sector'})[keep_cols]
                flux_err_avail = False
            
                for flux in flux_type:
                    LS_res,LS_periodogram_df = rot_tools.my_LS_multi_sector(lc_df = k2sff_lc,
                                                                              flux_type = flux,
                                                                              flux_err_avail=flux_err_avail,
                                                                              min_freq=min_freq)
                    AC_res,AC_periodogram = rot_tools.exo_acf_multi_sector(lc_df = k2sff_lc,
                                                                                flux_type = flux,
                                                                                flux_err_avail=flux_err_avail,
                                                                                max_per = 1/min_freq)
                    amp_df = rot_tools.amp_multi_sector(lc_df = k2sff_lc, flux_type = flux)
                
                    if flux == 'flux':
                        self.k2sff_rot_dict = {'LS_res':LS_res,'LS_periodogram':LS_periodogram_df,
                                             'AC_res':AC_res,'AC_periodogram':AC_periodogram,
                                             'amp_df':amp_df}
                        self.available_attributes.append('k2sff_rot_dict')
                print("Rotations added!")
            except:
                self.k2sff_rot_dict = {}
                print("Need to run 'download' function first!")
        else: 
            self.k2sff_rot_dict = {}
            print("Need to run 'download' function first!")
    
    def k2sff_rot_fig(self):
        if 'k2sff_rot_dict' in self.available_attributes:
            keep_cols = ['time','flux','sector']
            k2sff_lc = self.k2sff_lc.rename(columns = {'campaign':'sector'})[keep_cols]
            LS_res = self.k2sff_rot_dict['LS_res']
            LS_periodogram = self.k2sff_rot_dict['LS_periodogram']
            AC_res = self.k2sff_rot_dict['AC_res']
            AC_periodogram = self.k2sff_rot_dict['AC_periodogram']
            self.k2sff_rot_fig = rot_tools.period_graph(target_name = 'EPIC ' + str(self.epic),
                                          lc_df = k2sff_lc, flux_type = 'flux',
                                          LS_res = LS_res, LS_periodogram = LS_periodogram,
                                          AC_res = AC_res, AC_periodogram = AC_periodogram)
        else:
            print("Need to run_k2sff_rots first!")

    def run_spoc_rots(self,min_freq = 1/30):
        if 'spoc_lc' in self.available_attributes:
            flux_type = ['sap_flux','pdcsap_flux']
            flux_err_avail = True
            keep_cols = ['time','sap_flux','sap_flux_err','pdcsap_flux','pdcsap_flux_err',
                         'sector']
            spoc_lc = self.spoc_lc[keep_cols]
            try:
                for flux in flux_type:
                    LS_res,LS_periodogram_df = rot_tools.my_LS_multi_sector(lc_df = spoc_lc,
                                                                              flux_type = flux,
                                                                              flux_err_avail=flux_err_avail,
                                                                              min_freq=min_freq)
                    AC_res,AC_periodogram = rot_tools.exo_acf_multi_sector(lc_df = spoc_lc,
                                                                                flux_type = flux,
                                                                                flux_err_avail=flux_err_avail,
                                                                                max_per = 1/min_freq)
                    amp_df = rot_tools.amp_multi_sector(lc_df = spoc_lc, flux_type = flux)
                    
                    if flux == 'sap_flux':
                        self.tess_sap_rot_dict = {'LS_res':LS_res,'LS_periodogram':LS_periodogram_df,
                                              'AC_res':AC_res,'AC_periodogram':AC_periodogram,
                                              'amp_df':amp_df}
                        self.available_attributes.append('tess_sap_rot_dict')
                    if flux == 'pdcsap_flux':
                        self.tess_pdc_rot_dict = {'LS_res':LS_res,'LS_periodogram':LS_periodogram_df,
                                              'AC_res':AC_res,'AC_periodogram':AC_periodogram,
                                              'amp_df':amp_df}
                        self.available_attributes.append('tess_pdc_rot_dict')
                print("Rotations added!")
            except:
                print("Need to run 'download' function first!")
        
        if 'tpf_lc' in self.available_attributes:
            flux_type = ['flux']
            flux_err_avail = True
            try:
                for flux in flux_type:
                    LS_res,LS_periodogram_df = rot_tools.my_LS_multi_sector(lc_df = self.tpf_lc,
                                                                              flux_type = flux,
                                                                              flux_err_avail=flux_err_avail,
                                                                              min_freq=min_freq)
                    AC_res,AC_periodogram = rot_tools.exo_acf_multi_sector(lc_df = self.tpf_lc,
                                                                                flux_type = flux,
                                                                                flux_err_avail=flux_err_avail,
                                                                                max_per = 1/min_freq)
                    amp_df = rot_tools.amp_multi_sector(lc_df = self.tpf_lc, flux_type = flux)
                
                    if flux == 'flux':
                        self.tpf_rot_dict = {'LS_res':LS_res,'LS_periodogram':LS_periodogram_df,
                                             'AC_res':AC_res,'AC_periodogram':AC_periodogram,
                                             'amp_df':amp_df}
                        self.available_attributes.append('tpf_rot_dict')
                print("Rotations added!")
            except:
                self.tpf_rot_dict = {}
                print("Need to run 'download' function first!")  
                
    def run_kepler_rots(self,min_freq = 1/50):
        if 'kepler_lc' in self.available_attributes:
            flux_type = ['sap_flux','pdcsap_flux']
            flux_err_avail = True
            keep_cols = ['time','sap_flux','sap_flux_err','pdcsap_flux','pdcsap_flux_err',
                         'sector']
            kepler_lc = self.kepler_lc.rename(columns = {'quarter':'sector'})[keep_cols]
            try:
                for flux in flux_type:
                    LS_res,LS_periodogram_df = rot_tools.my_LS_multi_sector(lc_df = kepler_lc,
                                                                              flux_type = flux,
                                                                              flux_err_avail=flux_err_avail,
                                                                              min_freq=min_freq)
                    AC_res,AC_periodogram = rot_tools.exo_acf_multi_sector(lc_df = kepler_lc,
                                                                                flux_type = flux,
                                                                                flux_err_avail=flux_err_avail,
                                                                                max_per = 1/min_freq)
                    amp_df = rot_tools.amp_multi_sector(lc_df = kepler_lc, flux_type = flux)
                    
                    if flux == 'sap_flux':
                        self.kepler_sap_rot_dict = {'LS_res':LS_res,'LS_periodogram':LS_periodogram_df,
                                              'AC_res':AC_res,'AC_periodogram':AC_periodogram,
                                              'amp_df':amp_df}
                        self.available_attributes.append('kepler_sap_rot_dict')
                    if flux == 'pdcsap_flux':
                        self.kepler_pdc_rot_dict = {'LS_res':LS_res,'LS_periodogram':LS_periodogram_df,
                                              'AC_res':AC_res,'AC_periodogram':AC_periodogram,
                                              'amp_df':amp_df}
                        self.available_attributes.append('kepler_pdc_rot_dict')
                print("Rotations added!")
            except:
                print("Need to run 'download' function first!")
        
    def tpf_rot_fig(self):
        if 'tpf_rot_dict' in self.available_attributes:
            LS_res = self.tpf_rot_dict['LS_res']
            LS_periodogram = self.tpf_rot_dict['LS_periodogram']
            AC_res = self.tpf_rot_dict['AC_res']
            AC_periodogram = self.tpf_rot_dict['AC_periodogram']
            self.tpf_rot_fig = rot_tools.period_graph(target_name = 'TIC ' + str(self.tic),
                                          lc_df = self.tpf_lc, flux_type = 'flux',
                                          LS_res = LS_res, LS_periodogram = LS_periodogram,
                                          AC_res = AC_res, AC_periodogram = AC_periodogram)
        else:
            print("Need to run_spoc_rots first!")
            
    def spoc_rot_fig(self):
        flux_type = ['sap_flux','pdcsap_flux']
        self.target_info_query()
        if 'tess_pdc_rot_dict' in self.available_attributes:
            for flux in flux_type:
                if flux == 'sap_flux':
                    LS_res = self.tess_sap_rot_dict['LS_res']
                    LS_periodogram = self.tess_sap_rot_dict['LS_periodogram']
                    AC_res = self.tess_sap_rot_dict['AC_res']
                    AC_periodogram = self.tess_sap_rot_dict['AC_periodogram']
                    self.tess_sap_rot_fig = rot_tools.period_graph(target_name = 'TIC ' + str(self.tic),
                                                  lc_df = self.spoc_lc, flux_type = flux,
                                                  LS_res = LS_res, LS_periodogram = LS_periodogram,
                                                  AC_res = AC_res, AC_periodogram = AC_periodogram,
                                                  target_df = self.target_df)
                if flux == 'pdcsap_flux':
                    LS_res = self.tess_pdc_rot_dict['LS_res']
                    LS_periodogram = self.tess_pdc_rot_dict['LS_periodogram']
                    AC_res = self.tess_pdc_rot_dict['AC_res']
                    AC_periodogram = self.tess_pdc_rot_dict['AC_periodogram']
                    self.tess_pdc_rot_fig = rot_tools.period_graph(target_name = 'TIC ' + str(self.tic),
                                                  lc_df = self.spoc_lc, flux_type = flux,
                                                  LS_res = LS_res, LS_periodogram = LS_periodogram,
                                                  AC_res = AC_res, AC_periodogram = AC_periodogram,
                                                  target_df = self.target_df)
        else:
            print("Need to run_spoc_rots first!")
            
    def kepler_rot_fig(self):
        flux_type = ['sap_flux','pdcsap_flux']
        keep_cols = ['time','sap_flux','sap_flux_err','pdcsap_flux','pdcsap_flux_err',
                     'sector']
        kepler_lc = self.kepler_lc.rename(columns = {'quarter':'sector'})[keep_cols]
        self.target_info_query()
        if 'kepler_pdc_rot_dict' in self.available_attributes:
            for flux in flux_type:
                if flux == 'sap_flux':
                    LS_res = self.kepler_sap_rot_dict['LS_res']
                    LS_periodogram = self.kepler_sap_rot_dict['LS_periodogram']
                    AC_res = self.kepler_sap_rot_dict['AC_res']
                    AC_periodogram = self.kepler_sap_rot_dict['AC_periodogram']
                    self.kepler_sap_rot_fig = rot_tools.period_graph(target_name = 'KIC ' + str(self.kic),
                                                  lc_df = kepler_lc, flux_type = flux,
                                                  LS_res = LS_res, LS_periodogram = LS_periodogram,
                                                  AC_res = AC_res, AC_periodogram = AC_periodogram,
                                                  target_df = self.target_df)
                if flux == 'pdcsap_flux':
                    LS_res = self.kepler_pdc_rot_dict['LS_res']
                    LS_periodogram = self.kepler_pdc_rot_dict['LS_periodogram']
                    AC_res = self.kepler_pdc_rot_dict['AC_res']
                    AC_periodogram = self.kepler_pdc_rot_dict['AC_periodogram']
                    self.kepler_pdc_rot_fig = rot_tools.period_graph(target_name = 'KIC ' + str(self.kic),
                                                  lc_df = kepler_lc, flux_type = flux,
                                                  LS_res = LS_res, LS_periodogram = LS_periodogram,
                                                  AC_res = AC_res, AC_periodogram = AC_periodogram,
                                                  target_df = self.target_df)
        else:
            print("Need to run_kepler_rots first!")
        
        
