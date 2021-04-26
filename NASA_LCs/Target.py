# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 10:12:22 2021

@author: jlbus
"""
import pandas as pd
import numpy as np

import os

import matplotlib.pyplot as plt

import NASA_LCs.catalog_queries as catQ
import NASA_LCs.lk_interface as lk_int
import NASA_LCs.rotation_tools as rot_tools
import NASA_LCs.target_tools as tt

from tess_cpm.interface import cpm_interface as cpm_int

class Target:
    def __init__(self,tic = None, ra = None, dec = None):
        self.tic = tic
        self.ra = ra
        self.dec = dec
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
                
        ##add target TIC and Gaia info
        self.TIC_query,_ = catQ.get_TIC_data(ra = self.ra, dec = self.dec)
        self.TIC_query = self.TIC_query.rename(columns = {'ID':'tic','ra':'RA','dec':'DEC'})
        self.gaia_query = catQ.get_gaia_data(ra = self.ra, dec = self.dec, gaia_kwrgs = 'all')
        self.gaia_query['tic'] = self.tic
        #create target df
        self.target_df = self.TIC_query.merge(right = self.gaia_query, on = 'tic', how = 'left').drop_duplicates(subset = ['tic']).reset_index(drop=True)
                
    def add_spoc_LCs(self,tpf = True,lk_lc=True):
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
        
        
    def add_cpm_LC(self,cpm_kwargs=None):
        cpm_obj = cpm_int(tic = self.tic)
        cpm_obj.download_extract()
        self.cpm_lc = cpm_obj.lc_df
        self.median_cpm_im = cpm_obj.median_im
        self.cpm_im_header = cpm_obj.im_header
        if len(self.cpm_lc) > 0: self.available_attributes.append('cpm_lc')
        
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
            LS_res = self.cpm_rot_dict['LS_res']
            LS_periodogram = self.cpm_rot_dict['LS_periodogram']
            AC_res = self.cpm_rot_dict['AC_res']
            AC_periodogram = self.cpm_rot_dict['AC_periodogram']
            self.cpm_rot_fig = rot_tools.period_graph(target_name = 'TIC ' + str(self.tic),
                                          lc_df = self.cpm_lc, flux_type = 'cpm',
                                          LS_res = LS_res, LS_periodogram = LS_periodogram,
                                          AC_res = AC_res, AC_periodogram = AC_periodogram)
        else:
            print("Need to run_cpm_rots first!")

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
                        self.sap_rot_dict = {'LS_res':LS_res,'LS_periodogram':LS_periodogram_df,
                                              'AC_res':AC_res,'AC_periodogram':AC_periodogram,
                                              'amp_df':amp_df}
                        self.available_attributes.append('sap_rot_dict')
                    if flux == 'pdcsap_flux':
                        self.pdc_rot_dict = {'LS_res':LS_res,'LS_periodogram':LS_periodogram_df,
                                              'AC_res':AC_res,'AC_periodogram':AC_periodogram,
                                              'amp_df':amp_df}
                        self.available_attributes.append('pdc_rot_dict')
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
        else: 
            self.tpf_rot_dict = {}
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
        if 'spoc_rot_dict' in self.available_attributes:
            for flux in flux_type:
                if flux == 'sap_flux':
                    LS_res = self.sap_rot_dict['LS_res']
                    LS_periodogram = self.sap_rot_dict['LS_periodogram']
                    AC_res = self.sap_rot_dict['AC_res']
                    AC_periodogram = self.sap_rot_dict['AC_periodogram']
                    self.sap_rot_fig = rot_tools.period_graph(target_name = 'TIC ' + str(self.tic),
                                                  lc_df = self.spoc_lc, flux_type = flux,
                                                  LS_res = LS_res, LS_periodogram = LS_periodogram,
                                                  AC_res = AC_res, AC_periodogram = AC_periodogram)
                if flux == 'pdcsap_flux':
                    LS_res = self.pdc_rot_dict['LS_res']
                    LS_periodogram = self.pdc_rot_dict['LS_periodogram']
                    AC_res = self.pdc_rot_dict['AC_res']
                    AC_periodogram = self.pdc_rot_dict['AC_periodogram']
                    self.pdc_rot_fig = rot_tools.period_graph(target_name = 'TIC ' + str(self.tic),
                                                  lc_df = self.spoc_lc, flux_type = flux,
                                                  LS_res = LS_res, LS_periodogram = LS_periodogram,
                                                  AC_res = AC_res, AC_periodogram = AC_periodogram)
        else:
            print("Need to run_spoc_rots first!")
        
        
