# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 10:12:22 2021

@author: jlbus
"""
import pandas as pd
import numpy as np

import os

import NASA_LCs.catalog_queries as catQ

import NASA_LCs.lk_interface as lk_int

import NASA_LCs.rotation_tools as rot_tools

from tess_cpm.interface import cpm_interface as cpm_int

class Target:
    def __init__(self,tic = None, ra = None, dec = None):
        self.tic = tic
        self.ra = ra
        self.dec = dec
        self.available_attributes = []
        
        if self.tic is None:
            if (self.ra is not None) & (self.dec is not None):
                tic = catQ.get_tic(ra = self.ra, dec = self.dec)
            else:
                print("Please re-initialize Target object with TIC or RA/Dec.")
        if self.ra is None:
            if self.tic is not None:
                self.ra,self.dec = catQ.get_coord_from_ID(id_type = 'TIC',ID = self.tic)
            else:
                print("Please re-initialize Target object with TIC or RA/Dec.")
                
    def add_lk_LCs(self):
        self.all_LCs, self.spoc120_lc, self.lk_search_table = lk_int.get_lk_LCs(tic = self.tic)
        if len(self.spoc120_lc) > 0: self.available_attributes.append('spoc120_lc')
        if len(self.all_LCs) > 0: self.available_attributes.append('all_LCs')
        
        
    def add_cpm_LC(self,cpm_kwargs=None):
        cpm_obj = cpm_int(tic = self.tic)
        cpm_obj.download_extract()
        self.cpm_lc = cpm_obj.lc_df
        if len(self.cpm_lc) > 0: self.available_attributes.append('cpm_lc')
            
    def run_spoc_rots(self,min_freq = 1/30):
        flux_type = ['sap_flux','pdcsap_flux']
        flux_err_avail = True
        try:
            for flux in flux_type:
                LS_res,LS_periodogram_df = rot_tools.my_LS_multi_sector(lc_df = self.spoc120_lc,
                                                                          flux_type = flux,
                                                                          flux_err_avail=flux_err_avail,
                                                                          min_freq=min_freq)
                AC_res,AC_periodogram = rot_tools.exo_acf_multi_sector(lc_df = self.spoc120_lc,
                                                                            flux_type = flux,
                                                                            flux_err_avail=flux_err_avail,
                                                                            max_per = 1/min_freq)
                amp_df = rot_tools.amp_multi_sector(lc_df = self.spoc120_lc, flux_type = flux)
                rot_fig = rot_tools.period_graph(target_name = 'TIC ' + str(self.tic),
                                                      lc_df = self.spoc120_lc, flux_type = flux,
                                                      LS_res = LS_res, LS_periodogram = LS_periodogram_df,
                                                      AC_res = AC_res, AC_periodogram = AC_periodogram)
                
                if flux == 'sap_flux':
                    self.sap_rot_dict = {'LS_res':LS_res,'LS_periodogram':LS_periodogram_df,
                                         'AC_res':AC_res,'AC_periodogram':AC_periodogram,
                                         'amp_df':amp_df,'rot_fig':rot_fig}
                    self.available_attributes.append('sap_rot_dict')
                if flux == 'pdcsap_flux':
                    self.pdc_rot_dict = {'LS_res':LS_res,'LS_periodogram':LS_periodogram_df,
                                         'AC_res':AC_res,'AC_periodogram':AC_periodogram,
                                         'amp_df':amp_df,'rot_fig':rot_fig}
                    self.available_attributes.append('pdc_rot_dict')
            print("Rotations added!")
        except:
            print("Need to run 'download' function first!")
            
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
                    rot_fig = rot_tools.period_graph(target_name = 'TIC ' + str(self.tic),
                                                          lc_df = self.cpm_lc, flux_type = flux,
                                                          LS_res = LS_res, LS_periodogram = LS_periodogram_df,
                                                          AC_res = AC_res, AC_periodogram = AC_periodogram)
                    
                    if flux == 'cpm':
                        self.cpm_rot_dict = {'LS_res':LS_res,'LS_periodogram':LS_periodogram_df,
                                             'AC_res':AC_res,'AC_periodogram':AC_periodogram,
                                             'amp_df':amp_df,'rot_fig':rot_fig}
                        self.available_attributes.append('cpm_rot_dict')
                print("Rotations added!")
            except:
                self.cpm_rot_dict = {}
                print("Need to run 'download' function first!")
        else: 
            self.cpm_rot_dict = {}
            print("Need to run 'download' function first!")
