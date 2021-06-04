# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 15:09:29 2021

@author: jlbus
"""
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt

import os

import NASA_LCs.group_tools as gt
import NASA_LCs.catalog_queries as catQ

class Group:
    def __init__(self,name,group_df,group_toi_dict=None,group_info=None):
        self.name = str(name)
        self.group_df = group_df
        self.group_toi_dict = group_toi_dict
        self.group_info = None
        self.attributes_list = []
        
    def add_tics(self,ra_col_name = 'ra',dec_col_name = 'dec', tic_col_name = None):
        if tic_col_name is not None:
            self.tics = self.group_df[tic_col_name].to_numpy(dtype = 'str')
        else:
            self.tics = catQ.get_tic_bulk(query_df = self.group_df,
                                          ra_col_name = ra_col_name, dec_col_name = dec_col_name)
        self.attributes_list.append('tics')
    def add_TIC_info(self, ra_col_name = 'ra', dec_col_name = 'dec',append_tics = True):
        self.tics,self.TIC_query = catQ.get_TIC_data_bulk(query_df = self.group_df,
                                                          ra_col_name = ra_col_name, 
                                                          dec_col_name = dec_col_name,
                                                          append_tics = append_tics)
        self.TIC_query = self.TIC_query.rename(columns = {'ID':'tic'})
        self.attributes_list.append('tics')
        self.attributes_list.append('TIC_query')
    def add_gaia_info(self, ra_col_name = 'ra', dec_col_name = 'dec', gaia_kwrgs = 'all', id_col_name = None, galactic_coords = True, delta_pm = True):
        self.gaia_query = catQ.get_gaia_data_bulk(query_df = self.group_df,
                                                 ra_col_name = ra_col_name, dec_col_name = dec_col_name,
                                                 gaia_kwrgs = gaia_kwrgs, id_col_name = id_col_name)        
        
        if galactic_coords == True:
            ##update this function to take option reference tic
            if delta_pm == True:
                self.gaia_query = catQ.add_gaia_galactic_coords(gaia_query = self.gaia_query, tic = self.group_toi_dict['tic'])
            else:
                self.gaia_query = catQ.add_gaia_galactic_coords(gaia_query = self.gaia_query, tic = None)
            
                
        self.attributes_list.append('gaia_query')    
    def add_tess_LCs(self,download_dir = None, lc_types = ['spoc','cpm'], spoc_kwrgs = ['tpf','lk_lc']):
        ## need to expand to add spoc rots later
        tic_list = self.tics
        self.rots_dict_collection = gt.bulk_download(tic_list = tic_list, 
                                                     download_dir = download_dir, 
                                                     lc_types = lc_types,
                                                     spoc_kwrgs = spoc_kwrgs)
        
        ## below could be used to update best_rots selection
        
        # if 'cpm' in lc_types:
        #     self.cpm_rot_container = {}
        #     for tic in target_dict.keys():
        #         target = target_dict[tic]
        #         if 'cpm_rot_dict' in target.available_attributes:
        #             self.cpm_rot_container[tic] = target.cpm_rot_dict
        #         else:
        #             self.cpm_rot_container[tic] = None
                    
        # if 'spoc' in lc_types:
        #     self.sap_rot_container = {}
        #     self.pdc_rot_contatiner = {}
        #     for tic in target_dict.keys():
        #         target = target_dict[tic]
        #         if 'sap_rot_dict' in target.available_attributes:
        #             self.sap_rot_container[tic] = target.sap_rot_dict
        #         else:
        #             self.sap_rot_container[tic] = None
        #         if 'pdc_rot_dict' in target.available_attributes:
        #             self.pdc_rot_container[tic] = target.pdc_rot_dict
        #         else:
        #             self.pdc_rot_container[tic] = None 
        if 'spoc' in lc_types: self.attributes_list.append('spoc_LCs')
        if 'cpm' in lc_types: self.attributes_list.append('cpm_LCs')
        
        self.attributes_list.append('rots_dict_collection')
        
        #return(target_dict)
        
    def rot_summary(self,lc_types = ['spoc','cpm'],tmag_list = None, cont_thresh = 0.7):
        ## create best rots dict (needs updating to save memory and not rely on 
        ## target_dict. possibly can rely on a rots_target_dict)
        
        ##UPDATE: fixed to rely on rots_target_dict!!!
        self.best_rots_dict = gt.best_tess_rots(rots_dict_collection = self.rots_dict_collection,
                                                lc_types = lc_types)
        
        ## add tmag summary for each lc_type within kepler
        if tmag_list is not None: 
            self.tmag_summary_dict = {}
            for lc_type in lc_types:
                temp_rot_df = self.best_rots_dict[lc_type]
                temp_rot_df = temp_rot_df.merge(right = self.TIC_query, on = 'tic', how = 'left').drop_duplicates(subset = ['tic'])
                temp_tmag_table,temp_tmag_summary = gt.add_Tmag_rot_summary(best_rots = temp_rot_df,
                                                                            cont_thresh = cont_thresh,
                                                                            tmag_list = tmag_list)
                self.tmag_summary_dict[lc_type] = {'tmag_table':temp_tmag_table,
                                                   'tmag_summary':temp_tmag_summary}
        
        self.attributes_list.append('best_rots_dict')
        self.attributes_list.append('tmag_summary_dict')
        
    #def add_final_rots():
        #condition on rots_summary
    def save_plot_df(self,fn = None,lc_type = 'cpm'):
        ## create plot df
        best_rots_df = self.best_rots_dict[lc_type]
        plot_df = self.group_df.merge(right = best_rots_df, on = 'tic', how = 'left').drop_duplicates(subset = ['tic'])
        plot_df = plot_df.merge(right = self.TIC_query, on = 'tic', how = 'left').drop_duplicates(subset = ['tic'])
        gaia_query_df = self.gaia_query.drop(columns = ['ra','dec'])
        self.plot_df = plot_df.merge(right = gaia_query_df, on = 'tic', how = 'left').drop_duplicates(subset = ['tic']).reset_index(drop = True)
        
        if fn is not None:
            self.plot_df.to_csv(fn, index = False)
        
    
    def add_pc_seq_fig(self,group_type = 'ff', lc_type = 'cpm'):#flux_type = 'cpm', color = 'bp_rp', final_rots_col = None,color_bar_kwrgs=None)
        ## create plot df
        best_rots_df = self.best_rots_dict[lc_type]
        plot_df = self.group_df.merge(right = best_rots_df, on = 'tic', how = 'left').drop_duplicates(subset = ['tic'])
        plot_df = plot_df.merge(right = self.TIC_query, on = 'tic', how = 'left').drop_duplicates(subset = ['tic'])
        gaia_query_df = self.gaia_query.drop(columns = ['ra','dec'])
        plot_df = plot_df.merge(right = gaia_query_df, on = 'tic', how = 'left').drop_duplicates(subset = ['tic'])
        
        if group_type == 'ff':
            ## just ff pc seq first
            ## NEEDS UPDATE - add one with just Tmag summary as second subplot
            ## NEEDS UPDATE - add one with big pc_seq, xyz,pmdelta, and tmag_summary
            fig_ff_pc_seq = plt.figure(figsize = (9,9))
            ax = fig_ff_pc_seq.add_subplot()
            if 'toi' in self.group_toi_dict.keys():
                title = 'TOI ' + str(self.group_toi_dict['toi']) + ' - FF Rotations'
            else:
                title = None
            gt.ff_pc_seq(ax,plot_df = plot_df, group_toi_dict = self.group_toi_dict,
                                          title = title)
            self.ff_pc_seq = fig_ff_pc_seq
            plt.close(fig_ff_pc_seq)
            
            ## now add uvwxyz plot
            tmag_summary = self.tmag_summary_dict[lc_type]
            if 'toi' in self.group_toi_dict.keys():
                toi_label = 'TOI ' + str(self.group_toi_dict['toi'])
            else:
                toi_label = None
            self.ff_uvwxyz = gt.pc_seq_uvwxyz(plot_df = plot_df, tmag_summary_dict = tmag_summary,
                                              group_toi_dict = self.group_toi_dict, toi_label = toi_label)
            
            
