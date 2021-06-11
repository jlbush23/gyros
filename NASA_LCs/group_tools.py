# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 21:48:07 2021

@author: jlbus
"""
# import pandas as pd
# import numpy as np

import os
import shutil
import pickle as pkl

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib as mpl

from NASA_LCs.Target import Target
from NASA_LCs.Group import Group
import NASA_LCs.catalog_queries as catQ

def gg_run(group_name,group_df,group_fn,download_dir, lc_types = ['cpm'], 
           spoc_kwrgs = None, group_toi_dict = None, target_catQ = True,use_tic = True,use_kic = False,use_epic = False):
    ## run a general group given a group df with ra,dec columns    
    #create group
    if group_toi_dict is None:
        group = Group(name = group_name, group_df = group_df)
    else:
        group = Group(name = group_name, group_df = group_df, group_toi_dict = group_toi_dict)
        
    ## determine id_list for bulk_download
    if (use_tic & use_kic) | (use_tic & use_epic) | (use_epic & use_kic):
        print("Please decide between using TIC-ID, KIC-ID, or EPIC-ID for bulk download.")
        return()
    if use_tic == True:
        if 'tic' in group_df.columns.to_numpy(dtype = 'str'):
            append_tics = False
            group.tics = group_df['tic'].to_numpy(dtype = 'str')
            group.group_df['tic'] = group.group_df['tic'].to_numpy(dtype = 'str')
        else:
            append_tics = True
            group.add_TIC_info(append_tics = append_tics)
        save_group_object(group,group_fn)
        tic_list = group.tics 
    else:
        tic_list = None
    if use_kic == True:
        if 'kic' in group_df.columns.to_numpy(dtype = 'str'):
            group.kics = group_df['kic'].to_numpy(dtype = 'str')
        else:
            group.kics = catQ.get_kic_bulk(query_df = group.group_df)
        kic_list = group.kics
    else:
        kic_list = None  
    if use_epic == True:
        if 'epic' in group_df.columns.to_numpy(dtype = 'str'):
            group.epics = group_df['epic'].to_numpy(dtype = 'str')
        else:
            gorup.epics = catQ.get_epic_bulk(query_df = group.group_df)
        epic_list = group.epics 
    else:
        epic_list = None
        
    
    #add lcs, save group with rotation_dict_collection (individual targets rotation results)
    if os.path.exists(download_dir) == False: os.mkdir(download_dir)
    lc_download_dir = os.path.join(download_dir,'lc_pickles')
    if os.path.exists(lc_download_dir) == False: os.mkdir(lc_download_dir)
    
    #add LCs by group method (issue with saving progress)
    # group.add_tess_LCs(download_dir = lc_download_dir, lc_types = lc_types, spoc_kwrgs = spoc_kwrgs)
    

    
    #add LCs externally, to update/save group object after each target
    rots_dict_collection = bulk_download(tic_list = tic_list, 
                                         kic_list = kic_list,
                                         epic_list = epic_list,
                                         download_dir = lc_download_dir, 
                                         lc_types = lc_types,
                                         spoc_kwrgs = spoc_kwrgs,
                                         group_obj = group,
                                         group_fn = group_fn)
    #group.group_df = group.group_df.merge(right = tc_avail_df, on = 'tic', how = 'left')
    group.rots_dict_collection = rots_dict_collection
    save_group_object(group,group_fn)
    
    # add Gaia query
    if target_catQ == True:
        if group_toi_dict is None:    
            group.add_gaia_info(id_col_name = 'tic',galactic_coords = True,delta_pm = False)
        else:
            group.add_gaia_info(id_col_name = 'tic',galactic_coords = True, delta_pm = True)
    save_group_object(group,group_fn)
    
    #organize best_rots, Tmag summary
    group.rot_summary(lc_types = lc_types, tmag_list = [14,15,16,99])
    save_group_object(group,group_fn)
    
    # #save plotting results
    # plot_df_fn = os.path.join(ff_product_folder,targname.replace(" ","") + "_rots.csv")
    # group.save_plot_df(fn = plot_df_fn)
    # save_group_object(group,group_fn)
    
def load_theia(groups = True, stars = True):
    lit_rot_folder = os.path.join(os.path.expanduser("~"),'NASA_LCs','literature_rotations')
    
    theia_groups_fn = os.path.join(lit_rot_folder,'theia_groups.csv')
    theia_stars_fn = os.path.join(lit_rot_folder,'theia_stars.csv')
    
    if groups: theia_groups = pd.read_csv(theia_groups_fn)
    if stars: theia_stars = pd.read_csv(theia_stars_fn)
    
    if groups & stars:
        return(theia_groups,theia_stars)
    elif (groups == True) & (stars == False):
        return(theia_groups)
    elif (groups == False) & (stars == True):
        return(theia_stars)
    else:
        print("At least one of 'groups' or 'stars' needs to be 'True'.")
        return()

def theia_group(group_num):
    theia_groups,theia_stars = load_theia()
    
    group_name = 'Theia ' + str(group_num)
    group_df = theia_stars[theia_stars['theia'] == int(group_num)].reset_index(drop = True)
    group_info = theia_groups[theia_groups['Theia'] == int(group_num)].reset_index(drop = True)
    
    group = Group(name = group_name, group_df = group_df, group_info = group_info)
    
    return(group)

def theia_run(group_num, group_fn, download_dir, lc_types = ['cpm'], group_toi_dict = None, spoc_kwrgs = None, target_catQ = True):
    #create group
    group = theia_group(group_num)
    
    # add TIC catalog info
    if 'tic' in group.group_df.columns.to_numpy(dtype = 'str'):
        append_tics = False
        group.tics = group_df['tic'].to_numpy(dtype = 'str')
        group.group_df['tic'] = group.group_df['tic'].to_numpy(dtype = 'str')
    else:
        append_tics = True
    if target_catQ == True:
        group.add_TIC_info(append_tics = append_tics)
        
    save_group_object(group,group_fn)
    
    #add lcs, save group with rotation_dict_collection (individual targets rotation results)
    lc_download_dir = os.path.join(download_dir,'lc_pickles')
    if os.path.exists(lc_download_dir) == False: os.mkdir(lc_download_dir)
    
    #add LCs by group method (issue with saving progress)
    # group.add_tess_LCs(download_dir = lc_download_dir, lc_types = lc_types, spoc_kwrgs = spoc_kwrgs)
    
    #add LCs externally, to update/save group object after each target
    rots_dict_collection = bulk_download(tic_list = group.tics, 
                                                     download_dir = lc_download_dir, 
                                                     lc_types = lc_types,
                                                     spoc_kwrgs = spoc_kwrgs,
                                                     group_obj = group,
                                                     group_fn = group_fn)
    #group.group_df = group.group_df.merge(right = tc_avail_df, on = 'tic', how = 'left')
    group.rots_dict_collection = rots_dict_collection
    save_group_object(group,group_fn)
    
    # add Gaia query
    if target_catQ == True:
        if group_toi_dict is None:    
            group.add_gaia_info(id_col_name = 'tic',galactic_coords = True,delta_pm = False)
        else:
            group.add_gaia_info(id_col_name = 'tic',galactic_coords = True, delta_pm = True)
    save_group_object(group,group_fn)
    
    #organize best_rots, Tmag summary
    group.rot_summary(lc_types = lc_types, tmag_list = [14,15,16,99], spoc_kwrgs = spoc_kwrgs)
    save_group_object(group,group_fn)
         

def load_toi_catalog(augment=True):
    lit_rot_folder = os.path.join(os.path.expanduser("~"),'NASA_LCs','literature_rotations')
    toi_cat_fn = os.path.join(lit_rot_folder,'csv-file-toi-catalog.csv')
    
    toi_df = pd.read_csv(toi_cat_fn, skiprows = 4, dtype = {'TIC':np.str_,'Full TOI ID':np.str_})
    
    if augment == True:
        toi_df = toi_df.rename(columns = {'TIC':'tic',
                                          'TIC Right Ascension':'ra', 'TIC Declination':'dec',
                                          'Tmag Value':'Tmag'})
        def strip_TOI(data):
            TOIs = data['Full TOI ID'].split(".")[0]
            return(TOIs)
        
        toi_df.insert(loc = 0, column = 'TOI', value = toi_df.apply(func = strip_TOI, axis = 1))
    
    return(toi_df)

def create_target_dict(pickle_folder):
    target_list = os.listdir(pickle_folder)
    target_dict = {}
    for i,fn in enumerate(target_list):
        print("Adding object " + str(i+1) + "/" + str(len(target_list)) + ".")
        tic = fn.split(".")[0].split("tic")[1]
        file_path = os.path.join(pickle_folder,fn)
        with open(file_path,'rb') as infile:
            targ = pkl.load(infile)
        target_dict[tic] = targ

    return(target_dict)            
    
#def create_group_pc_seq_dict():

def save_group_object(group_object,filepath):
    with open(filepath,'wb') as outfile:
        pkl.dump(group_object,outfile)

def read_group_obj(filepath):
    with open(filepath,'rb') as infile:
        group_obj = pkl.load(infile)
    return(group_obj)

def save_target_object(target_object, download_dir, id_type = None):
    #saving the file
    if id_type is not None:
        if id_type == 'tic':
            target_fn = "tic" + target_object.tic + ".pkl"            
        if id_type == 'kic':
            target_fn = "kic" + target_object.kic + ".pkl"
    else:
        try:
            target_fn = "tic" + target_object.tic + ".pkl" 
        except:
            print("Can't save this target object, no proper ID specified.")
            return()
    
    target_path = os.path.join(download_dir,target_fn) 
    with open(target_path,'wb') as outfile:
        pkl.dump(target_object,outfile)      
    
    # if keep_fits == False:
    #     #saving the file
    #     save_path = os.path.join(target_object.parent_folder,target_fn)
    #     target_object.target_path = save_path
    #     with open(save_path,'wb') as outfile:
    #         pkl.dump(target_object,outfile)
    #     shutil.rmtree(target_object.folder_name)
    # else:
    #     with open(target_path,'wb') as outfile:
    #         pkl.dump(target_object,outfile)
    
def read_target_object(filepath):
    with open(filepath,'rb') as infile:
        target_obj = pkl.load(infile)
    return(target_obj)
        
def bulk_download(tic_list, download_dir, kic_list = None, epic_list = None, 
                  lc_types = ['spoc','cpm'],spoc_kwrgs = None,
                  run_rotations = True, min_freq = 1/30,
                  #rot_options = {'flux_type':['spoc','cpm'],'flux_err_avail':[True,False],'min_freq':1/30},
                  save_objects = True, keep_fits = False, group_obj = None, group_fn = None):
    
    if tic_list is not None:
        id_list = tic_list
        id_type = 'tic'
    elif kic_list is not None:
        id_list = kic_list
        id_type = 'kic'
    elif epic_list is not None:
        id_list = epic_list
        id_type = 'epic'
    else:
        print("Need to provide either TIC or KIC list of IDs.")
        return()
    
    rots_dict_collection = {}
    tc_avail_holder = []
    kepler_avail_holder = []
    k2sff_avail_holder = []
    for i,ID in enumerate(id_list):
        print("Working on object " + str(i+1) + "/" + str(len(id_list)) + ".")
        #print(tic)
        if (str(ID) == 'nan'):
            print("No TIC for this object. Moving to next.")
            continue
        
        if tic_list is not None:
            target_obj = Target(tic = ID)
        if kic_list is not None:
            target_obj = Target(kic = ID)
        if epic_lsit is not None:
            target_obj = Target(epic = ID)
        # try spoc
        if 'spoc' in lc_types:
            #handle spoc_kwrgs
            if spoc_kwrgs is None:
                tpf_lc = lk_lc = True
            else:
                if 'tpf' in spoc_kwrgs:
                    tpf_lc = True
                else:
                    tpf_lc = False
                if 'lk_lc' in spoc_kwrgs:
                    lk_lc = True
                else:
                    lk_lc =False                   
            try: 
                target_obj.add_spoc_LCs(tpf = tpf_lc,lk_lc = lk_lc)
                if (run_rotations == True) & ('spoc_lc' in target_obj.available_attributes): target_obj.run_spoc_rots(min_freq = min_freq)
            except:
                print("Data download error. Check if SPOC available on MAST.")
            
        if 'cpm' in lc_types: 
            try: 
                tc_avail = target_obj.add_cpm_LC()
                #target_obj.check_ffi_contamination()
                tc_avail_df = pd.DataFrame(data = {'tic':[str(target_obj.tic)],'tc_avail':[tc_avail]})
                tc_avail_holder.append(tc_avail_df)
                if (run_rotations == True) & ('cpm_lc' in target_obj.available_attributes): target_obj.run_cpm_rots(min_freq = min_freq)
            except:
                print("Data download error or error with running CPM rotations.")
                
        if 'kepler' in lc_types: 
            try: 
                kepler_avail = target_obj.add_kepler_LC()
                #target_obj.check_ffi_contamination()
                kepler_avail_df = pd.DataFrame(data = {'kic':[str(target_obj.kic)],'kepler_avail':[kepler_avail]})
                kepler_avail_holder.append(kepler_avail_df)
                if (run_rotations == True) & ('kepler_lc' in target_obj.available_attributes): target_obj.run_kepler_rots(min_freq = 1/50)
            except:
                print("Data download error or error with running Kepler rotations.")
        if 'k2sff' in lc_types: 
            try: 
                k2sff_avail = target_obj.add_k2sff_LC()
                #target_obj.check_ffi_contamination()
                k2sff_avail_df = pd.DataFrame(data = {'epic':[str(target_obj.epic)],'k2sff_avail':[k2sff_avail]})
                k2sff_avail_holder.append(k2sff_avail_df)
                if (run_rotations == True) & ('k2sff_lc' in target_obj.available_attributes): target_obj.run_k2sff_rots(min_freq = 1/50)
            except:
                print("Data download error or error with running Kepler rotations.")
        
        if save_objects == True: save_target_object(target_object = target_obj, download_dir = download_dir, id_type = id_type)
        ## store rotation dicts of each object, but remove 
        target_rots_dict = {}
        if 'spoc' in lc_types:
            if 'tess_sap_rot_dict' in target_obj.available_attributes:
                sap_rot_dict = target_obj.tess_sap_rot_dict
                if 'rot_fig' in sap_rot_dict.keys():
                    del sap_rot_dict['rot_fig']
            else:
                sap_rot_dict = {}
            if 'tess_pdc_rot_dict' in target_obj.available_attributes:
                pdc_rot_dict = target_obj.tess_pdc_rot_dict
                if 'rot_fig' in pdc_rot_dict.keys():    
                    del pdc_rot_dict['rot_fig']
            else:
                pdc_rot_dict = {}
            
            target_rots_dict['tess_sap'] = sap_rot_dict
            target_rots_dict['tess_pdc'] = pdc_rot_dict

            if 'tpf_rot_dict' in target_obj.available_attributes:
                tpf_rot_dict = target_obj.tpf_rot_dict
                if 'rot_fig' in tpf_rot_dict.keys():
                    del spoc_rot_dict['rot_fig']
                
                target_rots_dict['tpf'] = tpf_rot_dict
            else:
                target_rots_dict['tpf'] = {}
            
        if 'cpm' in lc_types:
            if 'cpm_rot_dict' in target_obj.available_attributes:
                cpm_rot_dict = target_obj.cpm_rot_dict
                if 'rot_fig' in cpm_rot_dict.keys():
                    del cpm_rot_dict['rot_fig']
                cpm_rot_dict['tc_avail'] = tc_avail
                
                target_rots_dict['cpm'] = cpm_rot_dict
            else:
                target_rots_dict['cpm'] = {}
                
        if 'k2sff' in lc_types:
            if 'k2sff_rot_dict' in target_obj.available_attributes:
                k2sff_rot_dict = target_obj.k2sff_rot_dict
                if 'rot_fig' in k2sff_rot_dict.keys():
                    del k2sff_rot_dict['rot_fig']
                k2sff_rot_dict['k2sff_avail'] = k2sff_avail
                
                target_rots_dict['k2sff'] = k2sff_rot_dict
            else:
                target_rots_dict['k2sff'] = {}                
                
        if 'kepler' in lc_types:
            if 'kepler_sap_rot_dict' in target_obj.available_attributes:
                sap_rot_dict = target_obj.kepler_sap_rot_dict
                if 'rot_fig' in sap_rot_dict.keys():
                    del sap_rot_dict['rot_fig']
                sap_rot_dict['kepler_avail'] = kepler_avail
            else:
                sap_rot_dict = {}
            if 'kepler_pdc_rot_dict' in target_obj.available_attributes:
                pdc_rot_dict = target_obj.kepler_pdc_rot_dict
                if 'rot_fig' in pdc_rot_dict.keys():    
                    del pdc_rot_dict['rot_fig']
                pdc_rot_dict['kepler_avail'] = kepler_avail
            else:
                pdc_rot_dict = {}
            
            target_rots_dict['kepler_sap'] = sap_rot_dict
            target_rots_dict['kepler_pdc'] = pdc_rot_dict

        
        rots_dict_collection[id_type + str(ID)] = target_rots_dict
        if group_obj is not None: 
            group_obj.rots_dict_collection = rots_dict_collection
            if group_fn is not None: 
                save_group_object(group_obj, group_fn)
            else:
                print("Can't save group, no save filepath provided.")
        
        del target_obj        
    
    if group_obj is not None:
        if 'cpm' in lc_types: 
            tc_avail_res = pd.concat(tc_avail_holder)
            #group_obj.group_df = group_obj.group_df.merge(right = tc_avail_res, on = 'tic', how = 'left')
            group_obj.tc_avail_df = tc_avail_res
        if 'kepler' in lc_types: 
            kepler_avail_res = pd.concat(kepler_avail_holder)
            #group_obj.group_df = group_obj.group_df.merge(right = kepler_avail_res, on = 'tic', how = 'left')
            group_obj.kepler_avail_df = kepler_avail_res
        if 'k2sff' in lc_types:
            k2sff_avail_res = pd.concat(k2sff_avail_holder)
            group_obj.k2sff_avail_df = k2sff_avail_res
        if group_fn is not None: 
            save_group_object(group_obj, group_fn)
        else:
            print("Can't save group, no save filepath provided.")    
    return(rots_dict_collection)

def best_tess_rots(rots_dict_collection,lc_types = ['spoc','cpm'], spoc_kwrgs = None):
    
    if (spoc_kwrgs is None) & ('spoc' in lc_types):
        print("Need to specificy 'spoc_kwrgs' for SPOC best rotations.")
        return
    elif (spoc_kwrgs is not None) & ('spoc' in lc_types):
        if 'tpf_lc' in spoc_kwrgs: 
            tpf_lc = True
        else:
            tpf_lc = False
        if 'lk_lc' in spoc_kwrgs: 
            lk_lc = True
        else:
            lk_lc = False
    else:
        tpf_lc = lk_lc = False
        
    ## somehow feed a dict of rot dicts instead?
    
    if 'spoc' in lc_types: 
        best_tess_sap_rots = []
        best_tess_pdc_rots = []
        best_tpf_rots = []
    if 'cpm' in lc_types: best_cpm_rots =[]
    if 'k2sff' in lc_types: best_k2sff_rots = []
    if 'kepler' in lc_types:
        best_kepler_sap_rots = []
        best_kepler_pdc_rots = []
    
    for tic in rots_dict_collection.keys():
        #print(tic)
        targ_rot_dict = rots_dict_collection[tic]
        
        if ('cpm' in targ_rot_dict.keys()) & ('cpm' in lc_types):
            cpm_rot_dict = targ_rot_dict['cpm']
            
            if 'LS_res' in cpm_rot_dict.keys():
                try:
                    cpm_LS_res = cpm_rot_dict['LS_res']
                except KeyError:
                    print("Key Error. LS_res is empty, moving to next target.")
                    continue 
                
                cpm_AC_res = cpm_rot_dict['AC_res']
                
                temp_res = cpm_LS_res.merge(right = cpm_AC_res, on = 'sector', how = 'outer')
                try:
                    best_idx = np.where(np.max(temp_res['LS_Power1']) == temp_res['LS_Power1'])[0][0]
                except IndexError:
                    print("Can't index it because it has no length. No period returned.")
                    print("Moving to next target.")
                    continue 
                
                best_sector_res = temp_res.iloc[best_idx,:].to_frame().transpose()
                best_sector_res.insert(loc = 0, column = 'tic', value = str(tic))
                best_cpm_rots.append(best_sector_res)
                
        if ('tpf' in targ_rot_dict.keys()) & ('spoc' in lc_types):
            tpf_rot_dict = targ_rot_dict['tpf']
            if 'LS_res' in tpf_rot_dict.keys():
                try:
                    tpf_LS_res = tpf_rot_dict['LS_res']
                except KeyError:
                    print("Key Error. LS_res is empty, moving to next target.")
                    continue 
                
                tpf_AC_res = tpf_rot_dict['AC_res']
                
                temp_res = tpf_LS_res.merge(right = tpf_AC_res, on = 'sector', how = 'outer')
                try:
                    best_idx = np.where(np.max(temp_res['LS_Power1']) == temp_res['LS_Power1'])[0][0]
                except IndexError:
                    print("Can't index it because it has no length. No period returned.")
                    print("Moving to next target.")
                    continue 
                
                best_sector_res = temp_res.iloc[best_idx,:].to_frame().transpose()
                best_sector_res.insert(loc = 0, column = 'tic', value = str(tic))
                best_tpf_rots.append(best_sector_res)

        if ('tess_sap' in targ_rot_dict.keys()) & ('spoc' in lc_types):
            sap_rot_dict = targ_rot_dict['tess_sap']
            if 'LS_res' in sap_rot_dict.keys():
                try:
                    sap_LS_res = sap_rot_dict['LS_res']
                except KeyError:
                    print("Key Error. LS_res is empty, moving to next target.")
                    continue 
                
                sap_AC_res = sap_rot_dict['AC_res']
                
                temp_res = sap_LS_res.merge(right = sap_AC_res, on = 'sector', how = 'outer')
                try:
                    best_idx = np.where(np.max(temp_res['LS_Power1']) == temp_res['LS_Power1'])[0][0]
                except IndexError:
                    print("Can't index it because it has no length. No period returned.")
                    print("Moving to next target.")
                    continue 
                best_sector_res = temp_res.iloc[best_idx,:].to_frame().transpose()
                best_sector_res.insert(loc = 0, column = 'tic', value = str(tic))
                best_tess_sap_rots.append(best_sector_res)
            
        if ('tess_pdc' in targ_rot_dict.keys()) & ('spoc' in lc_types):
            pdc_rot_dict = targ_rot_dict['tess_pdc']
            if 'LS_res' in pdc_rot_dict.keys():
                pdc_rot_dict = targ_rot_dict['tess_pdc']
                try:
                    pdc_LS_res = pdc_rot_dict['LS_res']
                except KeyError:
                    print("Key Error. LS_res is empty, moving to next target.")
                    continue
                pdc_AC_res = pdc_rot_dict['AC_res']
                
                temp_res = pdc_LS_res.merge(right = pdc_AC_res, on = 'sector', how = 'outer')
                try:
                    best_idx = np.where(np.max(temp_res['LS_Power1']) == temp_res['LS_Power1'])[0][0]
                except IndexError:
                    print("Can't index it because it has no length. No period returned.")
                    print("Moving to next target.")
                    continue
                best_sector_res = temp_res.iloc[best_idx,:].to_frame().transpose()
                best_sector_res.insert(loc = 0, column = 'tic', value = str(tic))
                best_tess_pdc_rots.append(best_sector_res) 
                
        if ('kepler_sap' in targ_rot_dict.keys()) & ('kepler' in lc_types):
            sap_rot_dict = targ_rot_dict['kepler_sap']
            if 'LS_res' in sap_rot_dict.keys():
                try:
                    sap_LS_res = sap_rot_dict['LS_res']
                except KeyError:
                    print("Key Error. LS_res is empty, moving to next target.")
                    continue 
                
                sap_AC_res = sap_rot_dict['AC_res']
                
                temp_res = sap_LS_res.merge(right = sap_AC_res, on = 'sector', how = 'outer')
                try:
                    best_idx = np.where(np.max(temp_res['LS_Power1']) == temp_res['LS_Power1'])[0][0]
                except IndexError:
                    print("Can't index it because it has no length. No period returned.")
                    print("Moving to next target.")
                    continue 
                best_sector_res = temp_res.iloc[best_idx,:].to_frame().transpose()
                best_sector_res.insert(loc = 0, column = 'tic', value = str(tic))
                best_kepler_sap_rots.append(best_sector_res)
            
        if ('kepler_pdc' in targ_rot_dict.keys()) & ('kepler' in lc_types):
            pdc_rot_dict = targ_rot_dict['kepler_pdc']
            if 'LS_res' in pdc_rot_dict.keys():
                #pdc_rot_dict = targ_rot_dict['pdc']
                try:
                    pdc_LS_res = pdc_rot_dict['LS_res']
                except KeyError:
                    print("Key Error. LS_res is empty, moving to next target.")
                    continue
                pdc_AC_res = pdc_rot_dict['AC_res']
                
                temp_res = pdc_LS_res.merge(right = pdc_AC_res, on = 'sector', how = 'outer')
                try:
                    best_idx = np.where(np.max(temp_res['LS_Power1']) == temp_res['LS_Power1'])[0][0]
                except IndexError:
                    print("Can't index it because it has no length. No period returned.")
                    print("Moving to next target.")
                    continue
                best_sector_res = temp_res.iloc[best_idx,:].to_frame().transpose()
                best_sector_res.insert(loc = 0, column = 'tic', value = str(tic))
                best_kepler_pdc_rots.append(best_sector_res)
    
    #augment best rots df's with perc_err, matches, and aliases
    best_rots_dict = {}
    if 'cpm' in lc_types: 
        best_cpm_rots_df = pd.concat(best_cpm_rots)
        best_cpm_rots_df['perc_err'] = np.divide(best_cpm_rots_df['LS_Per1'],best_cpm_rots_df['ac_period'])
        best_cpm_rots_df['perc_err_match'] = (best_cpm_rots_df['perc_err'] < 1.1) & (best_cpm_rots_df['perc_err'] > 0.9)
        best_cpm_rots_df['alias'] = ((best_cpm_rots_df['perc_err'] < 0.55) & (best_cpm_rots_df['perc_err'] > 0.45)) | ((best_cpm_rots_df['perc_err'] < 2.05) & (best_cpm_rots_df['perc_err'] > 1.95))
        alias_type = []
        for i,row in best_cpm_rots_df.iterrows():
            if ((row['perc_err'] < 0.56) & (row['perc_err'] > 0.44)): 
                alias_type.append("0.5x")
            elif ((row['perc_err'] < 2.06) & (row['perc_err'] > 1.94)): 
                alias_type.append("2x")
            else:
                alias_type.append(np.nan)
        best_cpm_rots_df['alias_type'] = alias_type
        best_cpm_rots_df['rot_avail'] = best_cpm_rots_df['perc_err_match'] | best_cpm_rots_df['alias']
        
        best_rots_dict['cpm'] = best_cpm_rots_df
        
    if 'spoc' in lc_types: 
        if tpf_lc == True:
            best_tpf_rots_df = pd.concat(best_tpf_rots)
            best_tpf_rots_df['perc_err'] = np.divide(best_tpf_rots_df['LS_Per1'],best_tpf_rots_df['ac_period'])
            best_tpf_rots_df['perc_err_match'] = (best_tpf_rots_df['perc_err'] < 1.1) & (best_tpf_rots_df['perc_err'] > 0.9)
            best_tpf_rots_df['alias'] = ((best_tpf_rots_df['perc_err'] < 0.55) & (best_tpf_rots_df['perc_err'] > 0.45)) | ((best_tpf_rots_df['perc_err'] < 2.05) & (best_tpf_rots_df['perc_err'] > 1.95))
            alias_type = []
            for i,row in best_tpf_rots_df.iterrows():
                if ((row['perc_err'] < 0.56) & (row['perc_err'] > 0.44)): 
                    alias_type.append("0.5x")
                elif ((row['perc_err'] < 2.06) & (row['perc_err'] > 1.94)): 
                    alias_type.append("2x")
                else:
                    alias_type.append(np.nan)
            best_tpf_rots_df['alias_type'] = alias_type
            best_tpf_rots_df['rot_avail'] = best_tpf_rots_df['perc_err_match'] | best_tpf_rots_df['alias']
            
            best_rots_dict['tpf'] = best_tpf_rots_df
            
        if lk_lc == True:
            best_sap_rots_df = pd.concat(best_tess_sap_rots)
            best_sap_rots_df['perc_err'] = np.divide(best_sap_rots_df['LS_Per1'],best_sap_rots_df['ac_period'])
            best_sap_rots_df['perc_err_match'] = (best_sap_rots_df['perc_err'] < 1.1) & (best_sap_rots_df['perc_err'] > 0.9)
            best_sap_rots_df['alias'] = ((best_sap_rots_df['perc_err'] < 0.55) & (best_sap_rots_df['perc_err'] > 0.45)) | ((best_sap_rots_df['perc_err'] < 2.05) & (best_sap_rots_df['perc_err'] > 1.95))
            alias_type = []
            for i,row in best_sap_rots_df.iterrows():
                if ((row['perc_err'] < 0.56) & (row['perc_err'] > 0.44)): 
                    alias_type.append("0.5x")
                elif ((row['perc_err'] < 2.06) & (row['perc_err'] > 1.94)): 
                    alias_type.append("2x")
                else:
                    alias_type.append(np.nan)
            best_sap_rots_df['alias_type'] = alias_type
            best_sap_rots_df['rot_avail'] = best_sap_rots_df['perc_err_match'] | best_sap_rots_df['alias']
    
            best_pdc_rots_df = pd.concat(best_tess_pdc_rots)
            best_pdc_rots_df['perc_err'] = np.divide(best_pdc_rots_df['LS_Per1'],best_pdc_rots_df['ac_period'])
            best_pdc_rots_df['perc_err_match'] = (best_pdc_rots_df['perc_err'] < 1.1) & (best_pdc_rots_df['perc_err'] > 0.9)
            best_pdc_rots_df['alias'] = ((best_pdc_rots_df['perc_err'] < 0.55) & (best_pdc_rots_df['perc_err'] > 0.45)) | ((best_pdc_rots_df['perc_err'] < 2.05) & (best_pdc_rots_df['perc_err'] > 1.95))
            alias_type = []
            for i,row in best_pdc_rots_df.iterrows():
                if ((row['perc_err'] < 0.56) & (row['perc_err'] > 0.44)): 
                    alias_type.append("0.5x")
                elif ((row['perc_err'] < 2.06) & (row['perc_err'] > 1.94)): 
                    alias_type.append("2x")
                else:
                    alias_type.append(np.nan)
            best_pdc_rots_df['alias_type'] = alias_type
            best_pdc_rots_df['rot_avail'] = best_pdc_rots_df['perc_err_match'] | best_pdc_rots_df['alias']
            
            best_rots_dict['tess_sap'] = best_sap_rots_df
            best_rots_dict['tess_pdc'] = best_pdc_rots_df

    if 'kepler' in lc_types:
        if len(best_kepler_sap_rots) > 0:
            best_sap_rots_df = pd.concat(best_kepler_sap_rots)
            best_sap_rots_df['perc_err'] = np.divide(best_sap_rots_df['LS_Per1'],best_sap_rots_df['ac_period'])
            best_sap_rots_df['perc_err_match'] = (best_sap_rots_df['perc_err'] < 1.1) & (best_sap_rots_df['perc_err'] > 0.9)
            best_sap_rots_df['alias'] = ((best_sap_rots_df['perc_err'] < 0.55) & (best_sap_rots_df['perc_err'] > 0.45)) | ((best_sap_rots_df['perc_err'] < 2.05) & (best_sap_rots_df['perc_err'] > 1.95))
            alias_type = []
            for i,row in best_sap_rots_df.iterrows():
                if ((row['perc_err'] < 0.56) & (row['perc_err'] > 0.44)): 
                    alias_type.append("0.5x")
                elif ((row['perc_err'] < 2.06) & (row['perc_err'] > 1.94)): 
                    alias_type.append("2x")
                else:
                    alias_type.append(np.nan)
            best_sap_rots_df['alias_type'] = alias_type
            best_sap_rots_df['rot_avail'] = best_sap_rots_df['perc_err_match'] | best_sap_rots_df['alias']
            
            best_rots_dict['kepler_sap'] = best_sap_rots_df
        else:
            best_rots_dict['kepler_sap'] = pd.DataFrame()
        if len(best_kepler_pdc_rots) > 0:
            best_pdc_rots_df = pd.concat(best_kepler_pdc_rots)
            best_pdc_rots_df['perc_err'] = np.divide(best_pdc_rots_df['LS_Per1'],best_pdc_rots_df['ac_period'])
            best_pdc_rots_df['perc_err_match'] = (best_pdc_rots_df['perc_err'] < 1.1) & (best_pdc_rots_df['perc_err'] > 0.9)
            best_pdc_rots_df['alias'] = ((best_pdc_rots_df['perc_err'] < 0.55) & (best_pdc_rots_df['perc_err'] > 0.45)) | ((best_pdc_rots_df['perc_err'] < 2.05) & (best_pdc_rots_df['perc_err'] > 1.95))
            alias_type = []
            for i,row in best_pdc_rots_df.iterrows():
                if ((row['perc_err'] < 0.56) & (row['perc_err'] > 0.44)): 
                    alias_type.append("0.5x")
                elif ((row['perc_err'] < 2.06) & (row['perc_err'] > 1.94)): 
                    alias_type.append("2x")
                else:
                    alias_type.append(np.nan)
            best_pdc_rots_df['alias_type'] = alias_type
            best_pdc_rots_df['rot_avail'] = best_pdc_rots_df['perc_err_match'] | best_pdc_rots_df['alias']
        
            best_rots_dict['kepler_pdc'] = best_pdc_rots_df
        else:
            best_rots_dict['kepler_pdc'] = pd.DataFrame()
        
        
    return(best_rots_dict)



def add_Tmag_rot_summary(best_rots, cont_thresh = 0.7, tmag_list = [14,15,16,99]):
    best_rots = best_rots.drop_duplicates(subset = ['tic'])
    match_col_name = 'Match' #+ '(' + str(sum(best_rots['perc_err_match'] == True)) + ')'
    alias_col_name = 'Alias' #+ '(' + str(sum(best_rots['alias'] == True)) + ')'
    cont_col_name = 'ContRatio > ' + str(cont_thresh) #+ ' (' + str(sum(best_rots['contratio'] > 0.5)) + ')'
    
    row_holder = []
    for i,tmag in enumerate(tmag_list):
        if i == 0:
            best_rots_tmag = best_rots[best_rots['Tmag'] < tmag]
            row_name = 'T < ' + str(tmag) #+ ' (' + str(len(best_rots_tmag)) + ')'
        elif i == len(tmag_list) - 1:
            best_rots_tmag = best_rots[best_rots['Tmag'] > tmag_list[i-1]]
            row_name = 'T > ' + str(tmag_list[i-1]) #+ ' (' + str(len(best_rots_tmag)) + ')'
        else:
            best_rots_tmag = best_rots[(best_rots['Tmag'] > tmag_list[i-1]) & (best_rots['Tmag'] < tmag)]
            row_name = str(tmag_list[i-1]) + '<T<' + str(tmag) #+ ' (' + str(len(best_rots_tmag)) + ')'
        
        tot_tmag_avail = len(best_rots_tmag)
        #count then remove contaminated
        num_contaminated = sum(best_rots_tmag['contratio'] > cont_thresh)
        num_match_contam = sum((best_rots_tmag['contratio'] > cont_thresh) & (best_rots_tmag['perc_err_match'] == True))
        num_alias_contam = sum((best_rots_tmag['contratio'] > cont_thresh) & (best_rots_tmag['alias'] == True))
        best_rots_tmag = best_rots_tmag[np.logical_not(best_rots_tmag['contratio'] > cont_thresh)]
        num_avail = len(best_rots_tmag)
        
        #count remaining matching/alias
        num_match = sum(best_rots_tmag['perc_err_match'] == True)
        num_alias = sum(best_rots_tmag['alias'] == True)
        num_no_rot_contam = num_contaminated - num_match_contam - num_alias_contam
        #calculate return percentage
        num_return = num_match + num_alias
        perc_return = round((num_return/num_avail)*100,2)
        
        
        row_df = pd.DataFrame(data = {'Tmag':[row_name],'Total Avail':[tot_tmag_avail],
                                      match_col_name:[num_match],'Match Contam':[num_match_contam],
                                      alias_col_name:[num_alias],'Alias Contam':[num_alias_contam],'No Rot Contam':[num_no_rot_contam],
                                      cont_col_name:[num_contaminated],'Returned':[num_return],
                                      'Available':[num_avail],'Return %':[perc_return]})
        row_holder.append(row_df)
    
    tmag_table = pd.concat(row_holder)
    
    tot_match = sum(tmag_table[match_col_name])
    tot_alias = sum(tmag_table[alias_col_name])
    tot_contam = sum(tmag_table[cont_col_name])
    tot_return = sum(tmag_table['Returned'])
    tot_avail = sum(tmag_table['Available'])
    return_perc = round((tot_match + tot_alias)*100/tot_avail,2)
    
    tmag_summary_dict = {'tot_match':tot_match, 'tot_alias':tot_alias, 'tot_contam':tot_contam,
                         'tot_return':tot_return,'tot_avail':tot_avail,'return_perc':return_perc}
    
    return(tmag_table, tmag_summary_dict)

def tmag_plot(ax,tmag_df, tmag_summary, width = 0.35):
    #fig, ax = plt.subplots(figsize = (13,10))
    ax.bar(tmag_df['Tmag'],tmag_df['Total Avail'], width = width+0.15, color = 'Grey', label = 'Available')
    ax.bar(tmag_df['Tmag'],tmag_df['No Rot Contam'], width = width +0.15, color = 'Grey',
           bottom = tmag_df['Total Avail'], hatch = '///', label = 'No Rot, Contam')
    ax.bar(tmag_df['Tmag'],tmag_df['Match'],width = width, 
           color = 'mediumseagreen',label = "Match")
    ax.bar(tmag_df['Tmag'],tmag_df['Match Contam'],width = width,
           bottom = tmag_df['Match'], 
           color = 'mediumseagreen', hatch = '///', label = "Match, Contam")
    ax.bar(tmag_df['Tmag'],tmag_df['Alias'],width = width, 
           bottom = tmag_df['Match'] + tmag_df['Match Contam'], 
           color = 'mediumturquoise', label = 'Alias')
    ax.bar(tmag_df['Tmag'],tmag_df['Alias Contam'],width = width, 
           bottom = tmag_df['Match'] + tmag_df['Match Contam'] + tmag_df['Alias'], 
           color = 'mediumturquoise', hatch = '////', label = "Alias, Contam")
    
    rects = ax.patches
    heights = []
    for rect in rects:
        heights.append(rect.get_height())
        
    ax.set_ylim(ymax = 1.10*np.max(heights))
    
    for rect,label in zip(rects,tmag_df['Return %']):
            height = rect.get_height()
            ax.text(x = rect.get_x() + rect.get_width() / 2, y = 1.05*height, s = str(label) + "%",
                    ha='center', va='bottom', fontsize = 14)
    
    ax.legend(loc = 'upper right', fontsize = 'large')
    
    tot_return = tmag_summary['tot_return']
    tot_avail = tmag_summary['tot_avail']
    return_perc = tmag_summary['return_perc']
    
    ax.set_title("Hit Rate: " + str(tot_return) + " Ret / " + str(tot_avail) + " Uncontam = " + str(return_perc) + "%")# + " (Rots " + str(len(best_rots)) + " / " + str(len(group.group_df.drop_duplicates(subset = ['tic']))) + " Mems)")

def ff_pc_seq(ax,plot_df,group_toi_dict,cont_thresh = 0.7, color_type = 'bp_rp', xlim = (0.05,4), ylim = None, title = None, lit_seq = ['Pleiades','Praesepe','Hyades']):
    tic = group_toi_dict['tic']
    if 'toi' in group_toi_dict.keys():
        toi = group_toi_dict['toi']
        toi_label = 'TOI ' + str(toi)
    else:
        toi_label = 'TIC ' + str(tic)
    
    #create plot_df sub data frames for plotting different labeled data sets 
    plot_toi = plot_df[plot_df['tic'] == tic]
    plot_df2 = plot_df[np.logical_not(plot_df['contratio'] > cont_thresh)]
    match = plot_df2[plot_df2['perc_err_match'] == True]
    alias2 = plot_df2[plot_df2['alias_type'] == '2x']
    alias_half = plot_df2[plot_df2['alias_type'] == '0.5x']
    
    vmin = 0
    vmax = np.max(plot_df2['Voff(km/s)'])
    
    #add literature sequences
    ## NEEDS UPDATE: to handle different literature pc_seq function parameters
    if 'Pleiades' in lit_seq: 
        pleiades_on = True
    else: 
        pleiades_on = False
    if 'Praesepe' in lit_seq: 
        praesepe_on = True
    else: 
        praesepe_on = False        
    if 'Hyades' in lit_seq: 
        hyades_on = True
    else: 
        hyades_on = False
    if 'Upper Sco' in lit_seq: 
        upper_sco_on = True
    else: 
        upper_sco_on = False        
    pc_seq_fig(ax = ax, xlim = xlim, 
               pleiades_on = pleiades_on, praesepe_on = praesepe_on, hyades_on = hyades_on, upper_sco_on = upper_sco_on)
    mapable = ax.scatter(match[color_type],match['LS_Per1'], #c='red',
                c = match['Voff(km/s)'], cmap = "autumn", vmin = vmin, vmax = vmax, 
                s = match['LS_Power1'].to_numpy(dtype = 'float')*200, 
                alpha = 0.7, edgecolors = 'black', 
                label = 'Match', marker = 'o')
    ax.scatter(alias2[color_type],alias2['LS_Per1'], #c='red',
                c = alias2['Voff(km/s)'], cmap = "autumn", vmin = vmin, vmax = vmax, 
                s = alias2['LS_Power1'].to_numpy(dtype = 'float')*200, 
                alpha = 0.7, edgecolors = 'black',
                label = 'Alias 2x', marker = 's')
    ax.scatter(alias_half[color_type],alias_half['LS_Per1'], #c='red',
                c = alias_half['Voff(km/s)'], cmap = "autumn", vmin = vmin, vmax = vmax, 
                s = alias_half['LS_Power1'].to_numpy(dtype = 'float')*200, 
                alpha = 0.7, edgecolors = 'black',
                label = 'Alias 0.5x', marker = '^')
    ax.scatter(plot_toi[color_type],plot_toi['LS_Per1'],
                c = 'springgreen', s = 200,
                alpha = 1, edgecolors = 'black',
                label = toi_label, marker = 'x')
    ax.legend(loc = 'lower left', fontsize = 'large', framealpha = 0.4)
    if ylim is not None:
        ax.set_ylim(ylim)
    
    if title is None:
        title = 'TIC ' + str(tic) + ' - FF Rotations'
    ax.set_title(title)
    cbar = plt.colorbar(mappable = mapable, ax = ax)
    cbar.ax.get_yaxis().labelpad = 25
    cbar.ax.set_ylabel('Vtan,off (km/s)', rotation=270)
    
    return(ax)
    
    
def pc_seq_fig(ax, color_type = 'bp_rp', pleiades_on = True, praesepe_on = True, hyades_on = True, upper_sco_on = False, xlim = (0.05,3.5), guidelines_on = False):
    #This stuff is everything, use it for any python plot to make it nicer.
    import matplotlib as mpl
    mpl.rcParams['lines.linewidth'] =3
    mpl.rcParams['axes.linewidth'] = 2
    mpl.rcParams['xtick.major.width'] =2
    mpl.rcParams['ytick.major.width'] =2
    mpl.rcParams['xtick.minor.width'] =1.5
    mpl.rcParams['ytick.minor.width'] =1.5
    mpl.rcParams['ytick.labelsize'] = 17
    mpl.rcParams['xtick.labelsize'] = 17
    mpl.rcParams['axes.labelsize'] = 17
    #mpl.rcParams['legend.numpoints'] = 1
    mpl.rcParams['axes.labelweight']='semibold'
    mpl.rcParams['mathtext.fontset']='stix'
    mpl.rcParams['font.weight'] = 'semibold'
    mpl.rcParams['axes.titleweight']='semibold'
    mpl.rcParams['axes.titlesize']=17
    
    lit_rot_folder = os.path.join(os.path.expanduser("~"),'NASA_LCs','literature_rotations')
    if hyades_on == True:
        hyades_full_fn = os.path.join(lit_rot_folder,'hyades_gaia.csv')
        hyades_df = pd.read_csv(hyades_full_fn, dtype = {'tic':np.str_})
    if pleiades_on == True:
        pleiades_fn = os.path.join(lit_rot_folder,'pleiades.csv')
        pleades_full = pd.read_csv(pleiades_fn)
        pleades_single = pleades_full[pleades_full['Per2'] == 0]
        pleades_single_color = pleades_single[color_type].to_numpy(dtype = 'float')
        pleades_single_period = pleades_single['Per1'].to_numpy(dtype = 'float')
    if praesepe_on == True:
        praes_fn = os.path.join(lit_rot_folder,'praesepe.csv')
        praes_full = pd.read_csv(praes_fn)
        praes_color = praes_full[color_type].to_numpy(dtype='float')
        praes_period = praes_full['Prot'].to_numpy(dtype = 'float')
    if upper_sco_on == True:
        upper_sco_fn = os.path.join(lit_rot_folder,'upper_sco.csv')
        upper_sco_df = pd.read_csv(upper_sco_fn, dtype = {'EPIC':np.str_})
    if guidelines_on == True:
        #no NAN pleiades
        PLE_noNAN = pleades_full.dropna()
        PLE_noNAN = PLE_noNAN[(PLE_noNAN['bp_rp'] > 0.52) & (PLE_noNAN['Per1'] > 0.2) & (PLE_noNAN['bp_rp'] < 5.0)]
        PLE_noNAN = PLE_noNAN[PLE_noNAN['Per2'] == 0]
        
        color_PLE_noNAN = (PLE_noNAN['bp_rp'].dropna()).to_numpy(dtype='float')
        per_PLE_noNAN = (PLE_noNAN['Per1'].dropna()).to_numpy(dtype='float')
        
        #no NAN praesepe
        praes_columns = (praes_full.columns).to_numpy(dtype = 'str')
        praes_columns = praes_columns[:-7]
        
        PRA_noNAN = praes_full[praes_columns]
        
        PRA_noNAN = PRA_noNAN.dropna()
        
        color_PRA_noNAN = (PRA_noNAN['bp_rp'].dropna()).to_numpy(dtype='float')
        per_PRA_noNAN = (PRA_noNAN['Prot'].dropna()).to_numpy(dtype='float')
        
        ### rolling median
        from scipy.signal import medfilt
        
        PLE_plot = (pd.DataFrame(data = np.transpose(np.array([color_PLE_noNAN,per_PLE_noNAN])), columns = ['color','per'], dtype = 'float')).sort_values(by = 'color', ascending = True)
        PRA_plot = (pd.DataFrame(data = np.transpose(np.array([color_PRA_noNAN,per_PRA_noNAN])), columns = ['color','per'], dtype = 'float')).sort_values(by = 'color', ascending = True)
        
        roll_med_PLE = medfilt(PLE_plot['per'].to_numpy(dtype='float'), kernel_size = 51)
        roll_med_PRA = medfilt(PRA_plot['per'].to_numpy(dtype='float'), kernel_size = 51)
        
        def fwhm2sigma(fwhm):
            return fwhm / np.sqrt(8 * np.log(2))
        
        
        FWHM = .05
        sigma = fwhm2sigma(FWHM)
        
        smoothed_PLE_roll = np.zeros(roll_med_PLE.shape)
        for i,x_position in enumerate(PLE_plot['color'].to_numpy(dtype='float')):
            kernel = np.exp(-(PLE_plot['color'].to_numpy(dtype='float') - x_position) ** 2 / (2 * sigma ** 2))
            kernel = kernel / sum(kernel)
            smoothed_PLE_roll[i] = sum(roll_med_PLE * kernel)
        
        smoothed_PRA_roll = np.zeros(roll_med_PRA.shape)
        for i,x_position in enumerate(PRA_plot['color'].to_numpy(dtype='float')):
            kernel = np.exp(-(PRA_plot['color'].to_numpy(dtype='float') - x_position) ** 2 / (2 * sigma ** 2))
            kernel = kernel / sum(kernel)
            smoothed_PRA_roll[i] = sum(roll_med_PRA * kernel)
    
    ####### YES and Qmarks Only
    #fig = plt.figure(figsize = (7,7))
    
    if pleiades_on == True:
        ax.scatter(pleades_single_color, pleades_single_period, c = 'deepskyblue',
                    edgecolors = 'black', s=50,marker = 'o',alpha = 0.35, 
                    label = 'Pleiades - 100 Myr')
    if praesepe_on == True:
        ax.scatter(praes_color, praes_period, c = 'purple',s=50, marker = 's', 
                    alpha = 0.25, label = 'Praesepe - 625 Myr')
    if guidelines_on == True:
        import matplotlib.patheffects as pe
        ax.plot(PLE_plot['color'],smoothed_PLE_roll, '--',c = 'deepskyblue', linewidth = 3, path_effects=[pe.Stroke(linewidth=6, foreground='k'), pe.Normal()])
        ax.plot(PRA_plot['color'],smoothed_PRA_roll, '--', c = 'purple', linewidth = 3, path_effects=[pe.Stroke(linewidth=6, foreground='k'), pe.Normal()])
    if hyades_on == True:
        ax.scatter(hyades_df['bp_rp'],hyades_df['rotation'], c = 'darkgreen', 
                    marker = '^', s = 50, alpha = 0.85, edgecolors = 'black', 
                    label = 'Hyades - 625 Myr')
    if upper_sco_on == True:
        ax.scatter(upper_sco_df['bp_rp'],upper_sco_df['P1'], c = 'grey',
                    marker = '^', s = 50, alpha = 0.8, edgecolors = 'black',
                    label = 'Upper Scorpius - 15 Myr')
    
    ax.set_xlim(xlim)
    ax.set_yscale('log')
    ax.set_yticks(ticks = [0.1,1,10,25])
    ax.set_yticklabels(labels = ['0.1','1','10','25'])
    ax.set_xlabel("Color (BP-RP)")
    ax.set_ylabel("Period (d)")
    ax.legend(loc = 'lower left', framealpha = 0.6, fancybox = False, fontsize=16)
    
    return(ax)
            
def pc_seq_uvwxyz(plot_df, tmag_summary_dict, group_toi_dict, cont_thresh = 0.7, toi_label = None, legend_locs_dict = None, lit_seq = ['Pleiades','Praesepe','Hyades']):
    ### NEEDS UPDATE: variable legend locs labeled by subplot indices
    
    ## create rot avail/not avail subset dataframes for labeled plotting
    tic = group_toi_dict['tic']
    if toi_label is None: toi_label = 'TIC ' + str(tic)
    
    plot_toi = plot_df[plot_df['tic'] == tic]
    plot_df2 = plot_df[np.logical_not(plot_df['contratio'] > cont_thresh)]
    rot_avail = plot_df2[plot_df2['rot_avail'] == True]
    rot_not_avail = plot_df2[plot_df2['rot_avail'] == False]
    
    fig, axs = plt.subplots(nrows = 3, ncols = 3, figsize = (20,20))
    
    ax1 = axs[0,0]
    ff_pc_seq(ax = ax1, plot_df = plot_df, group_toi_dict = group_toi_dict,
                 cont_thresh = cont_thresh, lit_seq = lit_seq)
    
    ax2 = axs[0,1]
    ax2.scatter(rot_not_avail['x'],rot_not_avail['z'], c = 'darkgrey', label = 'No Rot', alpha = 0.8, edgecolors = 'k')
    ax2.scatter(rot_avail['x'],rot_avail['z'], c = 'blue', label = 'Rot', alpha = 0.6, edgecolors = 'k')
    ax2.scatter(plot_toi['x'],plot_toi['z'],c = 'springgreen', s = 150,
                alpha = 1, edgecolors = 'black',
                label = toi_label, marker = 'x')
    ax2.set_xlabel('X [pc]')
    ax2.set_ylabel('Z [pc]')
    ax2.legend(loc = 'upper right', fontsize = 'x-large')
    
    ax3 = axs[0,2]
    ax3.scatter(rot_not_avail['y'],rot_not_avail['z'], c = 'darkgrey', label = 'No Rot', alpha = 0.8, edgecolors = 'k')
    ax3.scatter(rot_avail['y'],rot_avail['z'], c = 'blue', label = 'Rot', alpha = 0.6, edgecolors = 'k')
    ax3.scatter(plot_toi['y'],plot_toi['z'],c = 'springgreen', s = 150,
                alpha = 1, edgecolors = 'black',
                label = toi_label, marker = 'x')
    ax3.set_xlabel('Y [pc]')
    ax3.set_ylabel('Z [pc]')
    ax3.legend(loc = 'upper left', fontsize = 'x-large')
    
    ax4 = axs[1,0]
    ax4.scatter(rot_not_avail['u'],rot_not_avail['w'], c = 'darkgrey', label = 'No Rot', alpha = 0.8, edgecolors = 'k')
    ax4.scatter(rot_avail['u'],rot_avail['w'], c = 'blue', label = 'Rot', alpha = 0.6, edgecolors = 'k')
    ax4.scatter(plot_toi['u'],plot_toi['w'],c = 'springgreen', s = 150,
                alpha = 1, edgecolors = 'black',
                label = toi_label, marker = 'x')
    ax4.set_xlabel('U [km/s]')
    ax4.set_ylabel('W [km/s]')
    ax4.legend(loc = 'upper right', fontsize = 'x-large')
    
    ax5 = axs[1,1]
    ax5.scatter(rot_not_avail['delta_pmra'],rot_not_avail['delta_pmdec'], c = 'darkgrey', label = 'No Rot', alpha = 0.8, edgecolors = 'k')
    ax5.scatter(rot_avail['delta_pmra'],rot_avail['delta_pmdec'], c = 'blue', label = 'Rot', alpha = 0.6, edgecolors = 'k')
    ax5.scatter(plot_toi['delta_pmra'],plot_toi['delta_pmdec'],c = 'springgreen', s = 150,
                alpha = 1, edgecolors = 'black',
                label = toi_label, marker = 'x')
    ax5.set_xlabel(r'$\Delta \mu_{\delta}$ [mas/yr]')
    ax5.set_ylabel(r'$\Delta \mu_{\alpha}$ [mas/yr]')
    ax5.legend(loc = 'upper right', fontsize = 'x-large')
    
    ax6 = axs[1,2]
    ax6.scatter(rot_not_avail['x'],rot_not_avail['y'], c = 'darkgrey', label = 'No Rot', alpha = 0.8, edgecolors = 'k')
    ax6.scatter(rot_avail['x'],rot_avail['y'], c = 'blue', label = 'Rot', alpha = 0.6, edgecolors = 'k')
    ax6.scatter(plot_toi['x'],plot_toi['y'],c = 'springgreen', s = 150,
                alpha = 1, edgecolors = 'black',
                label = toi_label, marker = 'x')
    ax6.set_xlabel('X [pc]')
    ax6.set_ylabel('Y [pc]')
    ax6.legend(loc = 'lower left', fontsize = 'x-large')
    
    ax7 = axs[2,0]
    ax7.scatter(rot_not_avail['v'],rot_not_avail['w'], c = 'darkgrey', label = 'No Rot', alpha = 0.8, edgecolors = 'k')
    ax7.scatter(rot_avail['v'],rot_avail['w'], c = 'blue', label = 'Rot', alpha = 0.6, edgecolors = 'k')
    ax7.scatter(plot_toi['v'],plot_toi['w'],c = 'springgreen', s = 150,
                alpha = 1, edgecolors = 'black',
                label = toi_label, marker = 'x')
    ax7.set_xlabel('V [km/s]')
    ax7.set_ylabel('W [km/s]')
    ax7.legend(loc = 'upper right', fontsize = 'x-large')
    
    ax8 = axs[2,1]
    ax8.scatter(rot_not_avail['u'],rot_not_avail['v'], c = 'darkgrey', label = 'No Rot', alpha = 0.8, edgecolors = 'k')
    ax8.scatter(rot_avail['u'],rot_avail['v'], c = 'blue', label = 'Rot', alpha = 0.6, edgecolors = 'k')
    ax8.scatter(plot_toi['u'],plot_toi['v'],c = 'springgreen', s = 150,
                alpha = 1, edgecolors = 'black',
                label = toi_label, marker = 'x')
    ax8.set_xlabel('U [km/s]')
    ax8.set_ylabel('V [km/s]')
    ax8.legend(loc = 'upper right', fontsize = 'x-large')
    
    ax9 = axs[2,2]
    tmag_plot(ax = ax9,tmag_df = tmag_summary_dict['tmag_table'],
              tmag_summary = tmag_summary_dict['tmag_summary'])
    
    fig.tight_layout()
    return(fig)        
            
def comove_ff_rotations(query_df,download_dir,vlim=5,srad=25,ff_lists_loc = None):
    if ff_lists_loc is None:
        import Comove
        #query gaia for RVs
        gaia_query_res = catQ.get_gaia_data_bulk(query_df = query_df,id_col_name = 'tic').rename(columns = {'ra':'ra_gaia','dec':'dec_gaia'})
        query_df = query_df.merge(right = gaia_query_res, on = 'tic', how = 'left').drop_duplicates(subset=['tic']).reset_index(drop = True)
        
        vlim=vlim #*u.kilometer/u.second ##km/s
        srad=srad ##parsecs (spherical radius around target)
    else: 
        #get names of files
        ff_lists_fns = os.listdir(ff_lists_loc)
        TOIs_in_ff_lists = []
        for fn in ff_lists_fns:
            TOIs_in_ff_lists.append(fn.split(".")[0].split("toi")[1])

    for i,row in query_df.iterrows():
        toi = str(row['TOI']).split(".")[0]
        print("Working on TOI " + toi + ", object "  + str(i+1) + "/" + str(len(query_df)) + ".")
        tic = str(row['tic'])
        targname = "TOI " + toi
        
        #make folder for TOI , and name subfolder for toi comove products
        toi_folder = os.path.join(download_dir,'toi' + toi)
        if os.path.exists(toi_folder) == False: os.mkdir(toi_folder)
        
        ff_product_folder = os.path.join(toi_folder,'ff_products/')
        
        ## run Comove if ff_lists not provided
        if ff_lists_loc is None:
            
            rd = [row['ra']*(24/360),row['dec']]
            radvel = row['dr2_radial_velocity']
            if str(radvel) == 'nan': 
                print("No radial velocity for this object in Gaia. Trying next TOI.")
                continue
            
            try:
                output_location = Comove.findfriends(targname,radvel,velocity_limit=vlim,search_radius=srad,radec=rd,output_directory=ff_product_folder,verbose=False,showplots=False)
            except:
                print("Comove didn't find enough comving friends. On to next TOI.")
                continue
            ## get Comove results
            friends_fn = os.path.join(ff_product_folder,targname.replace(" ","") + ".txt")
            friends_df = pd.read_csv(friends_fn, delim_whitespace=True).rename(columns = {'RA':'ra','DEC':'dec'})
        else:
            #read ff_lists that match TOI name
            if os.path.exists(ff_product_folder) == False: os.mkdir(ff_product_folder)
            if toi in TOIs_in_ff_lists:
                friends_fn = os.path.join(ff_lists_loc,'toi' + toi + '.txt')
                friends_df = pd.read_csv(friends_fn, delim_whitespace=True).rename(columns = {'RA':'ra','DEC':'dec'})
            else:
                print("Could not find friend list for this TOI. Moving to next TOI.")
                continue
            
        
        
        ## limit number of friends if len > 800
        voff_list = [4.5,4,3.5,3,2.5,2]
        for voff in voff_list:
            if len(friends_df) < 800: break
            friends_df = friends_df[friends_df['Voff(km/s)'] < voff].reset_index(drop = True)
            
        ## run group rotations with 'ff' settings
        #create group toi dict
        group_toi_dict = {'tic':tic,'toi':toi}
        
        # run the run function
        ff_group_run(group_toi_dict = group_toi_dict, download_dir = toi_folder,
                     friends_df = friends_df, ff_product_folder = ff_product_folder)
        
def ff_group_run(group_toi_dict,download_dir,friends_df,ff_product_folder):
    toi = group_toi_dict['toi']
    tic = group_toi_dict['tic']
    targname = "TOI " + toi
    
    #create group
    group_fn = os.path.join(download_dir,'TOI' + toi + '_FFgroup.pkl')
    group = Group(name = "TOI " + toi +  " FF", group_df = friends_df, group_toi_dict = group_toi_dict)
    
    # add TIC catalog info
    group.add_TIC_info()
    save_group_object(group,group_fn)
    
    # #add lcs, save group with rotation_dict_collection
    lc_download_dir = os.path.join(download_dir,'lc_pickles')
    if os.path.exists(lc_download_dir) == False: os.mkdir(lc_download_dir)
    # group.add_tess_LCs(download_dir = lc_download_dir, lc_types = ['cpm'])
    # save_group_object(group,group_fn)
    
    #add LCs externally, to update/save group object after each target
    rots_dict_collection = bulk_download(tic_list = group.tics, 
                                                     download_dir = lc_download_dir, 
                                                     lc_types = ['cpm'],
                                                     #spoc_kwrgs = spoc_kwrgs,
                                                     group_obj = group,
                                                     group_fn = group_fn)
    group.rots_dict_collection = rots_dict_collection
    save_group_object(group,group_fn)
    
    # add Gaia query    
    group.add_gaia_info(id_col_name = 'tic')
    save_group_object(group,group_fn)
    
    #organize best_rots, Tmag summary
    group.rot_summary(lc_types = ['cpm'], tmag_list = [14,15,16,99])
    save_group_object(group,group_fn)
    
    #add pc_seq, save group and save rotation images with ff_products
    group.add_pc_seq_fig()
    save_group_object(group,group_fn)
    
    pc_seq_fn = os.path.join(ff_product_folder,targname.replace(" ","") + "_pc_seq.png")
    group.ff_pc_seq.savefig(pc_seq_fn)
    
    pc_seq_uvwxyz_fn = os.path.join(ff_product_folder,targname.replace(" ","") + "_pc_seq_uvwxyz.png")
    group.ff_uvwxyz.savefig(pc_seq_uvwxyz_fn)
    
    #save plotting results
    plot_df_fn = os.path.join(ff_product_folder,targname.replace(" ","") + "_rots.csv")
    group.save_plot_df(fn = plot_df_fn)
    save_group_object(group,group_fn)
    
      