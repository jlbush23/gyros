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

from NASA_LCs.Target import Target

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

def save_target_object(target_object, download_dir):
    #saving the file
    target_fn = "tic" + target_object.tic + ".pkl"            
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
        
def bulk_download(tic_list, download_dir, lc_types = ['spoc','cpm'],
                  run_rotations = True, min_freq = 1/30,
                  #rot_options = {'flux_type':['spoc','cpm'],'flux_err_avail':[True,False],'min_freq':1/30},
                  save_objects = True, keep_fits = False):
    target_dict = {}
    for i,tic in enumerate(tic_list):
        print("Working on object " + str(i+1) + "/" + str(len(tic_list)) + ".")
        #print(tic)
        if (str(tic) == 'nan'):
            print("No TIC for this object. Moving to next.")
            continue
        
        target_obj = Target(tic = tic)
        # try spoc
        if 'spoc' in lc_types:
            try: 
                target_obj.add_lk_LCs()
                if (run_rotations == True) & ('spoc120_lc' in target_obj.available_attributes): target_obj.run_spoc_rots(min_freq = min_freq)
            except:
                print("Data download error. Check if SPOC available on MAST.")
            
        if 'cpm' in lc_types: 
            try: 
                target_obj.add_cpm_LC()
                if (run_rotations == True) & ('cpm_lc' in target_obj.available_attributes): target_obj.run_cpm_rots(min_freq = min_freq)
            except:
                print("Data download error or error with running CPM rotations.")
        
        if save_objects == True: save_target_object(target_object = target_obj, download_dir = download_dir)
        target_dict[str(target_obj.tic)] = target_obj
        del target_obj        
        
    return(target_dict)

def best_tess_rots(target_dict,lc_types = ['spoc','cpm']):
    
    ## somehow feed a dict of rot dicts instead?
    
    if 'spoc' in lc_types: 
        best_sap_rots = []
        best_pdc_rots = []
    if 'cpm' in lc_types: best_cpm_rots =[]
    
    for tic in target_dict.keys():
        targ = target_dict[tic]
        
        if ('cpm_rot_dict' in targ.available_attributes) & ('cpm' in lc_types):
            cpm_LS_res = targ.cpm_rot_dict['LS_res']
            cpm_AC_res = targ.cpm_rot_dict['AC_res']
            
            temp_res = cpm_LS_res.merge(right = cpm_AC_res, on = 'sector', how = 'outer')
            best_idx = np.where(np.max(temp_res['LS_Power1']) == temp_res['LS_Power1'])[0][0]
            best_sector_res = temp_res.iloc[best_idx,:].to_frame().transpose()
            best_sector_res.insert(loc = 0, column = 'tic', value = str(targ.tic))
            best_cpm_rots.append(best_sector_res)

        if ('sap_rot_dict' in targ.available_attributes) & ('spoc' in lc_types):
            sap_LS_res = targ.sap_rot_dict['LS_res']
            sap_AC_res = targ.sap_rot_dict['AC_res']
            
            temp_res = sap_LS_res.merge(right = sap_AC_res, on = 'sector', how = 'outer')
            best_idx = np.where(np.max(temp_res['LS_Power1']) == temp_res['LS_Power1'])[0][0]
            best_sector_res = temp_res.iloc[best_idx,:].to_frame().transpose()
            best_sector_res.insert(loc = 0, column = 'tic', value = str(targ.tic))
            best_sap_rots.append(best_sector_res)
            
        if ('pdc_rot_dict' in targ.available_attributes) & ('spoc' in lc_types):
            pdc_LS_res = targ.pdc_rot_dict['LS_res']
            pdc_AC_res = targ.pdc_rot_dict['AC_res']
            
            temp_res = pdc_LS_res.merge(right = pdc_AC_res, on = 'sector', how = 'outer')
            best_idx = np.where(np.max(temp_res['LS_Power1']) == temp_res['LS_Power1'])[0][0]
            best_sector_res = temp_res.iloc[best_idx,:].to_frame().transpose()
            best_sector_res.insert(loc = 0, column = 'tic', value = str(targ.tic))
            best_pdc_rots.append(best_sector_res)            
    
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
        best_sap_rots_df = pd.concat(best_sap_rots)
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

        best_pdc_rots_df = pd.concat(best_pdc_rots)
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
        
        best_rots_dict['sap'] = best_sap_rots_df
        best_rots_dict['pdc'] = best_pdc_rots_df
        
        
    return(best_rots_dict)

def pc_seq_fig(ax, color_type = 'bp_rp', pleiades_on = True, praesepe_on = True, hyades_on = True, upper_sco_on = False, xlim = (0.05,3.5), guidelines_on = False):
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    import pandas as pd
    import numpy as np
    #This stuff is everything, use it for any python plot to make it nicer.
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
            
        
            
        