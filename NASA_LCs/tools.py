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

from NASA_LCs import mastQuery


def save_object(target_object, keep_fits = False):
    #saving the file
    target_fn = "tic" + target_object.tic + ".pkl"            
    target_path = os.path.join(target_object.folder_name,target_fn)            
    target_object.target_path = target_path
    
    if keep_fits == False:
        #saving the file
        save_path = os.path.join(target_object.parent_folder,target_fn)
        target_object.target_path = save_path
        with open(save_path,'wb') as outfile:
            pkl.dump(target_object,outfile)
        shutil.rmtree(target_object.folder_name)
    else:
        with open(target_path,'wb') as outfile:
            pkl.dump(target_object,outfile)
        
def bulk_download(tic_list, download_dir, products = ['LC'], run_rotations = True, 
                  rotation_options = {'flux_type':['SAP_flux'],'flux_err_avail':True,'min_freq':1/30},
                  save_objects = True, keep_fits = False):
    
    for i,tic in enumerate(tic_list):
        print("Working on object " + str(i+1) + "/" + str(len(tic_list)) + ".")
        print(tic)
        target_obj = mastQuery.mast_target(tic = tic)
        target_obj.download(products = products, download_dir = download_dir)
        if run_rotations == True: 
            target_obj.run_rotation_tools(flux_type = rotation_options['flux_type'],
                                          flux_err_avail = rotation_options['flux_err_avail'],
                                          min_freq = rotation_options['min_freq'])
        if save_objects == True: save_object(target_object = target_obj, keep_fits = keep_fits)
        
        del target_obj
        