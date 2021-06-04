# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 11:54:09 2020

@author: jlbus
"""
from astroquery.mast import Observations
from astroquery.mast import Catalogs
from astroquery.gaia import Gaia
#from astroquery.vizier import Vizier
#from astroquery.simbad import Simbad


from astropy import units as u
from astropy.coordinates import SkyCoord

import pandas as pd
import numpy as np

import os

from uvwxyz.xyzuvw import xyz,uvw


def get_coord_from_ID(id_type, ID):    
    if id_type == 'TIC':
        query_string = 'tic ' + str(ID)
        obs_table = Catalogs.query_object(query_string, radius = 0.002*u.deg, catalog = 'TIC')
        obs_table = obs_table[obs_table['ID'] == str(ID)]
        ra = obs_table['ra'][0]
        dec = obs_table['dec'][0]
        return(ra,dec)
    # else:
    #     obs_table = Simbad.query_object(id_label + str(row[id_label]))
    #     try:
    #         if len(obs_table) == 1:
    #             ra = ra2deg(obs_table['RA'][0])
    #             dec = dec2deg(obs_table['DEC'][0])
    #             print("Found RA = " + str(ra) + " & Dec = " + str(dec) + "!")
    #             return(ra,dec)
    #     except:
    #         return(np.nan,np.nan)
    
    

# def get_this_stuff(self,get = ['ra','dec','tic','epic','tmag','gaia','spoc_avail'], gaia_kwrgs = None, id_label = None, roll = True):
#     self.df_update = pd.DataFrame(columns = np.append(self.df_cols,get))
#     for i,row in self.df.iterrows():
#         print("Working on object " + str(i+1) + "/" + str(len(self.df)) + ".")
#         if 'ra' in get:
#             if (np.isnan(float(row['tic'])) == False):
#                 print("Working on coordinates for object " + str(i+1) + "/" + str(len(self.df)) + ".")
#                 row['ra'],row['dec'] = self.get_coord_from_identifier(row,id_label = id_label)
#             else:
#                 row['ra'] = np.nan
#                 row['dec'] = np.nan
#         if ('epic' in get) & (np.isnan(row['ra']) == False):
#             print("Working on EPIC for object " + str(i+1) + "/" + str(len(self.df)) + ".")
#             row['epic'] = self.get_epic(row)
#         if ('tic' in get) & (np.isnan(row['ra']) == False):
#             print("Working on TIC for object " + str(i+1) + "/" + str(len(self.df)) + ".")
#             if ('gaia' in self.df_cols) & (np.isnan(float(row['gaia'])) == False):
#                 row['tic'] = self.get_tic_w_gaia(row)
#             else:
#                 row['tic'] = self.get_tic(row)
#         if ('spoc_avail' in get):
#             if (np.isnan(float(row['tic'])) == False):
#                 row['spoc_avail'] = self.check_spoc(row)
#             else:
#                 row['spoc_avail'] = False
#         if ('kic' in get) & (np.isnan(row['ra']) == False):
#             print("Working on KIC for object " + str(i+1) + "/" + str(len(self.df)) + ".")
#             row['kic'] = self.get_kic(row)
#         if ('gaia' in get) & (np.isnan(row['ra']) == False):
#             print("Working on GAIA ID for object " + str(i+1) + "/" + str(len(self.df)) + ".")
#             row['gaia'] = self.get_gaia_id(row)
#         if (gaia_kwrgs is not None):
#             if ('gaia' in get) & (np.isnan(float(row['gaia'])) == False):
#                 row = row.append(self.get_gaia_data(row,gaia_kwrgs))
#             else:
#                 for col in gaia_kwrgs: row[col] = np.nan
#         if 'tmag' in get:
#             if (np.isnan(float(row['tic'])) == False):
#                 row['tmag'] = self.get_tmag(row)
#             else:
#                 row['tmag'] = np.nan
            
#         self.df_update = self.df_update.append(row, ignore_index = True)
#         self.df_update.to_csv(self.save_fn, index = False)
#     if roll:
#         self.df_update = self.roll_columns(self.df_update,num_roll = len(get))
#     self.df_update.to_csv(self.save_fn, index = False)
    

# #the complete epic catalogue is on Vizier
# def get_epic(self,row):
#     radii = np.linspace(start = 0.0001, stop = 0.001, num = 19)
#     #for i,row in self.df.iterrows():
#     epic_found = False
#     for rad in radii:
#         if epic_found == False:
#             query_string = str(row['ra']) + " " + str(row['dec']) # make sure to have a space between the strings!
#             obs_table_list = Vizier.query_region(coordinates = query_string, radius = rad*u.deg, catalog = "EPIC");
#             if len(obs_table_list) == 1:
#                 epic = obs_table_list[0]['ID'][0]
#                 epic_found = True
#                 continue
#     if epic_found == False:
#         epic = np.nan
#         print("Didn't find EPIC for this object.")
#     else:
#         print("Found EPIC " + str(epic) + "!")
#     return(epic)


def get_tic_bulk(query_df,ra_col_name='ra',dec_col_name='dec'):
    tics = []
    for i,row in query_df.iterrows():
        print("Finding TIC for object " + str(i+1) + "/" + str(len(query_df)) + ".")
        if (str(row[ra_col_name]) == 'nan') | (str(row[dec_col_name]) == 'nan'):
            print("No coordinates provided for this object.")
            tics.append(np.nan)
        else:
            tics.append(get_tic(ra = row[ra_col_name],dec = row[dec_col_name]))
    
    query_df.insert(loc = 0, column = 'tic', value = tics)
    return(tics)#,query_df)
    
#is the MAST website built on the TIC? 'ID' is the TIC ID for every query     
def get_tic(ra,dec):
    radii = np.linspace(start = 0.0001, stop = 0.001, num = 19)
    #for i,row in self.df.iterrows():
    tic_found = False
    try:
        for rad in radii:
            if tic_found == False:
                query_string = str(ra) + " " + str(dec) # make sure to have a space between the strings!#SkyCoord(ra = row['ra'], dec = row['dec'], frame = 'icrs') str(row['ra']) + " " + str(row['dec']) # make sure to have a space between the strings!
                obs_table = Catalogs.query_region(coordinates = query_string, radius = rad*u.deg, catalog = "TIC")
                obs_df = obs_table.to_pandas()
                if len(obs_table['ID']) == 1:
                    tic = obs_table['ID'][0]
                    tic_found = True
                    continue
                if len(obs_df[obs_df['GAIA'].to_numpy(dtype = 'str') != '']) == 1:
                    temp_obs_df = obs_df[obs_df['GAIA'].to_numpy(dtype = 'str') != '']
                    tic = temp_obs_df['ID'].iloc[0]
                    tic_found = True
                    continue
                if len(np.unique(obs_df[obs_df['HIP'].to_numpy(dtype = 'str') != '']['HIP'])) == 1:
                    tic = obs_table['ID'][0]
                    tic_found = True
                    continue
        if tic_found == False:
            tic = np.nan
            print("Didn't find TIC for this object.")
        else:
            print("Found TIC " + str(tic) + "!")
            #self.tic = tic
    except:
        print("Issue finding TIC for this object.")
        tic = np.nan
    return(tic)

#looks like the TIC has a KIC column... hopefully it includes the full KIC! 
def get_kic_bulk(query_df,ra_col_name='ra',dec_col_name='dec'):
    kics = []
    for i,row in query_df.iterrows():
        print("Finding KIC for object " + str(i+1) + "/" + str(len(query_df)) + ".")
        if (str(row[ra_col_name]) == 'nan') | (str(row[dec_col_name]) == 'nan'):
            print("No coordinates provided for this object.")
            kics.append(np.nan)
        else:
            kics.append(get_kic(ra = row[ra_col_name],dec = row[dec_col_name]))
    
    query_df.insert(loc = 0, column = 'tic', value = tics)
    return(tics)#,query_df)
  
def get_kic(ra,dec):
    radii = np.linspace(start = 0.0001, stop = 0.001, num = 19)
    #for i,row in self.df.iterrows():
    kic_found = False
    for rad in radii:
        if kic_found == False:
            query_string = str(ra) + " " + str(dec) # make sure to have a space between the strings!#SkyCoord(ra = row['ra'], dec = row['dec'], frame = 'icrs') str(row['ra']) + " " + str(row['dec']) # make sure to have a space between the strings!
            obs_table = Catalogs.query_region(coordinates = query_string, radius = rad*u.deg, catalog = "TIC")
            obs_df = obs_table.to_pandas()
            if len(obs_table['ID']) == 1:
                kic = obs_table['KIC'][0]
                kic_found = True
                continue
            if len(obs_df[obs_df['GAIA'].to_numpy(dtype = 'str') != '']) == 1:
                temp_obs_df = obs_df[obs_df['GAIA'].to_numpy(dtype = 'str') != '']
                kic = temp_obs_df['KIC'].iloc[0]
                kic_found = True
                continue
            if len(np.unique(obs_df[obs_df['HIP'].to_numpy(dtype = 'str') != '']['HIP'])) == 1:
                kic = obs_table['KIC'][0]
                kic_found = True
                continue
    if kic_found == False:
        kic = np.nan
        print("Didn't find KIC for this object.")
    else:
        print("Found KIC " + str(kic) + "!")
        #self.tic = tic
    return(kic)

def get_TIC_data_bulk(query_df,ra_col_name = 'ra', dec_col_name = 'dec',append_tics = True):
    tics = []
    TIC_dfs = []
    for i,row in query_df.iterrows():
        print("Finding TIC data for object " + str(i+1) + "/" + str(len(query_df)) + ".")
        if (str(row[ra_col_name]) == 'nan') | (str(row[dec_col_name]) == 'nan'):
            print("No coordinates provided for this object.")
            tics.append(np.nan)
            TIC_dfs.append(pd.DataFrame())
        else:
            try:
                temp_TIC_df, temp_tic = get_TIC_data(ra = row[ra_col_name],dec = row[dec_col_name])
                tics.append(temp_tic)
                TIC_dfs.append(temp_TIC_df)
            except:
                print("Issue getting TIC data for this object.")
                tics.append(np.nan)
                TIC_dfs.append(pd.DataFrame(data = {ra_col_name:[row[ra_col_name]],dec_col_name:[row[dec_col_name]]}))
    
    if append_tics == True: query_df.insert(loc = 0, column = 'tic', value = np.array(tics, dtype = 'str'))
    TIC_query = pd.concat(TIC_dfs)
    return(tics, TIC_query)#,query_df)

def get_TIC_data(ra,dec):
    radii = np.linspace(start = 0.0001, stop = 0.001, num = 19)
    #for i,row in self.df.iterrows():
    tic_found = False
    for rad in radii:
        if tic_found == False:
            query_string = str(ra) + " " + str(dec) # make sure to have a space between the strings!#SkyCoord(ra = row['ra'], dec = row['dec'], frame = 'icrs') str(row['ra']) + " " + str(row['dec']) # make sure to have a space between the strings!
            obs_table = Catalogs.query_region(coordinates = query_string, radius = rad*u.deg, catalog = "TIC")
            obs_df = obs_table.to_pandas()
            if len(obs_table['ID']) == 1:
                tic = obs_table['ID'][0]
                tic_found = True
                continue
            if len(obs_df[obs_df['GAIA'].to_numpy(dtype = 'str') != '']) == 1:
                temp_obs_df = obs_df[obs_df['GAIA'].to_numpy(dtype = 'str') != '']
                tic = temp_obs_df['ID'].iloc[0]
                tic_found = True
                continue
            if len(np.unique(obs_df[obs_df['HIP'].to_numpy(dtype = 'str') != '']['HIP'])) == 1:
                tic = obs_table['ID'][0]
                tic_found = True
                continue
    if tic_found == False:
        tic = np.nan
        print("Didn't find TIC for this object.")
    else:
        print("Found TIC Data for TIC " + str(tic) + "!")
        #self.tic = tic
    return(obs_df,tic)

# def get_tic_w_gaia(self,row):
#     rad = 0.005
#     #for i,row in self.df.iterrows():
#     tic_found = False
#     query_string = str(row['ra']) + " " + str(row['dec']) # make sure to have a space between the strings!#SkyCoord(ra = row['ra'], dec = row['dec'], frame = 'icrs') str(row['ra']) + " " + str(row['dec']) # make sure to have a space between the strings!
#     obs_table = Catalogs.query_region(coordinates = query_string, radius = rad*u.deg, catalog = "TIC")
#     obs_df = self.tab2df(obs_table)
#     obs_df = obs_df[obs_df['GAIA'] == row['gaia']]
#     if len(obs_df) == 1:
#         tic = str(obs_df['ID'].to_numpy(dtype = 'str')[0])
#         tic_found = True
    
#     if tic_found == False:
#         tic = np.nan
#         print("Didn't find TIC for this object.")
#     else:
#         print("Found TIC " + str(tic) + "!")
#         #self.tic = tic
#     return(tic)
    

# def check_spoc(self,row):
#     '''
#     NEEDS LIGHTKURVE UPDATE

#     Parameters
#     ----------
#     row : TYPE
#         DESCRIPTION.

#     Returns
#     -------
#     None.

#     '''
#     query_string = "tic " + row['tic']
 
#     try:
#         obs_table = Observations.query_object(query_string)#, radius = 0.0005*u.deg)        
#         SPOC_df = obs_table[obs_table['provenance_name'] == 'SPOC'].to_pandas()
#         timeseries_ids = SPOC_df[SPOC_df['dataproduct_type'] == 'timeseries']['obsid'].to_numpy(dtype='str')            
#         if len(timeseries_ids) > 0:
#             spoc_avail = True
#             print("SPOC is available!")
#         else:
#             spoc_avail = False
#             print("SPOC is NOT available :-(")
#     except:
#         spoc_avail = False
#         print("SPOC is NOT available :-(")
        
#     return(spoc_avail)
        
        
# def get_tmag(self,row):
#     #need tic!
#     query_string = 'tic ' + str(row['tic'])
#     obs_table = Catalogs.query_object(query_string, radius = 0.0005*u.degree, catalog = "TIC")
#     obs_df = obs_table.to_pandas()
#     obs_df = obs_df[obs_df['ID'] == str(row['tic'])]
#     return(obs_df['Tmag'][0])        

# #only partial KIC is on vizier. total KIC is on MAST
# def get_kic(self,row):
#     radii = np.linspace(start = 0.0001, stop = 0.001, num = 19)
#     #for i,row in self.df.iterrows():
#     kic_found = False
#     for rad in radii:
#         if kic_found == False:
#             query_string = str(row['ra']) + " " + str(row['dec']) # make sure to have a space between the strings!#SkyCoord(ra = row['ra'], dec = row['dec'], frame = 'icrs') str(row['ra']) + " " + str(row['dec']) # make sure to have a space between the strings!
#             obs_table = Catalogs.query_region(coordinates = query_string, radius = rad*u.deg, catalog = "TIC")
#             obs_df = self.tab2df(obs_table)
#             if len(obs_table['ID']) == 1:
#                 kic = obs_table['KIC'][0]
#                 if kic == '--':
#                     kic_found = False
#                     continue
#                 else:
#                     kic_found = True
#                     continue
#             if len(obs_df[obs_df['GAIA'].to_numpy(dtype = 'str') != '']) == 1:
#                 temp_obs_df = obs_df[obs_df['GAIA'].to_numpy(dtype = 'str') != '']
#                 kic = temp_obs_df['KIC'].iloc[0]
#                 if kic == '--':
#                     kic_found = False
#                     continue
#                 else:
#                     kic_found = True
#                     continue
#             if len(np.unique(obs_df[obs_df['HIP'].to_numpy(dtype = 'str') != '']['HIP'])) == 1:
#                 kic = obs_table['KIC'][0]
#                 if kic == '--':
#                     kic_found = False
#                     continue
#                 else:
#                     kic_found = True
#                     continue
#     if kic_found == False:
#         kic = np.nan
#         print("Didn't find KIC for this object.")
#     else:
#         print("Found KIC " + str(kic) + "!")
#     return(kic)

def get_gaia_data_bulk(query_df,ra_col_name='ra',dec_col_name='dec',gaia_kwrgs = 'all',id_col_name = None):
    
    gaia_df_list = []
    for i,row in query_df.iterrows():
        print("Working on Gaia data for object " + str(i+1) + "/" + str(len(query_df)) + ".")
        
        temp_query_df = get_gaia_data(ra = row[ra_col_name],dec = row[dec_col_name],gaia_kwrgs=gaia_kwrgs)
        if id_col_name is not None: temp_query_df.insert(loc = 0, column = id_col_name, value = row[id_col_name])
            
        gaia_df_list.append(temp_query_df)
        
    gaia_query_df = pd.concat(gaia_df_list)
    
    return(gaia_query_df)
    

#most? TICs have an associated GAIA. possibly remove catalog filter to get more gaia
def get_gaia_id(ra,dec):
    radii = np.linspace(start = 0.0001, stop = 0.001, num = 19)
    #for i,row in self.df.iterrows():
    gaia_id_found = False
    try:
        for rad in radii:
            if gaia_id_found == False:
                query_string = str(ra) + " " + str(dec) # make sure to have a space between the strings!#SkyCoord(ra = row['ra'], dec = row['dec'], frame = 'icrs') str(row['ra']) + " " + str(row['dec']) # make sure to have a space between the strings!
                obs_table = Catalogs.query_region(coordinates = query_string, radius = rad*u.deg, catalog = "TIC")
                obs_df = obs_table.to_pandas()
                if len(obs_table['ID']) == 1:
                    gaia_id = str(obs_table['GAIA'][0])
                    gaia_id_found = True
                    continue
                if len(obs_df[obs_df['GAIA'].to_numpy(dtype = 'str') != '']) == 1:
                    temp_obs_df = obs_df[obs_df['GAIA'].to_numpy(dtype = 'str') != '']
                    gaia_id = str(temp_obs_df['GAIA'].iloc[0])
                    gaia_id_found = True
                    continue
                if len(np.unique(obs_df[obs_df['HIP'].to_numpy(dtype = 'str') != '']['HIP'])) == 1:
                    gaia_id = str(obs_table['GAIA'][0])
                    gaia_id_found = True
                    continue
        if gaia_id_found == False:
            gaia_id = np.nan
            print("Didn't find GAIA ID for this object.")
        else:
            if gaia_id == '--':
                gaia_id = np.nan
                print("Found object, but no Gaia ID availabel through TIC.")
            else:
                print("Found GAIA ID " + str(gaia_id) + "!")
    except:
        print("Issue finding GAIA ID for this object.")
        gaia_id = np.nan
    return(gaia_id)

def get_gaia_data(ra,dec,gaia_kwrgs):
    
    #print("Working on GAIA data.")
    c = SkyCoord(ra=ra*u.degree,dec=dec*u.degree)
    # Pgaia = Gaia.query_object_async(coordinate=c, radius=(5.0*u.arcsec))
    sqltext = "SELECT * FROM gaiaedr3.gaia_source WHERE CONTAINS( \
               POINT('ICRS',gaiaedr3.gaia_source.ra,gaiaedr3.gaia_source.dec), \
               CIRCLE('ICRS'," + str(c.ra.value) +","+ str(c.dec.value) +","+ str(6.0/3600.0) +"))=1;"
    job = Gaia.launch_job_async(sqltext , dump_to_file=False)
    obs_table_gaia = job.get_results()
    #obs_table_gaia = Gaia.query_object(c,radius = 0.002*u.degree)
    
    gaia_id = get_gaia_id(ra=ra,dec=dec) #gets gaia id from TIC catalog on MAST to confirm correct gaia query
    if str(gaia_id) != 'nan':
        obs_table_gaia = obs_table_gaia[obs_table_gaia['source_id'] == int(gaia_id)]
    else:
        obs_table_gaia = []
    
    if len(obs_table_gaia) == 1: #query_list1.append(obs_table_gaia.to_pandas())
        print("Found GAIA data!")
        if gaia_kwrgs == 'all':
            return(obs_table_gaia.to_pandas())
        else:
            return(obs_table_gaia[gaia_kwrgs].to_pandas())
    else:
        print("Issue cross-matching Gaia data.")
        if gaia_kwrgs == 'all':
            all_cols = ['solution_id', 'designation', 'source_id', 'random_index', 'ref_epoch',
                        'ra', 'ra_error', 'dec', 'dec_error', 'parallax', 'parallax_error',
                        'parallax_over_error', 'pmra', 'pmra_error', 'pmdec', 'pmdec_error',
                        'ra_dec_corr', 'ra_parallax_corr', 'ra_pmra_corr', 'ra_pmdec_corr',
                        'dec_parallax_corr', 'dec_pmra_corr', 'dec_pmdec_corr',
                        'parallax_pmra_corr', 'parallax_pmdec_corr', 'pmra_pmdec_corr',
                        'astrometric_n_obs_al', 'astrometric_n_obs_ac',
                        'astrometric_n_good_obs_al', 'astrometric_n_bad_obs_al',
                        'astrometric_gof_al', 'astrometric_chi2_al', 'astrometric_excess_noise',
                        'astrometric_excess_noise_sig', 'astrometric_params_solved',
                        'astrometric_primary_flag', 'astrometric_weight_al',
                        'astrometric_pseudo_colour', 'astrometric_pseudo_colour_error',
                        'mean_varpi_factor_al', 'astrometric_matched_observations',
                        'visibility_periods_used', 'astrometric_sigma5d_max',
                        'frame_rotator_object_type', 'matched_observations',
                        'duplicated_source', 'phot_g_n_obs', 'phot_g_mean_flux',
                        'phot_g_mean_flux_error', 'phot_g_mean_flux_over_error',
                        'phot_g_mean_mag', 'phot_bp_n_obs', 'phot_bp_mean_flux',
                        'phot_bp_mean_flux_error', 'phot_bp_mean_flux_over_error',
                        'phot_bp_mean_mag', 'phot_rp_n_obs', 'phot_rp_mean_flux',
                        'phot_rp_mean_flux_error', 'phot_rp_mean_flux_over_error',
                        'phot_rp_mean_mag', 'phot_bp_rp_excess_factor', 'phot_proc_mode',
                        'bp_rp', 'bp_g', 'g_rp', 'radial_velocity', 'radial_velocity_error',
                        'rv_nb_transits', 'rv_template_teff', 'rv_template_logg',
                        'rv_template_fe_h', 'phot_variable_flag', 'l', 'b', 'ecl_lon',
                        'ecl_lat', 'priam_flags', 'teff_val', 'teff_percentile_lower',
                        'teff_percentile_upper', 'a_g_val', 'a_g_percentile_lower',
                        'a_g_percentile_upper', 'e_bp_min_rp_val',
                        'e_bp_min_rp_percentile_lower', 'e_bp_min_rp_percentile_upper',
                        'flame_flags', 'radius_val', 'radius_percentile_lower',
                        'radius_percentile_upper', 'lum_val', 'lum_percentile_lower',
                        'lum_percentile_upper', 'datalink_url', 'epoch_photometry_url', 'dist']
            gaia_nan = [np.repeat(a = np.nan, repeats = len(all_cols))]
            gaia_df = pd.DataFrame(data = gaia_nan, columns = all_cols)
        else:
            gaia_nan = [np.repeat(a = np.nan, repeats = len(gaia_kwrgs))]
            gaia_df = pd.DataFrame(data = gaia_nan, columns = gaia_kwrgs)
        return(gaia_df)
    
def add_gaia_galactic_coords(gaia_query, tic = None):
    ### UPDATE: this function to take optional reference target
    # gaia_query = gaia_query_df
    
    if tic is not None:
        ref_pmra = gaia_query[gaia_query['tic'] == tic]['pmra'][0]
        ref_pmdec = gaia_query[gaia_query['tic'] == tic]['pmdec'][0]
    
    #parallax is in mas, d = 1000/parallax
    def xyz_cols(data,col):
        x,y,z=  xyz(ra = data['ra'], dec = data['dec'], d = 1000/data['parallax'])
        if col == 'x': return(x)
        if col == 'y': return(y)
        if col == 'z': return(z)
        
    def uvw_cols(data,col):
        u,v,w = uvw(ra = data['ra'], dec = data['dec'], d = 1000/data['parallax'],
                    pmra = data['pmra'], pmde = data['pmdec'], rv = data['dr2_radial_velocity'])
        if col == 'u': return(u)
        if col == 'v': return(v)
        if col == 'w': return(w)
        
    def delta_pm(data,col,ref = (np.nan,np.nan)):
        delta_pmra = ref[0] - data['pmra']
        delta_pmdec = ref[1] - data['pmdec']
        
        if col == 'delta_pmra': return(delta_pmra)
        if col == 'delta_pmdec': return(delta_pmdec)
        
    def vtan(data):
        vtan = 4.74 * np.sqrt(data['pmra']**2 + data['pmdec']**2) * (1000/data['parallax'])
        vtan_km_s = vtan/1000
        return(vtan_km_s)        
        
    cols = ['x','y','z']
    for col in cols:
        gaia_query[col] = gaia_query.apply(func = xyz_cols, axis = 1, args = (col))
        
    cols = ['u','v','w']
    for col in cols:
        gaia_query[col] = gaia_query.apply(func = uvw_cols, axis = 1, args = (col))
        
    gaia_query['Vtan(km/s)'] = gaia_query.apply(func = vtan, axis = 1)
    
    if tic is not None:
        cols = ['delta_pmra','delta_pmdec']
        for col in cols:
            gaia_query[col] = gaia_query.apply(func = delta_pm, axis = 1, args = (col,(ref_pmra,ref_pmdec)))
        
    return(gaia_query)

# def roll_columns(self,df,num_roll):
#     cols_rolled = np.roll(df.columns.to_numpy(dtype = 'str'),shift = num_roll)
#     return(df[cols_rolled])    
    
def ra2deg(RA):
    hours,minutes,seconds = RA.split(' ',maxsplit=2)
    decimal_hours = float(hours) + float(minutes)/60 + float(seconds)/3600
    ra = decimal_hours*(360/24)
    return(ra)
def dec2deg(DEC):
    degrees,minutes,seconds = DEC.split(' ',maxsplit=2)
    pos_or_neg = degrees[0]
    degrees = degrees.split(pos_or_neg,maxsplit=1)[1]
    
    dec = abs(float(degrees)) + float(minutes)/60 + float(seconds)/3600
    if pos_or_neg == '-':
        dec = -1*dec
    return(dec) 
