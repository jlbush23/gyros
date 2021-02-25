# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 11:54:09 2020

@author: jlbus
"""
from astroquery.mast import Observations
from astroquery.mast import Catalogs
from astroquery.gaia import Gaia
from astroquery.vizier import Vizier
from astroquery.simbad import Simbad


from astropy import units as u
from astropy.coordinates import SkyCoord

import pandas as pd
import numpy as np

import os


class queryCATS:
    def __init__(self,df,fn):
        self.df = df
        self.df_cols = df.columns.to_numpy(dtype = 'str')
        self.fn = fn
        self.save_fn = os.path.splitext(fn)[0] + '_v2' + os.path.splitext(fn)[1]
        
        if 'ra' in self.df_cols:
            self.coords_avail = True
            if type(self.df['ra']) == str:
                self.coords_convert = True
            else:
                self.coords_convert = False
        else:
            self.coords_avail = False
        
        if 'tic' in self.df_cols:
            self.tic_avail = True
        else:
            self.tic_avail = False
            
        if 'epic' in self.df_cols:
            self.epic_avail = True
        else:
            self.epic_avail = False
            
        if 'bp_rp' in self.df_cols:
            self.gaia_avail = True
        else:
            self.gaia_avail = False
            
    def get_my_stuff(self,get = ['ra','dec','tic','epic','tmag','gaia','spoc_avail'], gaia_kwrgs = None, id_label = None, roll = True):
        self.df_update = pd.DataFrame(columns = np.append(self.df_cols,get))
        for i,row in self.df.iterrows():
            print("Working on object " + str(i+1) + "/" + str(len(self.df)) + ".")
            if 'ra' in get:
                if (np.isnan(float(row['tic'])) == False):
                    print("Working on coordinates for object " + str(i+1) + "/" + str(len(self.df)) + ".")
                    row['ra'],row['dec'] = self.get_coord_from_identifier(row,id_label = id_label)
                else:
                    row['ra'] = np.nan
                    row['dec'] = np.nan
            if ('epic' in get) & (np.isnan(row['ra']) == False):
                print("Working on EPIC for object " + str(i+1) + "/" + str(len(self.df)) + ".")
                row['epic'] = self.get_epic(row)
            if ('tic' in get) & (np.isnan(row['ra']) == False):
                print("Working on TIC for object " + str(i+1) + "/" + str(len(self.df)) + ".")
                if ('gaia' in self.df_cols) & (np.isnan(float(row['gaia'])) == False):
                    row['tic'] = self.get_tic_w_gaia(row)
                else:
                    row['tic'] = self.get_tic(row)
            if ('spoc_avail' in get):
                if (np.isnan(float(row['tic'])) == False):
                    row['spoc_avail'] = self.check_spoc(row)
                else:
                    row['spoc_avail'] = False
            if ('kic' in get) & (np.isnan(row['ra']) == False):
                print("Working on KIC for object " + str(i+1) + "/" + str(len(self.df)) + ".")
                row['kic'] = self.get_kic(row)
            if ('gaia' in get) & (np.isnan(row['ra']) == False):
                print("Working on GAIA ID for object " + str(i+1) + "/" + str(len(self.df)) + ".")
                row['gaia'] = self.get_gaia_id(row)
            if (gaia_kwrgs is not None):
                if ('gaia' in get) & (np.isnan(float(row['gaia'])) == False):
                    row = row.append(self.get_gaia_data(row,gaia_kwrgs))
                else:
                    for col in gaia_kwrgs: row[col] = np.nan
            if 'tmag' in get:
                if (np.isnan(float(row['tic'])) == False):
                    row['tmag'] = self.get_tmag(row)
                else:
                    row['tmag'] = np.nan
                
            self.df_update = self.df_update.append(row, ignore_index = True)
            self.df_update.to_csv(self.save_fn, index = False)
        if roll:
            self.df_update = self.roll_columns(self.df_update,num_roll = len(get))
        self.df_update.to_csv(self.save_fn, index = False)
        
    def get_coord_from_identifier(self,row,id_label):
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
        
        if id_label == 'tic':
            query_string = 'tic ' + row['tic']
            obs_table = Catalogs.query_object(query_string, radius = 0.002*u.deg, catalog = 'TIC')
            obs_table = obs_table[obs_table['ID'] == str(row['tic'])]
            ra = obs_table['ra'][0]
            dec = obs_table['dec'][0]
            return(ra,dec)
        else:
            obs_table = Simbad.query_object(id_label + str(row[id_label]))
            try:
                if len(obs_table) == 1:
                    ra = ra2deg(obs_table['RA'][0])
                    dec = dec2deg(obs_table['DEC'][0])
                    print("Found RA = " + str(ra) + " & Dec = " + str(dec) + "!")
                    return(ra,dec)
            except:
                return(np.nan,np.nan)
    
    
    #the complete epic catalogue is on Vizier
    def get_epic(self,row):
        radii = np.linspace(start = 0.0001, stop = 0.001, num = 19)
        #for i,row in self.df.iterrows():
        epic_found = False
        for rad in radii:
            if epic_found == False:
                query_string = str(row['ra']) + " " + str(row['dec']) # make sure to have a space between the strings!
                obs_table_list = Vizier.query_region(coordinates = query_string, radius = rad*u.deg, catalog = "EPIC");
                if len(obs_table_list) == 1:
                    epic = obs_table_list[0]['ID'][0]
                    epic_found = True
                    continue
        if epic_found == False:
            epic = np.nan
            print("Didn't find EPIC for this object.")
        else:
            print("Found EPIC " + str(epic) + "!")
        return(epic)
    
    #is the MAST website built on the TIC? 'ID' is the TIC ID for every query     
    def get_tic(self,row):
        radii = np.linspace(start = 0.0001, stop = 0.001, num = 19)
        #for i,row in self.df.iterrows():
        tic_found = False
        for rad in radii:
            if tic_found == False:
                query_string = str(row['ra']) + " " + str(row['dec']) # make sure to have a space between the strings!#SkyCoord(ra = row['ra'], dec = row['dec'], frame = 'icrs') str(row['ra']) + " " + str(row['dec']) # make sure to have a space between the strings!
                obs_table = Catalogs.query_region(coordinates = query_string, radius = rad*u.deg, catalog = "TIC")
                obs_df = self.tab2df(obs_table)
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
        return(tic)
    
    def get_tic_w_gaia(self,row):
        rad = 0.005
        #for i,row in self.df.iterrows():
        tic_found = False
        query_string = str(row['ra']) + " " + str(row['dec']) # make sure to have a space between the strings!#SkyCoord(ra = row['ra'], dec = row['dec'], frame = 'icrs') str(row['ra']) + " " + str(row['dec']) # make sure to have a space between the strings!
        obs_table = Catalogs.query_region(coordinates = query_string, radius = rad*u.deg, catalog = "TIC")
        obs_df = self.tab2df(obs_table)
        obs_df = obs_df[obs_df['GAIA'] == row['gaia']]
        if len(obs_df) == 1:
            tic = str(obs_df['ID'].to_numpy(dtype = 'str')[0])
            tic_found = True
        
        if tic_found == False:
            tic = np.nan
            print("Didn't find TIC for this object.")
        else:
            print("Found TIC " + str(tic) + "!")
            #self.tic = tic
        return(tic)
        
    
    def check_spoc(self,row):
        
        query_string = "tic " + row['tic']
     
        try:
            obs_table = Observations.query_object(query_string)#, radius = 0.0005*u.deg)        
            SPOC_df = obs_table[obs_table['provenance_name'] == 'SPOC'].to_pandas()
            timeseries_ids = SPOC_df[SPOC_df['dataproduct_type'] == 'timeseries']['obsid'].to_numpy(dtype='str')            
            if len(timeseries_ids) > 0:
                spoc_avail = True
                print("SPOC is available!")
            else:
                spoc_avail = False
                print("SPOC is NOT available :-(")
        except:
            spoc_avail = False
            print("SPOC is NOT available :-(")
            
        return(spoc_avail)
            
            
    def get_tmag(self,row):
        #need tic!
        query_string = 'tic ' + str(row['tic'])
        obs_table = Catalogs.query_object(query_string, radius = 0.0005*u.degree, catalog = "TIC")
        obs_df = obs_table.to_pandas()
        obs_df = obs_df[obs_df['ID'] == str(row['tic'])]
        return(obs_df['Tmag'][0])        
    
    #only partial KIC is on vizier. total KIC is on MAST
    def get_kic(self,row):
        radii = np.linspace(start = 0.0001, stop = 0.001, num = 19)
        #for i,row in self.df.iterrows():
        kic_found = False
        for rad in radii:
            if kic_found == False:
                query_string = str(row['ra']) + " " + str(row['dec']) # make sure to have a space between the strings!#SkyCoord(ra = row['ra'], dec = row['dec'], frame = 'icrs') str(row['ra']) + " " + str(row['dec']) # make sure to have a space between the strings!
                obs_table = Catalogs.query_region(coordinates = query_string, radius = rad*u.deg, catalog = "TIC")
                obs_df = self.tab2df(obs_table)
                if len(obs_table['ID']) == 1:
                    kic = obs_table['KIC'][0]
                    if kic == '--':
                        kic_found = False
                        continue
                    else:
                        kic_found = True
                        continue
                if len(obs_df[obs_df['GAIA'].to_numpy(dtype = 'str') != '']) == 1:
                    temp_obs_df = obs_df[obs_df['GAIA'].to_numpy(dtype = 'str') != '']
                    kic = temp_obs_df['KIC'].iloc[0]
                    if kic == '--':
                        kic_found = False
                        continue
                    else:
                        kic_found = True
                        continue
                if len(np.unique(obs_df[obs_df['HIP'].to_numpy(dtype = 'str') != '']['HIP'])) == 1:
                    kic = obs_table['KIC'][0]
                    if kic == '--':
                        kic_found = False
                        continue
                    else:
                        kic_found = True
                        continue
        if kic_found == False:
            kic = np.nan
            print("Didn't find KIC for this object.")
        else:
            print("Found KIC " + str(kic) + "!")
        return(kic)
    
    #most? TICs have an associated GAIA. possibly remove catalog filter to get more gaia
    def get_gaia_id(self,row):
        radii = np.linspace(start = 0.0001, stop = 0.001, num = 19)
        #for i,row in self.df.iterrows():
        gaia_id_found = False
        for rad in radii:
            if gaia_id_found == False:
                query_string = str(row['ra']) + " " + str(row['dec']) # make sure to have a space between the strings!#SkyCoord(ra = row['ra'], dec = row['dec'], frame = 'icrs') str(row['ra']) + " " + str(row['dec']) # make sure to have a space between the strings!
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
        return(gaia_id)
    
    def get_gaia_data(self,row,gaia_kwrgs):
        print("Working on GAIA data.")
        coord = SkyCoord(ra=row['ra']*u.degree,dec=row['dec']*u.degree)
        obs_table_gaia = Gaia.query_object(coord,radius = 0.002*u.degree)
        obs_table_gaia = obs_table_gaia[obs_table_gaia['source_id'] == int(row['gaia'])]
        
        if len(obs_table_gaia) == 1: #query_list1.append(obs_table_gaia.to_pandas())
            print("Found GAIA data!")
            if gaia_kwrgs == 'all':
                return(obs_table_gaia.to_pandas().squeeze())
            else:
                obs_pd_series = obs_table_gaia[gaia_kwrgs].to_pandas().squeeze()
                return(obs_pd_series)
        else:
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
                for i,col in enumerate(all_cols):
                    if i == 0:
                        gaia_series = pd.Series()
                    gaia_series[col] = np.nan
            
            else:
                for i,col in enumerate(gaia_kwrgs):
                    if i == 0:
                        gaia_series = pd.Series()
                    gaia_series[col] = np.nan
            
            return(gaia_series)
            
    def tab2df(self,table):    
        df = pd.DataFrame()
        for col in table.colnames: df[col] = table[col]
        return(df)    

    def roll_columns(self,df,num_roll):
        cols_rolled = np.roll(df.columns.to_numpy(dtype = 'str'),shift = num_roll)
        return(df[cols_rolled])    
        
        