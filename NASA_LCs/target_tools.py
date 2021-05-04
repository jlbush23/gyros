# -*- coding: utf-8 -*-
"""
Created on Sun Apr 11 21:02:53 2021

@author: jlbus
"""

import pandas as pd
import numpy as np

import os

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS

from astroquery.gaia import Gaia

from scipy.stats import binned_statistic as bin_stat

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl
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

def target_contam_gaia(ra,dec,srad=15.5*20):
    c = SkyCoord(ra=ra*u.degree,dec=dec*u.degree)
    # Pgaia = Gaia.query_object_async(coordinate=c, radius=(5.0*u.arcsec))
    sqltext = "SELECT * FROM gaiaedr3.gaia_source WHERE CONTAINS( \
               POINT('ICRS',gaiaedr3.gaia_source.ra,gaiaedr3.gaia_source.dec), \
               CIRCLE('ICRS'," + str(c.ra.value) +","+ str(c.dec.value) +","+ str(srad/3600.0) +"))=1;"
    job = Gaia.launch_job_async(sqltext , dump_to_file=False)
    obs_table_gaia = job.get_results()
    return(obs_table_gaia)

def contamination_plot(fig,gaia_contam,median_im,im_header,target_df):
    target_df = target_df.squeeze()
    
    def plot_cutout(image):
        """
        Plot image and add grid lines.
        """
        plt.imshow(image, origin = 'lower', cmap = plt.cm.YlGnBu_r, 
               vmax = np.nanpercentile(image, 98),
               vmin = np.nanpercentile(image, 2))
        cbar = plt.colorbar(fraction = 0.04, orientation = 'horizontal', pad = 0.12)
        cbar.ax.get_yaxis().labelpad = 15
        cbar.ax.set_xlabel('Flux (e-/s)')#, rotation=180)
        plt.grid(axis = 'both',color = 'white', ls = 'solid')
        
    wcs = WCS(im_header)
    
    fig.add_subplot(111, projection = wcs)
    plot_cutout(median_im)
    
    plt.xlabel('RA', fontsize = 18)
    plt.ylabel('Dec', fontsize = 18)    
    
    starloc = wcs.all_world2pix([[target_df['ra'],target_df['dec']]],0)  #Second is origin
    plt.scatter(starloc[0,0], starloc[0,1],s = 200,c = target_df['phot_g_mean_mag'],
                cmap = 'autumn', vmin = np.nanmin(gaia_contam['phot_g_mean_mag']),
                vmax = np.nanmax(gaia_contam['phot_g_mean_mag']), marker='*')
    
    # Plot nearby stars as well, which we created using our Catalog call above.
    contam_loc_array = gaia_contam.to_pandas()[['ra','dec']].to_numpy(dtype = 'float')
    nearbyLoc = wcs.all_world2pix(contam_loc_array,0)
    
    plt.scatter(nearbyLoc[:, 0], nearbyLoc[:, 1], 
                s = 25, c = gaia_contam['phot_g_mean_mag'], cmap = 'autumn')
    
    cbar = plt.colorbar(fraction = 0.04)
    cbar.ax.get_yaxis().labelpad = 35
    cbar.ax.set_ylabel('$m_G$', rotation = 270, fontsize = 30)
    
    plt.title('TIC ' + str(target_df['tic']) + ": (Tmag = " + str(round(target_df['Tmag'],2)) + 
                           ", ContRatio = " + str(round(target_df['contratio'],5)) + ")")
    
    return(fig)

def create_target_lc_dict(target):
    
    lc_dict = {}

    if 'spoc_lc' in target.available_attributes: lc_dict['spoc_lc'] = target.spoc_lc
    if 'tpf_lc' in target.available_attributes: lc_dict['tpf_lc'] = target.tpf_lc
    if 'cpm_lc' in target.available_attributes: lc_dict['cpm'] = target.cpm_lc
    
    return(lc_dict)

def plot_LCs(fig, spec, target,
                 # target_name, lc_dict, rot_dict,
                 # lc_df, flux_type, LS_res, LS_periodogram, AC_res, AC_periodogram
                 ):    
    lc_dict = create_target_lc_dict(target)
    
    ##extraction methods and sector list
    lc_types = list(lc_dict.keys())
    flux_types = []
    for lc_type in lc_types:
        if lc_type == 'cpm': flux_types.append('cpm')
        if lc_type == 'spoc':
            flux_types.append('sap_flux')
            flux_types.append('pdcsap_flux')
    
    sector_list = np.unique(lc_dict[lc_types[0]]['sector'].to_numpy(dtype = 'str'))
    
    #unpack target info
    target_name = 'TIC ' + str(target.tic)
    if 'target_df' in target.available_attributes:
        tmag_info = 'Tmag = ' + str(target.target_df['Tmag'])
        contratio_info = 'ContRatio = ' + str(target.target_df['contratio'])
    else:
        tmag_info = None
    
    # fig = plt.figure(figsize = (15,12))
    # outer = gridspec.GridSpec(2, 1,figure = fig,height_ratios = [1,3])
    # spec1 = gridspec.GridSpecFromSubplotSpec(1,len(sector_list),subplot_spec = outer[0],wspace=0.1)
    # spec2 = gridspec.GridSpecFromSubplotSpec(2,2,subplot_spec = outer[1], hspace = 0.2)
    
    axs = []
    for i in range(len(sector_list)):
        if i == 0: axs.append(fig.add_subplot(spec[0,i]))
        if i>0: axs.append(fig.add_subplot(spec[0,i],sharey=axs[i-1])) 
    # row2ax1 = fig.add_subplot(spec2[0,0])   
    # row2ax2 = fig.add_subplot(spec2[0,1])
    # row3ax1 = fig.add_subplot(spec2[1,0])
    # row3ax2 = fig.add_subplot(spec2[1,1])
    plt.tight_layout()
    
    d = .015 #subplot spacing constant

    #plot all sectors on all axes
    flux_max = []
    flux_min = []
    color_list = ['tab:blue','tab:orange','tab:green','tab:red','tab:purple','tab:brown',
                  'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan','darkslateblue', 'gold', 'deepskyblue',
                  'navy', 'orangered','springgreen','deeppink','dimgrey','darggoldenrod','violet','cadetblue']
    # colors = ['b', 'r', 'g', 'c']
    # cc = itertools.cycle(colors)
    
    for j,ax in enumerate(axs):
        #for iterate over flux_types
        for flux_type in flux_types:
            if flux_type == 'cpm': 
                lc_df = lc_dict['cpm']
            if (flux_type == 'sap_flux') or (flux_type == 'pdcsap_flux'):
                lc_df = lc_dict['spoc_lc']
            # create flux label
            if flux_type == 'sap_flux': flux_label = 'SAP'
            if flux_type == 'pdcsap_flux': flux_label = 'PDC'
            if flux_type == 'cpm': flux_label = 'CPM'
            # if flux_type == 'tpf':
            
            lc_df['sector'] = lc_df['sector'].to_numpy(dtype = 'str')
            sector_lines = []
            flux_lines = []
            for i,sector in enumerate(sector_list):
                # get sector specific temp LC 
                temp_lc_df = lc_df[lc_df['sector'] == sector]
                ### get time - NEED TO UPDATE TIME FORMAT !!
                if(len(temp_lc_df) == 0):
                    continue
                elif (type(temp_lc_df['time'].to_numpy()[0]) == np.datetime64):  
                    time = Time(temp_lc_df['time']).jd 
                else:
                    time = temp_lc_df['time'].to_numpy(dtype = 'float')
                # get flux from temp LC with flux_type
                flux = temp_lc_df[flux_type].to_numpy(dtype = 'float')
                
                #create binned flux for plot line
                bin_flux,bin_time,bin_idx = bin_stat(x=time,values=flux,statistic = 'median', bins = round(0.05*len(time)))
                bin_time = bin_time[0:len(bin_time)-1]
                      
                if flux_type != 'cpm':
                    # normalize and plot spoc lightcurve
                    bin_flux = bin_flux/np.nanmedian(bin_flux)
                    flux = flux/np.nanmedian(flux)
                    sector_line = ax.scatter(time,flux, s = 0.75, label = 'Sector ' + str(sector), c = color_list[i], alpha = 0.7, marker = '^')
                    flux_line = ax.scatter(time,flux, s = 0.75, label = flux_label, c = color_list[i], alpha = 0.7, marker = '^')
                    ax.plot(bin_time,bin_flux, c = 'black', linewidth = 1)#, c = c_sector[i])
                    #append lines to list for legend creation
                    sector_lines.append(sector_line)
                    flux_lines.append(flux_line)
                    
                    if j == 0:
                        flux_max.append(np.percentile(a = flux, q = 98))
                        flux_min.append(np.percentile(a = flux, q = 2))
                    
                if flux_type == 'cpm':
                    #plot cpm lightcurve
                    sector_line = ax.scatter(time,flux+1, s = 1, label = 'Sector ' + str(sector), c = color_list[i], alpha = 0.9, marker = 'o')#, c = c_sector[i])
                    flux_line = ax.scatter(time,flux+1, s = 1, label = flux_label, c = color_list[i], alpha = 0.9, marker = 'o')
                    #plt.plot(bin_time,bin_flux, c = 'black', linewidth = 1)
                    sector_lines.append(sector_line)
                    flux_lines.append(flux_line)                    
                
                # set title
                if (j== int(round(len(sector_list)/2)) - 1) or (len(sector_list) == 1):
                    title = target_name
                    if tmag_info is not None:
                        title = title + ', ' + tmag_info + ', ' + contratio_info
                    ax.set_title(title, loc = 'left')
                
        # set legends
        if j == 0:
            first_legend = ax.legend(handles = flux_lines, loc = 'upper left', fontsize = 'x-large', framealpha = 0.7)
            ax.add_artist(first_legend)
        if (j+1) == len(axs):
            second_legend = ax.legend(handles = flux_lines, loc = 'upper right', fontsize = 'x-large', framealpha = 0.7)
            ax.add_artist(second_legend)
                    
    for i,(ax,sector) in enumerate(zip(axs,sector_list)):
        temp_lc_df = lc_df[lc_df['sector'] == sector]
        if type(temp_lc_df['time'].to_numpy()[0]) == np.datetime64:  
            time = Time(temp_lc_df['time']).jd
        else:
            time = temp_lc_df['time'].to_numpy(dtype = 'float')
        ax.set_xlim([np.min(time), np.max(time)])
        #ax.set_ylim([np.nanmax(flux_min),np.nanmin(flux_max)])
        #if flux_type == 'cpm': ax.set_ylim((-0.25,0.25))
        if i==0: ax.set_ylabel("Normalized Flux")
        if i== int(round(len(sector_list)/2)): 
            ax.set_xlabel("Time (JD)", loc = 'left')
        
        if (i+1) < len(axs) :
            axs[i].spines['right'].set_visible(False)
            axs[i+1].tick_params(left=False,right = False)
            axs[i+1].spines['left'].set_visible(False)        
            #axs[i+1].yaxis.tick_right()
            axs[i+1].tick_params(labelleft=False)  # don't put tick labels on the right
            #axs[i+1].set_yticks([])
            kwargs = dict(transform=axs[i].transAxes, color='k', clip_on=False)
            axs[i].plot((1 - d, 1 + d), (-d, +d), **kwargs)
            kwargs.update(transform=axs[i+1].transAxes)
            axs[i+1].plot((-d, +d), (-d, +d), **kwargs)