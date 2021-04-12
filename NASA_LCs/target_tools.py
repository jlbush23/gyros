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

import matplotlib.pyplot as plt

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
               vmax = np.nanpercentile(image, 92),
               vmin = np.nanpercentile(image, 5))
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