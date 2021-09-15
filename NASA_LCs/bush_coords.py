# -*- coding: utf-8 -*-
"""
Created on Wed Sep 15 12:27:08 2021

@author: jlbus
"""
import numpy as np
import pandas as pd

import astropy.units as u
# from astroquery.gaia import Gaia
from astropy.coordinates import SkyCoord

import galpy.util.bovy_coords as bc #bc for bovy coords

def comove_coords(t,lit_gaia):
    ###could add other outputs like Vr, pred, in addition to sep,sep3d,and Vtan off
    
    
    Pcoord = t.sc
    
    
    # # Query Gaia with search radius and parallax cut
    # # Note, a cut on parallax_error was added because searches at low galactic latitude 
    # # return an overwhelming number of noisy sources that scatter into the search volume - ALK 20210325
    # print('Querying Gaia for neighbors')
    # if (searchradpc < Pcoord.distance):
    #     sqltext = "SELECT * FROM gaiaedr3.gaia_source WHERE CONTAINS( \
    #         POINT('ICRS',gaiaedr3.gaia_source.ra,gaiaedr3.gaia_source.dec), \
    #         CIRCLE('ICRS'," + str(Pcoord.ra.value) +","+ str(Pcoord.dec.value) +","+ str(searchraddeg.value) +"))\
    #         =1 AND parallax>" + str(minpar.value) + " AND parallax_error<0.5;"
    # if (searchradpc >= Pcoord.distance):
    #     sqltext = "SELECT * FROM gaiaedr3.gaia_source WHERE parallax>" + str(minpar.value) + " AND parallax_error<0.5;"
    #     print('Note, using all-sky search')
    # if verbose == True:
    #     print(sqltext)
    #     print()
    
    # job = Gaia.launch_job_async(sqltext , dump_to_file=False)
    # r = job.get_results()
       
    # if verbose == True: print('Number of records: ',len(r['ra']))
    
    
    # # Construct coordinates array for all stars returned in cone search
    
    # gaiacoord = SkyCoord( ra=r['ra'] , dec=r['dec'] , distance=(1000.0/r['parallax'])*u.parsec , \
    #                      frame='icrs' , \
    #                      pm_ra_cosdec=r['pmra'] , pm_dec=r['pmdec'] )
        
    lit_sc = SkyCoord(ra = lit_gaia.ra.to_numpy(dtype = 'float')*u.deg, 
                 dec = lit_gaia.dec.to_numpy(dtype = 'float')*u.deg,
                 pm_ra_cosdec = lit_gaia.pmra.to_numpy(dtype = 'float')*u.mas / u.yr,
                 pm_dec = lit_gaia.pmdec.to_numpy(dtype = 'float')*u.mas / u.yr,
                 distance = u.pc* (1000./lit_gaia.parallax.to_numpy(dtype = 'float')))
    
    sep = lit_sc.separation(Pcoord)#in degrees
    sep3d = lit_sc.separation_3d(Pcoord)#in parsec
    
    
    Pllbb     = bc.radec_to_lb(Pcoord.ra.value , Pcoord.dec.value , degree=True)
    Ppmllpmbb = bc.pmrapmdec_to_pmllpmbb( Pcoord.pm_ra_cosdec.value , Pcoord.pm_dec.value , \
                                         Pcoord.ra.value , Pcoord.dec.value , degree=True )
    Pvxvyvz   = bc.vrpmllpmbb_to_vxvyvz(Pcoord.radial_velocity.value , Ppmllpmbb[0] , Ppmllpmbb[1] , \
                                   Pllbb[0] , Pllbb[1] , Pcoord.distance.value/1000.0 , XYZ=False , degree=True)
    

    Gllbb = bc.radec_to_lb(lit_sc.ra.value , lit_sc.dec.value , degree=True)
    Gxyz = bc.lbd_to_XYZ( Gllbb[:,0] , Gllbb[:,1] , lit_sc.distance/1000.0 , degree=True)
    Gvrpmllpmbb = bc.vxvyvz_to_vrpmllpmbb( \
                    Pvxvyvz[0]*np.ones(len(Gxyz[:,0])) , Pvxvyvz[1]*np.ones(len(Gxyz[:,1])) , Pvxvyvz[2]*np.ones(len(Gxyz[:,2])) , \
                    Gxyz[:,0] , Gxyz[:,1] , Gxyz[:,2] , XYZ=True)
    Gpmrapmdec = bc.pmllpmbb_to_pmrapmdec( Gvrpmllpmbb[:,1] , Gvrpmllpmbb[:,2] , Gllbb[:,0] , Gllbb[:,1] , degree=True)
    
    # Code in case I want to do chi^2 cuts someday
    Gvtanerr = 1.0 * np.ones(len(Gxyz[:,0]))
    Gpmerr = Gvtanerr * 206265000.0 * 3.154e7 / (lit_sc.distance.value * 3.086e13)
    
    
    Gchi2 = ( (Gpmrapmdec[:,0]-lit_sc.pm_ra_cosdec.value)**2 + (Gpmrapmdec[:,1]-lit_sc.pm_dec.value)**2 )**0.5
    vtanoff = Gchi2 / Gpmerr #this is reported Vtan,off(km/s)
    
    ##vr pred
    vr_pred = Gvrpmllpmbb[:,0] 
    
    #create results dataframe
    res = pd.DataFrame(data = {'tic':[lit_gaia.tic.to_numpy(dtype = 'str')],
                               'designation':[lit_gaia.designation.to_numpy(dtype = 'str')],
                               'ra':[lit_sc.ra.value],
                               'dec':[lit_sc.dec.value],
                               'sep2D(deg)':[sep.value],
                               'sep3D(pc)':[sep3d.value],
                               'Vtan,off(km/s)':[vtanoff],
                               'Vr,pred(km/s)':[vr_pred]
                               }
                       )
    
    return(res)
    
    
    
    
    
    
    