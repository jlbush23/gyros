# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 19:31:41 2021

@author: jlbus
"""
import pandas as pd
import numpy as np

from astropy.timeseries import LombScargle

from astropy.table import Table

import exoplanet as xo

import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.gridspec as gridspec

from scipy.stats import binned_statistic as bin_stat


def my_LS_multi_sector(lc_df,flux_type,flux_err_avail = True, min_freq = 1/30):
    periodogram_list = []
    LS_results_list = []
    
    sector_list = np.unique(lc_df['sector'].to_numpy())
    
    for sector in sector_list:
        time = lc_df[lc_df['sector'] == sector]['time']
        flux = lc_df[lc_df['sector'] == sector][flux_type]
        if flux_err_avail == True:
            flux_err = lc_df[lc_df['sector'] == sector][flux_type + '_err']
        else:
            flux_err = None
        
        temp_LS_results_df,temp_periodogram_df = my_LS(time = time, flux = flux, flux_err = flux_err, max_per = 1/min_freq)
        
        temp_sector = np.repeat(a=str(sector),repeats = len(temp_LS_results_df))
        temp_LS_results_df['sector'] = temp_sector
        
        temp_sector = np.repeat(a=str(sector),repeats = len(temp_periodogram_df))
        temp_periodogram_df['sector'] = temp_sector
        
        temp_LS_results_df['sector'] = temp_LS_results_df['sector'].to_numpy(dtype = 'str')
    
        periodogram_list.append(temp_periodogram_df)
        LS_results_list.append(temp_LS_results_df)
        
    try:
        periodogram_LS = pd.concat(objs = periodogram_list, ignore_index = True)
        LS_results = pd.concat(objs = LS_results_list, ignore_index = True)
        periodogram_LS_available = True
        LS_results_available = True
        return(LS_results,periodogram_LS)
    except: 
        #might need to add a return statement down here
        periodogram_LS_available = False
        LS_results_available = False  
        LS_results = periodogram_LS = None
        return(LS_results,periodogram_LS)

def my_LS(time,flux,flux_err = None, max_per = 30): 
    min_freq = 1/max_per          
    
    if flux_err is None:                
        freq,power = LombScargle(time,flux).autopower(minimum_frequency = min_freq, maximum_frequency = 10.0)
    else:
        freq,power = LombScargle(time,flux,flux_err).autopower(minimum_frequency = min_freq, maximum_frequency = 10.0)
    
    threshhold = 0.005 
    peak_locs = []
    pow_peaks = []
    per_peaks = [] 
    
    for i in range(len(power)-2):
        if (power[i] > power[i+1]) and (power[i]>power[i-1]) and (power[i]>threshhold):
            if (power[i] > power[i-2]) and (power[i] > power[i+2]):
                peak_locs.append(i)
                pow_peaks.append(power[i])
                per_peaks.append((1/freq)[i])
            
    #take top 3 periods/powers
    sorted_pow_peaks = np.sort(pow_peaks)[::-1]
    top3 = sorted_pow_peaks[0:3]
    per_top3 = []
    for temp_power in top3:
        per_top3.append((per_peaks)[np.where(pow_peaks == temp_power)[0][0]])
        
    period = 1/freq
      
    if len(top3) == 0:
        Per1 = np.nan
        Power1 = np.nan
        Per2 = np.nan
        Power2 = np.nan
        Per3 = np.nan
        Power3 = np.nan
    
    if len(top3) == 1:
        Per1 = per_top3[0]
        Power1 = top3[0]
        Per2 = np.nan
        Power2 = np.nan
        Per3 = np.nan
        Power3 = np.nan
    
    if len(top3) == 2:
        if per_top3[0] != max_per:
            Per1 = per_top3[0]
            Power1 = top3[0]
            Per2 = per_top3[1]
            Power2 = top3[1]
            Per3 = np.nan
            Power3 = np.nan
        else:
            Per1 = per_top3[1]
            Power1 = top3[1]
            Per2 = np.nan
            Power2 = np.nan
            Per3 = np.nan
            Power3 = np.nan
                
    if len(top3) == 3:
        if per_top3[0] != max_per:
            Per1 = per_top3[0]
            Power1 = top3[0]
            Per2 = per_top3[1]
            Power2 = top3[1]
            Per3 = per_top3[2]
            Power3 = top3[2]
        else:
            Per1 = per_top3[1]
            Power1 = top3[1]
            Per2 = per_top3[2]
            Power2 = top3[2]
            Per3 = np.nan
            Power3 = np.nan
        
        
    temp_periodogram_df = (pd.DataFrame(data = np.transpose(np.array([period,power])), columns = ['period','power'], dtype = 'float')).sort_values(by = 'period', ascending = True)
    
    temp_LS_results_tab = Table([[Per1],[Per2],[Per3],[Power1],[Power2],[Power3]], 
                                names = ['LS_Per1','LS_Per2','LS_Per3','LS_Power1','LS_Power2','LS_Power3'])
    temp_LS_results_df = temp_LS_results_tab.to_pandas()
    return(temp_LS_results_df,temp_periodogram_df)


def exo_acf_multi_sector(lc_df, flux_type, flux_err_avail = False, max_per = 29):
    
    sector_list = np.unique(lc_df['sector'].to_numpy())
    
    result_list = []
    acf_periodograms = {}
    
    for sector in sector_list:
        temp_lc_df = lc_df[lc_df['sector'] == sector]
        
        time = temp_lc_df['time'].to_numpy(dtype = 'float')
        flux = temp_lc_df[flux_type].to_numpy(dtype = 'float')
        if flux_err_avail == True:
            flux_err = temp_lc_df[flux_type + '_err'].to_numpy(dtype = 'float')
        else:
            flux_err = None
        
        ac_period,temp_periodogram = exo_acf(time = time, flux = flux, flux_err = flux_err, max_per = 29)
        temp_result = pd.DataFrame(data = {'sector':[str(sector)],'ac_period':[ac_period]}, 
                          columns = ['sector','ac_period'])
        
        temp_sector = np.repeat(a=str(sector),repeats = len(temp_periodogram))
        temp_periodogram['sector'] = temp_sector
        
        acf_periodograms[sector] = temp_periodogram
        result_list.append(temp_result)
        
    acf_result = pd.concat(result_list)
    #id_repeats = np.repeat(self.id_num, len(acf_result))
    #acf_result.insert(loc = 0, column = id_label, value = id_repeats)
    
    return(acf_result,acf_periodograms)
        
def exo_acf(time,flux,flux_err = None, max_per = 29):        
    
    mu = np.mean(flux)
    flux = (flux / mu - 1)
    if flux_err is not None: flux_err = flux_err / mu
    
    try:
        #AC Results
        if flux_err is not None:
            ac_results = xo.autocorr_estimator(x=time, y=flux, yerr=flux_err, min_period=0.1, max_period=max_per)
        else:
            ac_results = xo.autocorr_estimator(x=time, y=flux, min_period=0.1, max_period=max_per)
        
        ac_peak1 = ac_results['peaks'][0]
        ac_period = ac_peak1['period']
        
        temp_periodogram = pd.DataFrame(data = np.transpose(ac_results['autocorr']), columns = ['Period','AC Power'])
        #ac_error = ac_peak1['period_uncert']
        
    except:
        ac_period = np.nan
        temp_periodogram = pd.DataFrame(columns = ['Period','AC Power'])
    
    # temp_result = pd.DataFrame(data = {'sector':str(sector),'ls_period1':[ls_period],'ls_logpower':[ls_power],
    #                               'ls_error':[ls_error],'ac_period':[ac_period]}, 
    #                       columns = ['sector','ls_period1','ls_logpower','ls_error','ac_period'])
    
    return(ac_period,temp_periodogram)


def amp_multi_sector(lc_df,flux_type):
    
    sector_list = np.unique(lc_df['sector'].to_numpy())
    
    for i,sector in enumerate(sector_list):        
        temp_lc = lc_df[lc_df['sector'] == sector]
        temp_flux = temp_lc[flux_type]
        amp = measure_amp(flux = temp_flux)
        if i == 0: 
            amp_df = pd.DataFrame(data = {'amp':[amp],'sector':[sector]})
        else:
            temp_df = pd.DataFrame(data = {'amp':[amp],'sector':[sector]})
            amp_df= pd.concat([amp_df,temp_df])      
    
    return(amp_df)

def measure_amp(flux):
    amplitude = np.percentile(a = flux, q = 95) - np.percentile(a = flux, q = 5)
    return(amplitude)

def best_period(rot_df):
    temp_Power1 = (rot_df['LS_Power1'].dropna()).to_numpy(dtype='float')        
    if len(temp_Power1) > 0:        
        row_idx = np.where(temp_Power1 == np.max(temp_Power1))[0][0]            
        best_df = rot_df.iloc[row_idx,:]
        best_per = best_df['LS_Per1']
        best_pow = best_df['LS_Power1']
        best_sect = str(best_df['sector'])
        return((best_per,best_pow,best_sect))
    else:
        return((np.nan,np.nan,np.nan))
    

def period_graph(target_name, lc_df, flux_type, LS_res, LS_periodogram, AC_res, AC_periodogram):
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
    
    # AC_power = self.periodogram_AC['power']
    
    # AC_height = np.max(AC_power) - np.min(AC_power)
    
    #make plots

    fig = plt.figure(figsize=(15,12))
    gridspec.GridSpec(3,2)
    
    ### PLOT FULL LC
    plt.subplot2grid((3,2), (0,0), colspan = 2)
    
    sector_list = np.unique(lc_df['sector'].to_numpy())
    for i,sector in enumerate(sector_list):
        
        time = lc_df[lc_df['sector'] == sector]['time']
        flux = lc_df[lc_df['sector'] == sector][flux_type]
        
        bin_flux,bin_time,bin_idx = bin_stat(x=time,values=flux,statistic = 'median', bins = round(0.05*len(time)))
        bin_time = bin_time[0:len(bin_time)-1]
              
        if flux_type != 'cpm':
            bin_flux = bin_flux/np.nanmedian(bin_flux)
            flux = flux/np.median(flux) 
            plt.scatter(time,flux, s = 0.75, label = sector)
            plt.plot(bin_time,bin_flux, c = 'black', linewidth = 1)#, c = c_sector[i])
        if flux_type == 'cpm':
            plt.scatter(time,flux, s = 1, label = sector)#, c = c_sector[i])
            #plt.plot(bin_time,bin_flux, c = 'black', linewidth = 1)
    
    if flux_type == 'cpm':
        plt.ylim((-0.15,0.15))
    plt.xlabel("Time (days)")
    plt.ylabel("Normalized Flux")
    plt.legend(loc = 'upper right', fontsize = 'xx-large')
    plt.title(str(flux_type) + "Light Curve of " + str(target_name))
    
    ### PLOT LS periodograms

    plt.subplot2grid((3,2), (1,0),colspan=1)
    
    for i,sector in enumerate(sector_list):
        LS_results = LS_res[LS_res['sector'] == sector]
        
        period = LS_periodogram[LS_periodogram['sector'] == sector]['period']
        power = LS_periodogram[LS_periodogram['sector'] == sector]['power']
        
#            ymax1 = LS_results['LS_Power1'].to_numpy()[0]/(1.05*LS_results['LS_Power1'].to_numpy()[0])
#            ymax2 = LS_results['LS_Power2'].to_numpy()[0]/(1.05*LS_results['LS_Power1'].to_numpy()[0])
#            ymax3 = LS_results['LS_Power3'].to_numpy()[0]/(1.05*LS_results['LS_Power1'].to_numpy()[0])
#            
        plt.plot(period,power)#,c = c_sector[i])
        plt.scatter(LS_results['LS_Per1'],LS_results['LS_Power1'], c = 'r')
        plt.scatter(LS_results['LS_Per2'],LS_results['LS_Power2'], c = 'r')
        plt.scatter(LS_results['LS_Per3'],LS_results['LS_Power3'], c = 'r')
#            plt.axvline(x=LS_results['LS_Per1'].to_numpy()[0], ymin=0.01,ymax = ymax1, c= 'r', linewidth = 4)
#            plt.axvline(x=LS_results['LS_Per2'].to_numpy()[0], ymin=0.01,ymax = ymax2, c= 'r', linewidth = 2)
#            plt.axvline(x=LS_results['LS_Per3'].to_numpy()[0], ymin=0.01,ymax = ymax3, c= 'r', linewidth = 1)
        
        
    temp_Power1 = (LS_res['LS_Power1'].dropna()).to_numpy(dtype='float')

    if len(temp_Power1) > 0:        
        row_idx = np.where(temp_Power1 == np.max(temp_Power1))[0][0]            
        best_per = LS_res['LS_Per1'].iloc[row_idx]            
        best_sector = LS_res['sector'].iloc[row_idx]
    else:
        best_per = best_sector = row_idx = np.nan
    
    plt.xlim((0,15))
    #plt.ylim((0,1))
    plt.xlabel("Period (days)")
    plt.ylabel("LS Power")
    plt.title("LS Period = " + str(round(best_per,4)) + " - Found in Sector" + str(best_sector))

    ### PLOT AC PERIODOGRAMS
    
    plt.subplot2grid((3,2), (1,1),colspan=1)
    
    for key in AC_periodogram:
        # nan_test = sum(np.isnan(self.lc_df[self.lc_df['sector'] == sector][self.flux_type].to_numpy(dtype = 'float')))
        # flux_len = len(self.lc_df[self.lc_df['sector'] == sector])
        # if nan_test != flux_len:
        AC_results = AC_res[AC_res['sector'] == str(key)]
        period = AC_periodogram[key]['Period']
        power = AC_periodogram[key]['AC Power']
        #print(AC_results)
        plt.plot(period,power)#, c = c_sector[i])
#             plt.scatter(AC_results['AC_Per1'],AC_results['AC_Power1'], c = 'r')
#             plt.scatter(AC_results['AC_Per2'],AC_results['AC_Power2'], c = 'r')
#             plt.scatter(AC_results['AC_Per3'],AC_results['AC_Power3'], c = 'r')
        plt.axvline(x=AC_results['ac_period'][0],color = 'red', ymin=0.01, ymax = 0.99, linewidth = 4)
#            plt.axvline(x=AC_results['AC_Per2'],color = 'red', ymin=0.01, ymax = (AC_results['AC_Power2'] - np.min(AC_power))/(AC_height + 0.05*AC_height), linewidth = 2)
#            plt.axvline(x=AC_results['AC_Per3'],color = 'red', ymin=0.01, ymax = (AC_results['AC_Power3'] - np.min(AC_power))/(AC_height + 0.05*AC_height), linewidth = 1)
    
    # ac periodogram title
    ac_title = ''
    for i,key in enumerate(AC_periodogram):
        AC_results = AC_res[AC_res['sector'] == str(key)]
        if i>0: ac_title = ac_title + ', '
        ac_title = ac_title + 'Sector ' + str(key) + ' : ' + str(round(AC_results['ac_period'][0],2)) + ' d'
    
    if len(AC_res['ac_period']) >= 1: 
        assoc_AC_per = AC_res[AC_res['sector'] == best_sector]['ac_period'][0]
    else:
        assoc_AC_per = np.nan
    # if np.isnan(row_idx) == False:
    #     assoc_AC_per = self.AC_results['AC_Per1'].iloc[row_idx]
    # else:
    #     assoc_AC_per = np.nan
    
    plt.xlim((0,15))
    #plt.ylim((0,1))
    plt.xlabel("Period (days)")
    plt.ylabel("AC Power")
    plt.title(ac_title)
    
    ### PLOT PHASE PHOLDED CURVES FOR BEST SECTOR
    
    time_best_sector = lc_df[lc_df['sector'] == best_sector]['time'].to_numpy(dtype = 'float')
    flux_best_sector = lc_df[lc_df['sector'] == best_sector][flux_type].to_numpy(dtype = 'float')
    
    bin_flux,bin_time,bin_idx = bin_stat(x=time_best_sector,values=flux_best_sector,statistic = 'median', bins = round(0.05*len(time)))
    bin_time = bin_time[0:len(bin_time)-1]
    bin_flux = bin_flux/np.nanmedian(bin_flux)
    flux = flux/np.median(flux)    
    
    #normalized phases arrays
   
    if (flux_type == 'cpm') | (flux_type == 'fcor'):
        phase_norm_LS = (time_best_sector % best_per)/best_per #from LS method 
        phase_norm_AC = (time_best_sector % assoc_AC_per)/assoc_AC_per #from AC
    else:
        phase_norm_LS = (bin_time % best_per)/best_per #from LS method        
        phase_norm_AC = (bin_time % assoc_AC_per)/assoc_AC_per #from AC
    #detect changes in phase

    c1 = []
    for i in range(len(phase_norm_LS)-1):
        if phase_norm_LS[i+1] - phase_norm_LS[i] < 0:
            c1.append(i+1)

    c2 = []
    for i in range(len(phase_norm_AC)-1):
        if phase_norm_AC[i+1] - phase_norm_AC[i] < 0:
            c2.append(i+1)
    
    # LS phase fold

    plt.subplot2grid((3,2), (2,0),colspan=1)
    first=0
    if (flux_type == 'cpm') | (flux_type == 'fcor'):
        for change in c1:
            plt.scatter(phase_norm_LS[first:change],flux_best_sector[first:change])
            first=change+1
    else:
        for change in c1:
            plt.scatter(phase_norm_LS[first:change],bin_flux[first:change])
            first=change+1
    plt.xlabel("Fraction of Period")
    plt.ylabel("Normalized Flux")
    plt.title("LS Phase-Folded - Sector " + str(best_sector))
    
    # AC phase fold

    plt.subplot2grid((3,2), (2,1),colspan=1)
    first = 0
    if (flux_type == 'cpm') | (flux_type == 'fcor'):
        for change in c2:
            plt.scatter(phase_norm_AC[first:change],flux_best_sector[first:change])
            first=change+1
    else:
        for change in c2:
            plt.scatter(phase_norm_AC[first:change],bin_flux[first:change])
            first=change+1        
    plt.xlabel("Fraction of Period")
    plt.ylabel("Normalized Flux")
    plt.title("AC Phase-Folded")
    
    fig.tight_layout()
    plt.close(fig=fig)
    
    
    return(fig)
        
        
    