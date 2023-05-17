# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 19:31:41 2021

@author: jlbus
"""
import pandas as pd
import numpy as np

from astropy.timeseries import LombScargle

from astropy.table import Table
from astropy.time import Time

import exoplanet as xo

import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['font.sans-serif'] = "DejaVu Sans"
mpl.rcParams['font.family'] = "sans-serif"
mpl.rcParams['text.usetex'] = False
import matplotlib.gridspec as gridspec

from scipy.stats import binned_statistic as bin_stat

from gatspy import periodic

import starspot as ss

import lightkurve as lk

def ss_rots_multi_sector(lc_df, flux_type, flux_err_avail= False,
                         min_per = 0.1, max_per = 29):
    
    sector_list = np.unique(lc_df['sector'].to_numpy())
    
    result_list = []
    
    for sector in sector_list:
        temp_lc_df = lc_df[lc_df['sector'] == sector]
        temp_lc_df = temp_lc_df.dropna(subset = ['time',flux_type])
        if len(temp_lc_df) == 0: continue
        
        if type(temp_lc_df['time'].to_numpy()[0]) == np.datetime64:  
            time = Time(temp_lc_df['time']).jd
        else:
            time = temp_lc_df['time'].to_numpy(dtype = 'float')
        flux = temp_lc_df[flux_type].to_numpy(dtype = 'float')
        if flux_err_avail == True:
            flux_err = temp_lc_df[flux_type + '_err'].to_numpy(dtype = 'float')
        else:
            flux_err = None
        
        ss_ls, ss_acf = starspot_rots(time = time, flux = flux, flux_err = flux_err, 
                                      flux_type = flux_type, min_per = min_per, max_per = max_per)
        temp_result = pd.DataFrame(data = {'sector':[str(sector)],'ss_ls':[ss_ls],'ss_acf':[ss_acf]}, 
                                   columns = ['sector','ss_ls','ss_acf'])
        
        # temp_sector = np.repeat(a=str(sector),repeats = len(temp_periodogram))
        # temp_periodogram['sector'] = temp_sector
        
        # acf_periodograms[sector] = temp_periodogram
        result_list.append(temp_result)
        
    ss_result = pd.concat(result_list)
    #id_repeats = np.repeat(self.id_num, len(acf_result))
    #acf_result.insert(loc = 0, column = id_label, value = id_repeats)
    
    return(ss_result)#,acf_periodograms)    

def starspot_rots(time,flux,flux_err=None,flux_type = None, min_per=0.1,max_per=30):
    
    rotate = ss.RotationModel(time,flux,flux_err = flux_err)

    # Calculate the Lomb Scargle periodogram period (highest peak in the periodogram).
    lomb_scargle_period = rotate.ls_rotation(min_period = min_per,max_period=max_per)
    
    # Calculate the autocorrelation function (ACF) period (highest peak in the ACF).
    # This is for evenly sampled data only -- time between observations is 'interval'.
    acf_period = rotate.acf_rotation(interval=np.diff(time)[0])
    
    return(lomb_scargle_period,acf_period)

def gatspy_LS_multi_sector(lc_df, flux_type, flux_err_avail = False, 
                           min_per = 0.1, max_per = 29):
    
    sector_list = np.unique(lc_df['sector'].to_numpy())
    
    result_list = []
    
    for sector in sector_list:
        temp_lc_df = lc_df[lc_df['sector'] == sector]
        temp_lc_df = temp_lc_df.dropna(subset = ['time',flux_type])
        if len(temp_lc_df) == 0: continue
        
        if type(temp_lc_df['time'].to_numpy()[0]) == np.datetime64:  
            time = Time(temp_lc_df['time']).jd
        else:
            time = temp_lc_df['time'].to_numpy(dtype = 'float')
        flux = temp_lc_df[flux_type].to_numpy(dtype = 'float')
        if flux_err_avail == True:
            flux_err = temp_lc_df[flux_type + '_err'].to_numpy(dtype = 'float')
        else:
            flux_err = None
        
        gatspy_period = gatspy_LS(time = time, flux = flux, flux_err = flux_err, 
                                  flux_type = flux_type, min_per = min_per, max_per = max_per)
        temp_result = pd.DataFrame(data = {'sector':[str(sector)],'gatspy_per':[gatspy_period]}, 
                                   columns = ['sector','gatspy_per'])
        
        # temp_sector = np.repeat(a=str(sector),repeats = len(temp_periodogram))
        # temp_periodogram['sector'] = temp_sector
        
        # acf_periodograms[sector] = temp_periodogram
        result_list.append(temp_result)
        
    gatspy_result = pd.concat(result_list)
    #id_repeats = np.repeat(self.id_num, len(acf_result))
    #acf_result.insert(loc = 0, column = id_label, value = id_repeats)
    
    return(gatspy_result)#,acf_periodograms)

def gatspy_LS(time,flux,flux_type,flux_err = None, min_per = 0.1, max_per = 30):
    
    model = periodic.LombScargleFast(fit_period=True,center_data = False)
    model.optimizer.period_range = (min_per, max_per)
    model.fit(time, flux)

    rot = model.best_period;
    
    return rot

def my_LS_multi_sector(lc_df,flux_type,flux_err_avail = True, min_freq = 1/30):
    periodogram_list = []
    LS_results_list = []
    
    sector_list = np.unique(lc_df['sector'].to_numpy())
    
    for sector in sector_list:
        temp_lc_df = lc_df[lc_df['sector'] == sector]
        temp_lc_df = temp_lc_df.dropna(subset = ['time',flux_type])
        if len(temp_lc_df) == 0: continue
    
        if type(temp_lc_df['time'].to_numpy()[0]) == np.datetime64:  
            time = Time(temp_lc_df['time']).jd
        else:
            time = temp_lc_df['time'].to_numpy(dtype = 'float')
        
        flux = temp_lc_df[flux_type].to_numpy(dtype = 'float')
        if flux_err_avail == True:
            flux_err = temp_lc_df[flux_type + '_err'].to_numpy(dtype = 'float')
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

def my_LS(time,flux,flux_err = None, min_per = 0.1, max_per = 20,
          normalization = 'standard', samples_per_peak = 5,
          threshold = 0.005): 
    max_freq = 1/min_per
    min_freq = 1/max_per          
    
    if flux_err is None:  
        ls = LombScargle(time,flux)              
        freq,power = ls.autopower(minimum_frequency = min_freq, maximum_frequency = max_freq,
                                  samples_per_peak = samples_per_peak,
                                  normalization = normalization)
    else:
        ls = LombScargle(time,flux,flux_err)
        freq,power = ls.autopower(minimum_frequency = min_freq, maximum_frequency = max_freq,
                                  samples_per_peak = samples_per_peak,
                                  normalization = normalization)
    
    # threshhold = 0.005 
    peak_locs = []
    pow_peaks = []
    per_peaks = [] 
    
    for i in range(len(power)-2):
        if (power[i] > power[i+1]) and (power[i]>power[i-1]) and (power[i]>threshold):
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
        Per1 = Power1 = Per2 = Power2 = Per3 = Power3 =  np.nan
        false_prob1 = false_prob2 = false_prob3 = np.nan
    
    if len(top3) == 1:
        Per1 = per_top3[0]
        Power1 = top3[0]
        false_prob1 = ls.false_alarm_probability(Power1)
        Per2 = Power2 = Per3 = Power3 =  np.nan
        false_prob2 = false_prob3 = np.nan
 
    
    if len(top3) == 2:
        if per_top3[0] != max_per:
            Per1 = per_top3[0]
            Power1 = top3[0]
            false_prob1 = ls.false_alarm_probability(Power1)
            Per2 = per_top3[1]
            Power2 = top3[1]
            false_prob2 = ls.false_alarm_probability(Power2)
            Per3 = Power3 = false_prob3 = np.nan
        else:
            Per1 = per_top3[1]
            Power1 = top3[1]
            false_prob1 = ls.false_alarm_probability(Power1)
            Per2 = Power2 = Per3 = Power3 =  np.nan
            false_prob2 = false_prob3 = np.nan
                
    if len(top3) == 3:
        if per_top3[0] != max_per:
            Per1 = per_top3[0]
            Power1 = top3[0]
            false_prob1 = ls.false_alarm_probability(Power1)
            Per2 = per_top3[1]
            Power2 = top3[1]
            false_prob2 = ls.false_alarm_probability(Power2)
            Per3 = per_top3[2]
            Power3 = top3[2]
            false_prob3 = ls.false_alarm_probability(Power3)
        else:
            Per1 = per_top3[1]
            Power1 = top3[1]
            false_prob1 = ls.false_alarm_probability(Power1)
            Per2 = per_top3[2]
            Power2 = top3[2]
            false_prob2 = ls.false_alarm_probability(Power2)
            Per3 = Power3 =  false_prob3 = np.nan
        
        
    temp_periodogram_df = (pd.DataFrame(data = np.transpose(np.array([period,power])), columns = ['period','power'], dtype = 'float')).sort_values(by = 'period', ascending = True)
    
    temp_LS_results_tab = Table([[Per1],[Per2],[Per3],[false_prob1],[false_prob2],[false_prob3], 
                                 [Power1],[Power2],[Power3]], 
                                names = ['LS_Per1','LS_Per2','LS_Per3',
                                         'False_Prob1','False_Prob2','False_Prob3',
                                         'LS_Power1','LS_Power2','LS_Power3'])
    temp_LS_results_df = temp_LS_results_tab.to_pandas()
    return(temp_LS_results_df,temp_periodogram_df)


def exo_acf_multi_sector(lc_df, flux_type, flux_err_avail = False, max_per = 29):
    
    sector_list = np.unique(lc_df['sector'].to_numpy())
    
    result_list = []
    acf_periodograms = {}
    
    for sector in sector_list:
        temp_lc_df = lc_df[lc_df['sector'] == sector]
        temp_lc_df = temp_lc_df.dropna(subset = ['time',flux_type])
        if len(temp_lc_df) == 0: continue
        
        if type(temp_lc_df['time'].to_numpy()[0]) == np.datetime64:  
            time = Time(temp_lc_df['time']).jd
        else:
            time = temp_lc_df['time'].to_numpy(dtype = 'float')
        flux = temp_lc_df[flux_type].to_numpy(dtype = 'float')
        if flux_err_avail == True:
            flux_err = temp_lc_df[flux_type + '_err'].to_numpy(dtype = 'float')
        else:
            flux_err = None
        
        ac_period,temp_periodogram = exo_acf(time = time, flux = flux, flux_err = flux_err, 
                                             flux_type = flux_type, max_per = max_per)
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
        
def exo_acf(time,flux,flux_type, flux_err = None, max_per = 29): 
    if flux_err is not None: lc_df = pd.DataFrame(data = {'time':time,'flux':flux, 'flux_err':flux_err})    
    if flux_err is None: lc_df = pd.DataFrame(data = {'time':time,'flux':flux})
    lc_df = lc_df.dropna()
    # lc_df = lc_df.sort_values(by = 'time',axis = 0)
    # lc_df = lc_df.drop_duplicates(subset = ['time'])
    
    time = lc_df['time'].to_numpy(dtype = 'float')
    flux = lc_df['flux'].to_numpy(dtype = 'float')
    if flux_err is not None: flux_err = lc_df['flux_err'].to_numpy(dtype = 'float')
    
    #probably could skip all of the above and just use np.nanmean below, but whatever
    if flux_type != 'cpm':
        mu = np.mean(flux)
        flux = (flux / mu) - 1
        if flux_err is not None: flux_err = flux_err / mu
    if flux_type == 'sap_poly3_fit_div':
        flux = flux - 1
    
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
        if flux_type != 'cpm':
            temp_flux = temp_flux/np.nanmedian(temp_flux)
        amplitude = amp(flux = temp_flux)
        if i == 0: 
            amp_df = pd.DataFrame(data = {'amp':[amplitude],'sector':[sector]})
        else:
            temp_df = pd.DataFrame(data = {'amp':[amplitude],'sector':[sector]})
            amp_df= pd.concat([amp_df,temp_df])      
    
    return(amp_df)

def amp(flux, percent = 5):
    amplitude = np.nanpercentile(a = flux, q = 100-percent) - np.nanpercentile(a = flux, q = percent)
    return(amplitude)

def cdpp(time,flux,flux_err = None, 
         transit_duration=13, savgol_window=101, savgol_polyorder=2, sigma=5.0):
    lcobj = lk.LightCurve(time = time, flux = 1+flux, flux_err = flux_err
                          )
    cdpp = lcobj.estimate_cdpp(transit_duration = transit_duration,
                          savgol_window = savgol_window,
                          savgol_polyorder = savgol_polyorder,
                          sigma = sigma)
    return(cdpp)

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
    

def period_graph(target_name, lc_df, flux_type, LS_res, LS_periodogram, AC_res, AC_periodogram, target_df):
    #This stuff is everything, use it for any python plot to make it nicer.
    import matplotlib as mpl
    mpl.rcParams['font.sans-serif'] = "DejaVu Sans"
    mpl.rcParams['font.family'] = "sans-serif"
    mpl.rcParams['text.usetex'] = False
    mpl.rcParams['lines.linewidth'] =3
    mpl.rcParams['axes.linewidth'] = 2
    mpl.rcParams['xtick.major.width'] =2
    mpl.rcParams['ytick.major.width'] =2
    mpl.rcParams['xtick.minor.width'] =1.5
    mpl.rcParams['ytick.minor.width'] =1.5
    mpl.rcParams['ytick.labelsize'] = 10
    mpl.rcParams['xtick.labelsize'] = 10
    mpl.rcParams['axes.labelsize'] = 10
    #mpl.rcParams['legend.numpoints'] = 1
    mpl.rcParams['axes.labelweight']='semibold'
    mpl.rcParams['mathtext.fontset']='stix'
    mpl.rcParams['font.weight'] = 'semibold'
    mpl.rcParams['axes.titleweight']='semibold'
    mpl.rcParams['axes.titlesize']=10
    
    # AC_power = self.periodogram_AC['power']
    
    # AC_height = np.max(AC_power) - np.min(AC_power)
    
    #make plots
    sector_list = np.unique(lc_df['sector'].to_numpy())
    
    ### PLOT FULL LC
    fig = plt.figure(figsize = (7.5,6), constrained_layout = True)
    outer = gridspec.GridSpec(2, 1,figure = fig,height_ratios = [1,3])
    spec1 = gridspec.GridSpecFromSubplotSpec(1,len(sector_list),subplot_spec = outer[0],wspace=0.1)
    spec2 = gridspec.GridSpecFromSubplotSpec(3,2,subplot_spec = outer[1], hspace = 0.2)
    
    axs = []
    for i in range(len(sector_list)):
        if i == 0: axs.append(fig.add_subplot(spec1[0,i]))
        if i>0: axs.append(fig.add_subplot(spec1[0,i],sharey=axs[i-1])) 
    r2ax1 = fig.add_subplot(spec2[0,0])
    r2ax2 = fig.add_subplot(spec2[0,1], frame_on = False, xticks = [], yticks = [])
    row2ax1 = fig.add_subplot(spec2[1,0])   
    row2ax2 = fig.add_subplot(spec2[1,1])
    row3ax1 = fig.add_subplot(spec2[2,0])
    row3ax2 = fig.add_subplot(spec2[2,1])
    #plt.tight_layout()
    
    d = .015

    #plot all sectors on all axes
    flux_max = []
    flux_min = []
    for j,ax in enumerate(axs):
        #for lc in lc_dict
        for i,sector in enumerate(sector_list):  
            temp_lc_df = lc_df[lc_df['sector'] == sector]
            if type(temp_lc_df['time'].to_numpy()[0]) == np.datetime64:  
                time = Time(temp_lc_df['time']).jd
            else:
                time = temp_lc_df['time'].to_numpy(dtype = 'float')
            flux = temp_lc_df[flux_type].to_numpy(dtype = 'float')
            
            bin_flux,bin_time,bin_idx = bin_stat(x=time,values=flux,statistic = 'median', bins = round(0.075*len(time)))
            bin_time = bin_time[1:len(bin_time)]
                  
            if flux_type != 'cpm':
                bin_flux = bin_flux/np.nanmedian(bin_flux)
                flux = flux/np.nanmedian(flux)
                ax.scatter(time,flux, s = 0.75, label = 'Sector ' + str(sector))
                ax.plot(bin_time,bin_flux, c = 'black', linewidth = 1)#, c = c_sector[i])
                
                if j == 0:
                    flux_max.append(np.percentile(a = flux, q = 98))
                    flux_min.append(np.percentile(a = flux, q = 2))
                
            if flux_type == 'cpm':
                ax.scatter(time,flux, s = 0.3, label = sector)#, c = c_sector[i])
                ax.plot(bin_time,bin_flux, c = 'black', linewidth = 0.75)
            if j== int(round(len(sector_list)/2)) - 1:    
                if flux_type == 'sap_flux':
                    ax.set_title("SAP Light Curve of " + str(target_name), loc = 'left')
                if flux_type == 'pdcsap_flux':
                    ax.set_title("PDCSAP Light Curve of " + str(target_name))
                if flux_type == 'cpm':
                    ax.set_title("CPM Light Curve of " + str(target_name))
            if j == int(round(len(sector_list))) - 1:
                ax.legend(loc = 'upper right', fontsize = 'x-large')
                
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
        

    #plt.legend(loc = 'upper right', fontsize = 'x-large')    
    
    ### PLOT LS periodograms
    
    #row2ax1 = fig.add_subplot(spec2[1,0:int(round(len(sector_list)/2))])
    
    #plt.subplot2grid((3,2), (1,0),colspan=1)
    #fig.add_subplot(2,2,1)
    
    for i,sector in enumerate(sector_list):
        LS_results = LS_res[LS_res['sector'] == sector]
        
        period = LS_periodogram[LS_periodogram['sector'] == str(sector)]['period']
        power = LS_periodogram[LS_periodogram['sector'] == str(sector)]['power']
        
#            ymax1 = LS_results['LS_Power1'].to_numpy()[0]/(1.05*LS_results['LS_Power1'].to_numpy()[0])
#            ymax2 = LS_results['LS_Power2'].to_numpy()[0]/(1.05*LS_results['LS_Power1'].to_numpy()[0])
#            ymax3 = LS_results['LS_Power3'].to_numpy()[0]/(1.05*LS_results['LS_Power1'].to_numpy()[0])
#            
        row2ax1.plot(period,power,linewidth=0.5)#,c = c_sector[i])
        row2ax1.scatter(LS_results['LS_Per1'],LS_results['LS_Power1'], c = 'r',s=0.3)
        row2ax1.scatter(LS_results['LS_Per2'],LS_results['LS_Power2'], c = 'r',s=0.3)
        row2ax1.scatter(LS_results['LS_Per3'],LS_results['LS_Power3'], c = 'r',s=0.3)
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
    
    row2ax1.set_xlim((0,15))
    #plt.ylim((0,1))
    row2ax1.set_xlabel("Period (days)", loc = 'right')
    row2ax1.set_ylabel("LS Power")
    row2ax1.set_title("LS Period = " + str(round(best_per,4)) + " - Found in Sector" + str(best_sector))

    ### PLOT AC PERIODOGRAMS
    
    #plt.subplot2grid((3,2), (1,1),colspan=1)
    #fig.add_subplot(2,2,2)
    #row2ax2 = fig.add_subplot(spec2[1,int(round(len(sector_list)/2)):])
    
    
    for key in AC_periodogram:
        # nan_test = sum(np.isnan(self.lc_df[self.lc_df['sector'] == sector][self.flux_type].to_numpy(dtype = 'float')))
        # flux_len = len(self.lc_df[self.lc_df['sector'] == sector])
        # if nan_test != flux_len:
        AC_results = AC_res[AC_res['sector'] == str(key)]
        period = AC_periodogram[key]['Period']
        power = AC_periodogram[key]['AC Power']
        #print(AC_results)
        row2ax2.plot(period,power,linewidth = 0.5)#, c = c_sector[i])
#             plt.scatter(AC_results['AC_Per1'],AC_results['AC_Power1'], c = 'r')
#             plt.scatter(AC_results['AC_Per2'],AC_results['AC_Power2'], c = 'r')
#             plt.scatter(AC_results['AC_Per3'],AC_results['AC_Power3'], c = 'r')
        row2ax2.axvline(x=AC_results['ac_period'][0],color = 'red', ymin=0.01, ymax = 0.99, linewidth = 0.5)
#            plt.axvline(x=AC_results['AC_Per2'],color = 'red', ymin=0.01, ymax = (AC_results['AC_Power2'] - np.min(AC_power))/(AC_height + 0.05*AC_height), linewidth = 2)
#            plt.axvline(x=AC_results['AC_Per3'],color = 'red', ymin=0.01, ymax = (AC_results['AC_Power3'] - np.min(AC_power))/(AC_height + 0.05*AC_height), linewidth = 1)
    
    # ac periodogram title
    ac_title = ''
    for i,key in enumerate(AC_periodogram):
        AC_results = AC_res[AC_res['sector'] == str(key)]
        if i>0: ac_title = ac_title + ', '
        ac_title = ac_title + 'Sector ' + str(key) + ' : ' + str(round(AC_results['ac_period'][0],2)) + ' d'
    
    if len(AC_res['ac_period']) >= 1: 
        if str(best_sector) == 'nan':
            assoc_AC_per = AC_res['ac_period'].iloc[0]
        else:            
            assoc_AC_per = AC_res[AC_res['sector'] == best_sector]['ac_period'][0]
    else:
        assoc_AC_per = np.nan
    # if np.isnan(row_idx) == False:
    #     assoc_AC_per = self.AC_results['AC_Per1'].iloc[row_idx]
    # else:
    #     assoc_AC_per = np.nan
    
    ac_title = 'Associated AC Period: ' + str(round(assoc_AC_per,3)) + ' d'
    
    row2ax2.set_xlim((0,15))
    #plt.ylim((0,1))
    row2ax2.set_xlabel("Period (days)", loc = 'right')
    row2ax2.set_ylabel("AC Power")
    row2ax2.set_title(ac_title)
    
    ### PLOT PHASE PHOLDED CURVES FOR BEST SECTOR
    if str(best_sector) == 'nan': best_sector = sector_list[0]
    if flux_type != 'cpm':
        time_best_sector = lc_df[lc_df['sector'] == int(best_sector)]['time']
        flux_best_sector = lc_df[lc_df['sector'] == int(best_sector)][flux_type].to_numpy(dtype = 'float')
    else:
        time_best_sector = lc_df[lc_df['sector'] == str(best_sector)]['time']
        flux_best_sector = lc_df[lc_df['sector'] == str(best_sector)][flux_type].to_numpy(dtype = 'float')
    
    if type(time_best_sector.to_numpy()[0]) == np.datetime64: time_best_sector = Time(time_best_sector).jd

    bin_flux,bin_time,bin_idx = bin_stat(x=time_best_sector,values=flux_best_sector,
                                         statistic = 'median', bins = round(0.075*len(time_best_sector)))
    bin_time = bin_time[1:len(bin_time)]
    #bin_flux = bin_flux/np.nanmedian(bin_flux)
    #flux = flux/np.median(flux)    
    
    #normalized phases arrays
   
    if (flux_type == 'cpm') | (flux_type == 'fcor'):
        phase_norm_LS = (time_best_sector % best_per)/best_per #from LS method 
        phase_norm_LS = phase_norm_LS.to_numpy()
        phase_norm_AC = (time_best_sector % assoc_AC_per)/assoc_AC_per#from AC
        phase_norm_AC = phase_norm_AC.to_numpy() 
    else:
        phase_norm_LS = (bin_time % best_per)/best_per #from LS method  
        #phase_norm_LS = phase_norm_LS.to_numpy()
        phase_norm_AC = (bin_time % assoc_AC_per)/assoc_AC_per #from AC
        #phase_norm_AC = phase_norm_AC.to_numpy()
    #detect changes in phase
    
    c1 = [i+1 for i in range(len(phase_norm_LS) -1) if (phase_norm_LS[i+1] - phase_norm_LS[i]) < 0]
    c2 = [i+1 for i in range(len(phase_norm_AC) -1) if (phase_norm_AC[i+1] - phase_norm_AC[i]) < 0]
    # c1 = []
    # for i in range(len(phase_norm_LS)-1):
    #     if phase_norm_LS[i+1] - phase_norm_LS[i] < 0:
    #         c1.append(i+1)

    # c2 = []
    # for i in range(len(phase_norm_AC)-1):
    #     if phase_norm_AC[i+1] - phase_norm_AC[i] < 0:
    #         c2.append(i+1)
    
    # LS phase fold

    #plt.subplot2grid((3,2), (2,0),colspan=1)
    #fig.add_subplot(3,2,1)
    #row3ax1 = fig.add_subplot(spec[2,0:int(round(len(sector_list)/2))])
    
    first=0
    if (flux_type == 'cpm') | (flux_type == 'fcor'):
        for change in c1:
            row3ax1.scatter(phase_norm_LS[first:change],flux_best_sector[first:change],s=0.3)
            first=change+1
    else:
        for change in c1:
            row3ax1.scatter(phase_norm_LS[first:change],bin_flux[first:change],s=0.3)
            first=change+1
    row3ax1.set_xlabel("Fraction of Period")
    row3ax1.set_ylabel("Normalized Flux")
    row3ax1.set_title("LS Phase-Folded - Sector " + str(best_sector), loc = 'left')
    
    # AC phase fold

    #plt.subplot2grid((3,2), (2,1),colspan=1)
    #fig.add_subplot(3,2,2)
    #row3ax2 = fig.add_subplot(spec[2,int(round(len(sector_list)/2)):])
    
    first = 0
    if (flux_type == 'cpm') | (flux_type == 'fcor'):
        for change in c2:
            row3ax2.scatter(phase_norm_AC[first:change],flux_best_sector[first:change],s=0.3)
            first=change+1
    else:
        for change in c2:
            row3ax2.scatter(phase_norm_AC[first:change],bin_flux[first:change],s=0.3)
            first=change+1        
    row3ax2.set_xlabel("Fraction of Period")
    row3ax2.set_ylabel("Normalized Flux")
    row3ax2.set_title("AC Phase-Folded", loc = 'left')
    
    #fig.tight_layout()
    # plt.subplots_adjust(hspace = 0.2)
    # plt.close(fig=fig)
    
    ## plot best sector
    # bin_flux,bin_time,bin_idx = bin_stat(x=time_best_sector,values=flux_best_sector,
    #                                      statistic = 'median', bins = round(0.05*len(time_best_sector)))
    # bin_time = bin_time[0:len(bin_time)-1]
                  
    
    r2ax1.scatter(time_best_sector,flux_best_sector, c = 'black', alpha = 0.7, s = 0.3)
    r2ax1.plot(bin_time,bin_flux, c = 'orange', linewidth = 0.75)
    r2ax1.set(xlabel = 'Time (BJD -2457000)',
              ylabel = 'Normalized Flux', title = 'Best Sector Close-Up')
    
    ## plot available info
    r2ax2.annotate("TIC " + str(target_df['tic'][0]), xy = (-0.25,0.6), xycoords = "axes fraction",
             fontsize = 14)
    
    r2ax2.annotate("Tmag = " + str(target_df['Tmag'][0]), xy = (-0.25,0.35), xycoords = "axes fraction",
                 fontsize = 14)#, color = "green")
    
    ls_div_ac = best_per / assoc_AC_per
    
    r2ax2.annotate("LS / ACF = " + str(round(ls_div_ac,4)), xy = (-0.25,0.1), xycoords = "axes fraction",
                   fontsize = 14)

    ratio_type,ratio_color = ratio_res(rot_ratio = ls_div_ac)
    
    r2ax2.annotate("Match Status: " + str(ratio_type), xy = (-0.25,0.9), xycoords = "axes fraction",
                   fontsize = 16, color = ratio_color)
            
    
    #fig.tight_layout()
    plt.close(fig = fig)
    
    
    
    
    return(fig)

def ratio_res(rot_ratio, perc_err = 0.1):
    if (rot_ratio < 1 + perc_err) & (rot_ratio > 1 - perc_err):
        ratio_type = "Match"
        ratio_color = "green"
    elif (rot_ratio < 2 + perc_err) & (rot_ratio > 2 - perc_err):
        ratio_type = "0.5x Alias"
        ratio_color = "#B2A700"  
    elif (rot_ratio > 0.5 - perc_err) & (rot_ratio < 0.5 + perc_err):
        ratio_type = "2x Alias"
        ratio_color = "#B2A700"
    else:
        ratio_type = "No Match"
        ratio_color = "red"  
    return(ratio_type,ratio_color)
        
        
    