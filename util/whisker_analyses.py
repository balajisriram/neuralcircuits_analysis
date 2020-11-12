import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy.optimize import curve_fit

from .kwik_utils import get_unclustered_unit_details

def exponential_func(x, a, b, c):
    return a*np.exp(-b*x)+c

def get_r_squared(data,fits):
    sse = np.sum(np.power(data-fits,2))
    sst = np.sum(np.power(data-data.mean(),2))
    return 1-sse/sst

def analyze_mua_by_channel_multistim(loc,type='raster',show_plot=True,min_z=5,z_zoom=1,interval=[-1,3],stim_time='timestamps.np',common_details=None,chans=range(16)):
    # kwik_file
    kwik_file = [f for f in os.listdir(loc) if f.endswith('.kwik')]
    kwik_file = kwik_file[0]
    unit_details = get_unclustered_unit_details(os.path.join(loc,kwik_file))
    units = []
    # timestamps
    timestamps = np.fromfile(os.path.join(loc,stim_time))
    istimi = np.mean(np.diff(timestamps))
    
    # is it single stim or multi_stim
    if np.abs(istimi-14)<1 or np.abs(istimi-7)<1:
        n_stim=20
    elif np.abs(istimi-10)<1 or np.abs(istimi-5)<1:
        n_stim=1
    else:
        breakpoint()
    # get every other timestamp
    ratios = 1
    fig, ax = plt.subplots(16,4,figsize=(20,20))
    plt.subplots_adjust(wspace=0.1, hspace=0.1)

    cmap = matplotlib.cm.get_cmap('Dark2')
    max_z_lim = -np.inf
    max_summed_z_lim = -np.inf
        
    for i,unit in enumerate(unit_details['units']):
        spike_times = unit['spike_time']
        bin_interval = 0.02 # 20ms
        binned_time= np.linspace(interval[0],interval[1],np.int((interval[1]-interval[0])/bin_interval+1))
        firing_rate = np.empty((timestamps[0::ratios].size,binned_time.size-1))
        unit_color = cmap(float(i)/len(unit_details['units']))
        for j,t in enumerate(timestamps[0::ratios]):
            # get the relevant spike timestamps
            spike_that_trial = spike_times-t
            spike_that_trial = spike_that_trial[np.bitwise_and(spike_that_trial>=interval[0],spike_that_trial<=interval[1])]
            temp = np.histogram(spike_that_trial,binned_time)
            firing_rate[j,:] = temp[0]
            ax[15-i,0].vlines(spike_that_trial,ymin=j,ymax=j+0.8,color=unit_color)
            
        m = np.mean(firing_rate,axis=0)/bin_interval
        sd = np.std(firing_rate,axis=0)/bin_interval
        m_prestim = np.mean(np.mean(firing_rate[:,binned_time[:-1]<0],axis=0))/bin_interval
        sd_prestim = np.mean(np.std(firing_rate[:,binned_time[:-1]<0],axis=0))/bin_interval
        z_score = np.divide(m-m_prestim,sd_prestim)
        norm_score = np.divide(m-m_prestim,m_prestim)
        
        if np.max(z_score)>max_z_lim: max_z_lim = np.max(z_score)
        offset = 0
        ax[15-i,2].plot([binned_time[0],binned_time[-1]],[offset,offset],'k--')
        ax[15-i,3].plot([binned_time[0],binned_time[-1]],[offset,offset],'k--')
        print('chan for subject',chans)
        if np.max(z_score)>min_z and i in chans:
            this_unit = {}
            ax[15-i,1].plot(binned_time[:-1],m,color=unit_color,linewidth=3,alpha=0.5)
            ax[15-i,2].plot(binned_time[:-1],offset+z_zoom*z_score,color=unit_color,linewidth=3,alpha=0.5)
            max_z_for_unit = []
            summed_z_for_unit = []
            time_base = []
            # get sequential max_values
            for k in range(n_stim):
                that_time_bin = np.bitwise_and(binned_time[:-1]>=0.1*k,binned_time[:-1]<=0.1*(k+1))
                max_z = np.max(z_score[that_time_bin])
                max_z_ind = np.argmax(z_score[that_time_bin])
                total_z = np.sum(z_score[that_time_bin])
                if np.max(total_z)>max_summed_z_lim: max_summed_z_lim = np.max(total_z)
                ax[15-i,2].plot(binned_time[:-1][that_time_bin][max_z_ind],offset+z_zoom*max_z,color=unit_color,marker='*')
                ax[15-i,3].plot(0.1*k+0.05,offset+total_z,color=unit_color,marker='x')
                max_z_for_unit.append(max_z)
                summed_z_for_unit.append(total_z)
                time_base.append(0.1*k+0.05)
            this_unit['mean_activity_prestim'] = m_prestim
            this_unit['mean_activity_unnormalized'] = m
            this_unit['std_activity_unnormalized'] = sd
            this_unit['std_activity_prestim'] = sd_prestim
            this_unit['max_z'] = np.asarray(max_z_for_unit)
            this_unit['summed_z'] = np.asarray(summed_z_for_unit)
            if n_stim==20:
                this_unit['stable_response_ms'] = np.mean(this_unit['max_z'][-4:])/this_unit['max_z'][0]
                this_unit['first_pulse_ratio_ms'] = this_unit['max_z'][1]/this_unit['max_z'][0]
                this_unit['second_pulse_ratio_ms'] = this_unit['max_z'][2]/this_unit['max_z'][0]
            else:
                this_unit['stable_response_ms'] = np.nan
                this_unit['first_pulse_ratio_ms'] = np.nan
                this_unit['second_pulse_ratio_ms'] = np.nan
            this_unit['max_response'] = np.max(max_z_for_unit)
            this_unit['channel_number'] = i
            this_unit['z_score'] = z_score
            this_unit['norm_score'] = norm_score
            this_unit['istimi'] = np.mean(np.diff(timestamps))
            this_unit['n_stim'] = n_stim
            
            time_base = np.asarray(time_base)
            
            # look at off responses
            if n_stim==20:
                time_bin_off = np.bitwise_and(binned_time[:-1]>=2.0,binned_time[:-1]<=2.5) # 2-2.5 s
            else:
                time_bin_off = np.bitwise_and(binned_time[:-1]>=0.2,binned_time[:-1]<=1) # 0.2 - 1 sec
            total_off_z = np.sum(z_score[time_bin_off])
            this_unit['off_response'] = total_off_z
            
            # fit to exponential_func
            try:
                if n_stim==20:
                    max_z = this_unit['max_z']
                    popt, pcov = curve_fit(exponential_func, time_base, max_z, p0=(max_z[0], 0.5, 1),bounds=([max_z[0]*0.95,1./800,0],[max_z[0]*1.05,np.inf,np.inf]))
                    xx = np.linspace(0, 2, 1000)
                    yy = exponential_func(xx, *popt)
                    ax[15-i,2].plot(xx,yy+offset,color=unit_color,linewidth=2)
                    this_unit['max_z_t_ms'] = 1/popt[1]
                    #ax[i,1].text(xx[-1],offset,'{0:.2f}'.format(1000/popt[1]),horizontalalignment='center',verticalalignment='bottom',fontsize=20)
                    this_unit['max_z_t_quality_ms'] = get_r_squared(max_z,exponential_func(time_base, *popt))
                else: 
                    this_unit['max_z_t_ms'] = np.nan
                    this_unit['max_z_t_quality_ms'] = np.nan
            except RuntimeError:
                print('Runtime')
                this_unit['max_z_t_ms'] = np.nan
                this_unit['max_z_t_quality_ms'] = np.nan
            except ValueError: 
                this_unit['max_z_t_ms'] = np.nan
                this_unit['max_z_t_quality_ms'] = np.nan
            
            try:
                if n_stim==20:
                    summed_z = this_unit['summed_z']
                    popt, pcov = curve_fit(exponential_func, time_base, summed_z, p0=(summed_z[0], 0.5, 1),bounds=([summed_z[0]*0.95,1./800,0],[summed_z[0]*1.05,np.inf,np.inf]))
                    xx = np.linspace(0, 2, 1000)
                    yy = exponential_func(xx, *popt)
                    ax[15-i,3].plot(xx,yy+offset,color=unit_color,linewidth=2)
                    this_unit['total_z_t_ms'] = 1/popt[1]
                    ax[15-i,3].text(xx[-1],2,'{0:.1f}'.format(1000/popt[1]),horizontalalignment='center',verticalalignment='bottom',fontsize=20)
                    this_unit['total_z_t_quality_ms'] = get_r_squared(max_z,exponential_func(time_base, *popt))
                else:
                    this_unit['max_z_t_ms'] = np.nan
                    this_unit['max_z_t_quality_ms'] = np.nan
            except RuntimeError:
                print('Runtime')
                this_unit['total_z_t_ms'] = np.nan
                this_unit['total_z_t_quality_ms'] = np.nan
            except ValueError: 
                this_unit['total_z_t_ms'] = np.nan
                this_unit['total_z_t_quality_ms'] = np.nan
                
            # add the common details
            if common_details:
                this_unit['subject_id'] = common_details['subject_id']
                this_unit['date'] = common_details['date']
                this_unit['time'] = common_details['time']
                this_unit['genotype'] = common_details['genotype']
                this_unit['stim_location'] = common_details['stim_location']
                this_unit['unit_id'] = common_details['subject_id'] + '_' + common_details['genotype'] + '_' + str(i)
                this_unit['unit_and_tank_id'] = common_details['tank_name'] + '_' + str(i)
                this_unit['repeat'] = common_details['repeat']
                if 'drug_state' in common_details: this_unit['drug_state'] = common_details['drug_state']
                
            units.append(this_unit)
        else: ax[15-i,2].plot(binned_time[:-1],offset+z_zoom*z_score,color=(0.5,0.5,0.5))
        
        # get the ylim for ax[0]
        ylim0 = ax[15-i,1].get_ylim()
        ylim1 = ax[15-i,2].get_ylim()
        for j in range(n_stim):
            stim1 = matplotlib.patches.Rectangle((j*0.1,0.0),0.05,ylim0[1],color='blue',alpha=0.15)
            stim2 = matplotlib.patches.Rectangle((j*0.1,0.0),0.05,ylim1[1],color='blue',alpha=0.15)
            ax[15-i,1].add_patch(stim1)
            ax[15-i,2].add_patch(stim2)

    fig.suptitle(loc)
    ax[0,0].set_title('raster',fontsize=20)
    ax[0,2].set_title('z-score',fontsize=20)
    ax[0,3].set_title('z-score summed',fontsize=20)
    if max_summed_z_lim==-np.inf: max_summed_z_lim=10
    try:
        for i in range(0,16):
            # ax[15-i,1].set_ylim(-2,10)
            # ax[15-i,1].set_yticklabels([])
            ax[15-i,2].set_xticklabels([])
            # ax[15-i,2].set_ylim(-2,20)
            # ax[15-i,2].set_yticklabels([])
            ax[15-i,3].set_xticklabels([])
    except:
        breakpoint()
            
    ax[0,1].set_xticks([0,2])
    ax[0,1].set_xticklabels([0,2])  
    ax[0,2].set_xticks([0,2])
    ax[0,2].set_xticklabels([0,2])    
    # ax[0,1].set_yticks([0,10])
    # ax[0,1].set_yticklabels([0,10])
    ax[0,3].set_xticks([0,2])
    ax[0,3].set_xticklabels([0,2])    
    # ax[0,2].set_yticks([0,20])
    # ax[0,2].set_yticklabels([0,20])
        
    if show_plot: plt.show()
    return fig,units

def analyze_mua_by_channel_nostim(loc,common_details=None,chans=range(16)):
    # kwik_file
    kwik_file = [f for f in os.listdir(loc) if f.endswith('.kwik')]
    kwik_file = kwik_file[0]
    unit_details = get_unclustered_unit_details(os.path.join(loc,kwik_file))
    units = []
    
    # is it single stim or multi_stim
    n_stim=0
    for i,unit in enumerate(unit_details['units']):
        spike_times = unit['spike_time']
        bin_interval = 0.02 # 20ms
        binned_time= np.linspace(np.min(spike_times),np.max(spike_times),np.int((interval[1]-interval[0])/bin_interval+1))
        firing_rate = np.empty((timestamps[0::ratios].size,binned_time.size-1))
        unit_color = cmap(float(i)/len(unit_details['units']))
        for j,t in enumerate(timestamps[0::ratios]):
            # get the relevant spike timestamps
            spike_that_trial = spike_times-t
            spike_that_trial = spike_that_trial[np.bitwise_and(spike_that_trial>=interval[0],spike_that_trial<=interval[1])]
            temp = np.histogram(spike_that_trial,binned_time)
            firing_rate[j,:] = temp[0]
            ax[15-i,0].vlines(spike_that_trial,ymin=j,ymax=j+0.8,color=unit_color)
            
        m = np.mean(firing_rate,axis=0)/bin_interval
        sd = np.std(firing_rate,axis=0)/bin_interval
        m_prestim = np.mean(np.mean(firing_rate[:,binned_time[:-1]<0],axis=0))/bin_interval
        sd_prestim = np.mean(np.std(firing_rate[:,binned_time[:-1]<0],axis=0))/bin_interval
        z_score = np.divide(m-m_prestim,sd_prestim)
        norm_score = np.divide(m-m_prestim,m_prestim)
        
        if np.max(z_score)>max_z_lim: max_z_lim = np.max(z_score)
        offset = 0
        ax[15-i,2].plot([binned_time[0],binned_time[-1]],[offset,offset],'k--')
        ax[15-i,3].plot([binned_time[0],binned_time[-1]],[offset,offset],'k--')
        print('chan for subject',chans)
        if np.max(z_score)>min_z and i in chans:
            this_unit = {}
            ax[15-i,1].plot(binned_time[:-1],m,color=unit_color,linewidth=3,alpha=0.5)
            ax[15-i,2].plot(binned_time[:-1],offset+z_zoom*z_score,color=unit_color,linewidth=3,alpha=0.5)
            max_z_for_unit = []
            summed_z_for_unit = []
            time_base = []
            # get sequential max_values
            for k in range(n_stim):
                that_time_bin = np.bitwise_and(binned_time[:-1]>=0.1*k,binned_time[:-1]<=0.1*(k+1))
                max_z = np.max(z_score[that_time_bin])
                max_z_ind = np.argmax(z_score[that_time_bin])
                total_z = np.sum(z_score[that_time_bin])
                if np.max(total_z)>max_summed_z_lim: max_summed_z_lim = np.max(total_z)
                ax[15-i,2].plot(binned_time[:-1][that_time_bin][max_z_ind],offset+z_zoom*max_z,color=unit_color,marker='*')
                ax[15-i,3].plot(0.1*k+0.05,offset+total_z,color=unit_color,marker='x')
                max_z_for_unit.append(max_z)
                summed_z_for_unit.append(total_z)
                time_base.append(0.1*k+0.05)
            this_unit['mean_activity_prestim'] = m_prestim
            this_unit['mean_activity_unnormalized'] = m
            this_unit['std_activity_unnormalized'] = sd
            this_unit['std_activity_prestim'] = sd_prestim
            this_unit['max_z'] = np.asarray(max_z_for_unit)
            this_unit['summed_z'] = np.asarray(summed_z_for_unit)
            if n_stim==20:
                this_unit['stable_response_ms'] = np.mean(this_unit['max_z'][-4:])/this_unit['max_z'][0]
                this_unit['first_pulse_ratio_ms'] = this_unit['max_z'][1]/this_unit['max_z'][0]
                this_unit['second_pulse_ratio_ms'] = this_unit['max_z'][2]/this_unit['max_z'][0]
            else:
                this_unit['stable_response_ms'] = np.nan
                this_unit['first_pulse_ratio_ms'] = np.nan
                this_unit['second_pulse_ratio_ms'] = np.nan
            this_unit['max_response'] = np.max(max_z_for_unit)
            this_unit['channel_number'] = i
            this_unit['z_score'] = z_score
            this_unit['norm_score'] = norm_score
            this_unit['istimi'] = np.mean(np.diff(timestamps))
            this_unit['n_stim'] = n_stim
            
            time_base = np.asarray(time_base)
            
            # look at off responses
            if n_stim==20:
                time_bin_off = np.bitwise_and(binned_time[:-1]>=2.0,binned_time[:-1]<=2.5) # 2-2.5 s
            else:
                time_bin_off = np.bitwise_and(binned_time[:-1]>=0.2,binned_time[:-1]<=1) # 0.2 - 1 sec
            total_off_z = np.sum(z_score[time_bin_off])
            this_unit['off_response'] = total_off_z
            
            # fit to exponential_func
            try:
                if n_stim==20:
                    max_z = this_unit['max_z']
                    popt, pcov = curve_fit(exponential_func, time_base, max_z, p0=(max_z[0], 0.5, 1),bounds=([max_z[0]*0.95,1./800,0],[max_z[0]*1.05,np.inf,np.inf]))
                    xx = np.linspace(0, 2, 1000)
                    yy = exponential_func(xx, *popt)
                    ax[15-i,2].plot(xx,yy+offset,color=unit_color,linewidth=2)
                    this_unit['max_z_t_ms'] = 1/popt[1]
                    #ax[i,1].text(xx[-1],offset,'{0:.2f}'.format(1000/popt[1]),horizontalalignment='center',verticalalignment='bottom',fontsize=20)
                    this_unit['max_z_t_quality_ms'] = get_r_squared(max_z,exponential_func(time_base, *popt))
                else: 
                    this_unit['max_z_t_ms'] = np.nan
                    this_unit['max_z_t_quality_ms'] = np.nan
            except RuntimeError:
                print('Runtime')
                this_unit['max_z_t_ms'] = np.nan
                this_unit['max_z_t_quality_ms'] = np.nan
            except ValueError: 
                this_unit['max_z_t_ms'] = np.nan
                this_unit['max_z_t_quality_ms'] = np.nan
            
            try:
                if n_stim==20:
                    summed_z = this_unit['summed_z']
                    popt, pcov = curve_fit(exponential_func, time_base, summed_z, p0=(summed_z[0], 0.5, 1),bounds=([summed_z[0]*0.95,1./800,0],[summed_z[0]*1.05,np.inf,np.inf]))
                    xx = np.linspace(0, 2, 1000)
                    yy = exponential_func(xx, *popt)
                    ax[15-i,3].plot(xx,yy+offset,color=unit_color,linewidth=2)
                    this_unit['total_z_t_ms'] = 1/popt[1]
                    ax[15-i,3].text(xx[-1],2,'{0:.1f}'.format(1000/popt[1]),horizontalalignment='center',verticalalignment='bottom',fontsize=20)
                    this_unit['total_z_t_quality_ms'] = get_r_squared(max_z,exponential_func(time_base, *popt))
                else:
                    this_unit['max_z_t_ms'] = np.nan
                    this_unit['max_z_t_quality_ms'] = np.nan
            except RuntimeError:
                print('Runtime')
                this_unit['total_z_t_ms'] = np.nan
                this_unit['total_z_t_quality_ms'] = np.nan
            except ValueError: 
                this_unit['total_z_t_ms'] = np.nan
                this_unit['total_z_t_quality_ms'] = np.nan
                
            # add the common details
            if common_details:
                this_unit['subject_id'] = common_details['subject_id']
                this_unit['date'] = common_details['date']
                this_unit['time'] = common_details['time']
                this_unit['genotype'] = common_details['genotype']
                this_unit['stim_location'] = common_details['stim_location']
                this_unit['unit_id'] = common_details['subject_id'] + '_' + common_details['genotype'] + '_' + str(i)
                this_unit['repeat'] = common_details['repeat']
                if 'drug_detail' in common_details: this_unit['drug_detail'] = common_details['drug_detail']
                
            units.append(this_unit)
        else: ax[15-i,2].plot(binned_time[:-1],offset+z_zoom*z_score,color=(0.5,0.5,0.5))
        
        # get the ylim for ax[0]
        ylim0 = ax[15-i,1].get_ylim()
        ylim1 = ax[15-i,2].get_ylim()
        for j in range(n_stim):
            stim1 = matplotlib.patches.Rectangle((j*0.1,0.0),0.05,ylim0[1],color='blue',alpha=0.15)
            stim2 = matplotlib.patches.Rectangle((j*0.1,0.0),0.05,ylim1[1],color='blue',alpha=0.15)
            ax[15-i,1].add_patch(stim1)
            ax[15-i,2].add_patch(stim2)

    fig.suptitle(loc)
    ax[0,0].set_title('raster',fontsize=20)
    ax[0,2].set_title('z-score',fontsize=20)
    ax[0,3].set_title('z-score summed',fontsize=20)
    if max_summed_z_lim==-np.inf: max_summed_z_lim=10
    try:
        for i in range(0,16):
            # ax[15-i,1].set_ylim(-2,10)
            # ax[15-i,1].set_yticklabels([])
            ax[15-i,2].set_xticklabels([])
            # ax[15-i,2].set_ylim(-2,20)
            # ax[15-i,2].set_yticklabels([])
            ax[15-i,3].set_xticklabels([])
    except:
        breakpoint()
            
    ax[0,1].set_xticks([0,2])
    ax[0,1].set_xticklabels([0,2])  
    ax[0,2].set_xticks([0,2])
    ax[0,2].set_xticklabels([0,2])    
    # ax[0,1].set_yticks([0,10])
    # ax[0,1].set_yticklabels([0,10])
    ax[0,3].set_xticks([0,2])
    ax[0,3].set_xticklabels([0,2])    
    # ax[0,2].set_yticks([0,20])
    # ax[0,2].set_yticklabels([0,20])
        
    if show_plot: plt.show()
    return fig,units

    
def analyze_mua_by_channel_singlestim(loc,type='raster',show_plot=True,min_z=5,z_zoom=5,interval=[-1,1]):
    # kwik_file
    kwik_file = [f for f in os.listdir(loc) if f.endswith('.kwik')]
    kwik_file = kwik_file[0]
    unit_details = get_unclustered_unit_details(os.path.join(loc,kwik_file))
    units = []
    # timestamps
    timestamps = np.fromfile(os.path.join(loc,'timestamps.np'))
    ratios = 1
    n_stim = 1
    fig, ax = plt.subplots(16,2,figsize=(13,20))
    plt.subplots_adjust(wspace=0.1, hspace=0.1)
    cmap = matplotlib.cm.get_cmap('Dark2')
    for i,unit in enumerate(unit_details['units']):
        spike_times = unit['spike_time']
        bin_interval = 0.02 # 20ms
        binned_time= np.linspace(interval[0],interval[1],np.int((interval[1]-interval[0])/bin_interval+1))
        firing_rate = np.empty((timestamps[0::ratios].size,binned_time.size-1))
        unit_color = cmap(float(i)/len(unit_details['units']))
        for j,t in enumerate(timestamps[0::ratios]):
            # get the relevant spike timestamps
            spike_that_trial = spike_times-t
            spike_that_trial = spike_that_trial[np.bitwise_and(spike_that_trial>=interval[0],spike_that_trial<=interval[1])]
            temp = np.histogram(spike_that_trial,binned_time)
            firing_rate[j,:] = temp[0]
            ax[15-i,0].vlines(spike_that_trial,ymin=j,ymax=j+0.8,color=unit_color)
            
        m = np.mean(firing_rate,axis=0)/bin_interval
        sd = np.std(firing_rate,axis=0)/bin_interval
        m_prestim = np.mean(np.mean(firing_rate[:,binned_time[:-1]<0],axis=0))/bin_interval
        sd_prestim = np.mean(np.std(firing_rate[:,binned_time[:-1]<0],axis=0))/bin_interval
        z_score = np.divide(m-m_prestim,sd_prestim)
        offset = 0
        ax[15-i,1].plot([binned_time[0],binned_time[-1]],[offset,offset],'k--')
        if np.max(z_score)>min_z:
            this_unit = {}
            ax[15-i,1].plot(binned_time[:-1],offset+z_zoom*z_score,color=unit_color,linewidth=3,alpha=0.5)
            max_z_for_unit = []
            summed_z_for_unit = []
            time_base = []
            # get sequential max_values
            for k in range(n_stim):
                that_time_bin = np.bitwise_and(binned_time[:-1]>=0.1*k,binned_time[:-1]<=0.1*(k+1))
                max_z = np.max(z_score[that_time_bin])
                max_z_ind = np.argmax(z_score[that_time_bin])
                total_z = np.sum(z_score[that_time_bin])
                ax[15-i,1].plot(binned_time[:-1][that_time_bin][max_z_ind],offset+z_zoom*max_z,color=unit_color,marker='*')
                max_z_for_unit.append(max_z)
                summed_z_for_unit.append(total_z)
                time_base.append(0.1*k+0.05)
            this_unit['mean_activity_ss'] = np.asarray(max_z_for_unit)
            this_unit['std_activity_ss'] = np.asarray(max_z_for_unit)
            this_unit['max_z_ss'] = np.asarray(max_z_for_unit)
            this_unit['summed_z_ss'] = np.asarray(summed_z_for_unit)
            this_unit['channel_number'] = i
            time_base = np.asarray(time_base)
            
            # look at off responses
            time_bin_off = np.bitwise_and(binned_time[:-1]>=0.1,binned_time[:-1]<=1.0)
            total_off_z = np.mean(z_score[time_bin_off])
            this_unit['off_response_ss'] = total_off_z
            print(total_off_z)

            units.append(this_unit)
        else: ax[15-i,1].plot(binned_time[:-1],offset+z_zoom*z_score,color=(0.85,0.85,0.85))
        
        # get the ylim for ax[0]
        ylim = ax[15-i,0].get_ylim()
        ax[15-i,0].plot([0.1,0.1],ylim,'k--')
        ax[15-i,0].plot([0.2,0.2],ylim,'k--')
        ax[15-i,0].set_xticks([0,0.1])
        ax[15-i,0].set_xticklabels([])
        ax[15-i,1].set_ylim(0,50)
        ax[15-i,1].plot([0.1,0.1],[0,50],'k--')
        ax[15-i,1].plot([0.2,0.2],[0,50],'k--')
        ax[15-i,1].set_yticks([0,50])
        ax[15-i,1].set_xticks([0,0.1])
        ax[15-i,1].set_yticklabels([])
        ax[15-i,1].set_xticklabels([])
        
        if 15-i==0:
            ax[15-i,1].set_yticklabels([0,50])
    
    
    for i in range(n_stim):
        for k in range(16):
            ylim1 = ax[k,0].get_ylim()
            ylim2 = ax[k,1].get_ylim()
            stim1 = matplotlib.patches.Rectangle((i*0.1,0.0),0.05,ylim1[1],color='blue',alpha=0.15)
            stim2 = matplotlib.patches.Rectangle((i*0.1,0.0),0.05,ylim2[1],color='blue',alpha=0.15)
            ax[k,0].add_patch(stim1)
            ax[k,1].add_patch(stim2)
    fig.suptitle(loc)
    ax[0,0].set_title('raster')
    ax[0,1].set_title('z-score')
    if show_plot: plt.show()
    return fig,units

def plot_unit(unit,timestamps,interval,ax,type):
    all_ts = []
    for i,ts in enumerate(timestamps):
        relevant_ts = unit[np.bitwise_and((unit-ts)>=interval[0],(unit-ts)<=interval[1])]-ts
        if type=='raster': ax.vlines(relevant_ts,i+0.1,i+0.8,)
        all_ts.extend(relevant_ts)
    if type=='psth':
        ax.hist(all_ts,bins=300)
        return None
    if type=='density':
        return None
    if type=='spikecount_100ms':
        t_s = np.array([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])
        hist_c = np.histogram(all_ts,t_s)
        #print('histogram::',hist_c)
        ax.bar(t_s[0:-1],hist_c[0],width = 0.08)
        try:
            popt, pcov = curve_fit(exponential_func, t_s[0:-1], hist_c[0], p0=(1, 1e-6, 1))
            xx = np.linspace(0, 1, 1000)
            yy = exponential_func(xx, *popt)
            ax.plot(xx,yy,'r',linewidth=2)
            return popt[1]
        except RuntimeError:
            return None
    return None
    
def lighten_color(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    color = colorsys.rgb_to_hls(*mc.to_rgb(color))
    return colorsys.hls_to_rgb(color[0], 1 - amount * (1 - color[1]), color[2])
    
def plot_all_units(loc,type='raster',plot_on=True):
    # kwik_file
    kwik_file = [f for f in os.listdir(loc) if f.endswith('.kwik')]
    kwik_file = kwik_file[0]
    unit_details = get_unit_details(os.path.join(loc,kwik_file))
    
    # timestamps
    timestamps = np.fromfile(os.path.join(loc,'timestamps.np'))
    
    interval = [-1,2]
    num_units = len(unit_details['units'])
    nx_subplot = np.ceil(np.sqrt(num_units))
    ny_subplot = np.ceil(np.sqrt(num_units))
    f = plt.figure()
    time_scales = []
    for i,unit in enumerate(unit_details['units']):
        spike_times = unit['spike_time'].T[0]
        ax = plt.subplot(nx_subplot,ny_subplot,i+1)
        # ax = plt.axes()
        out = plot_unit(spike_times,timestamps,interval,ax,type)
        time_scales.append(out)
        ax.set_xticks([-1,0,1,2])
    # print(time_scales)
    if plot_on: plt.show()
    else: plt.close(f)
    return time_scales