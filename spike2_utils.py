from __future__ import print_function
import pprint
import pickle
import sys
import os
import numpy as np
import h5py
from klusta.kwik.model import KwikModel
import matplotlib.pyplot as plt
import matplotlib
from scipy.optimize import curve_fit
import matplotlib.colors as mc
import colorsys
import scipy.stats as stats

ppr = pprint.PrettyPrinter(indent=2).pprint

def create_session_dir_if_needed(base,filter):
    output_loc = os.path.join(base,filter)
    if not os.path.exists(output_loc):        
        os.makedirs(output_loc)
    return output_loc

def load_and_pack_WB(base_loc,sess,output_loc=None,recorded_sequence=np.arange(1,17),dref = None,):
    """
        base_loc: folder of dates containing phys data
        sess: session details: (subj, date, session, rig)
        output_loc: folder where output lives. if None, standard value is created at base_loc
        n_recorded : number of channels recorded
        dref: digital referencing. 'ave','med',None
        
    """
    print("")
    assert dref in ['ave','med',None], 'only None, average and median referencing allowed'
    
    (subj,recording_loc,loc_name,modality,stim_name) = sess
    filename_filter = '{0}_{1}_{2}_{3}'.format(subj,loc_name,modality,stim_name)
    print(filename_filter)
    mat_file = [f for f in os.listdir(os.path.join(base_loc)) if f.endswith('.mat') and f.startswith(filename_filter)]
    mat_file_path = os.path.join(base_loc,mat_file[0])
    if not output_loc:
        output_loc = create_session_dir_if_needed(base_loc,filename_filter)
    
    # channel names are unique - get the underlying fieldname
    with h5py.File(mat_file_path,'r') as f:
        fieldname = f.keys()[0]
        fieldbase = fieldname.split('_')
        fieldname = fieldbase[0]
        for base in fieldbase[1:-1]:
            fieldname = fieldname+'_'
            fieldname = fieldname+base
        print('Field name base is {0}'.format(fieldname))
        
        # get the max values
        max_value = np.NINF
        n_samples = np.inf
        for i,recorded_chan in enumerate(recorded_sequence):
            print('%d.' %recorded_chan ,end='')
            sys.stdout.flush()
            channel_name = '{0}_Ch{1}'.format(fieldname,recorded_chan)
            temp = f[channel_name]['values'][0]
            temp = np.asarray(temp)
            if np.max(np.abs(temp))>max_value:
                max_value = np.max(np.abs(temp))
            if temp.size<n_samples:
                n_samples = temp.size
        print("")
        print("maximum absolute value found = %s" % max_value)
        print("n_samples = {0}".format(n_samples))
        
        # create the data array
        data_array = np.zeros([n_samples, recorded_sequence.size], np.int16)
        
        
        print("loading, normalizing and packing")
        for i,recorded_chan in enumerate(recorded_sequence):
            print('%d.' %recorded_chan ,end='')
            sys.stdout.flush()
            channel_name = '{0}_Ch{1}'.format(fieldname,recorded_chan)
            temp = f[channel_name]['values'][0]
            temp = np.asarray(temp)
            temp = np.int16(temp/max_value*32767)
            try: data_array[:,i] = temp[0:n_samples]
            except: 
                import pdb
                pdb.set_trace()
        print("")
        
        voltage_scale=max_value/32767
    
        if dref:
            if dref == 'ave':
                print('Digital referencing to average of all channels.')
                reference = np.mean(data_array,1)
                for i in range(data_array.shape[1]):
                    data_array[:,i] = data_array[:,i] - reference
                dat_filename = subj+'_'+date+'_CAR.dat'
            elif dref == 'med':
                print ('Digital referencing to average of all channels.')
                reference = np.median(data_array,1)
                for i in range(data_array.shape[1]):
                    data_array[:,i] = data_array[:,i] - reference
                dat_filename = subj+'_'+date+'_CMR.dat'
        else:
            dat_filename = 'data.dat'
        
        data_array.tofile(os.path.join(output_loc,dat_filename))
        voltage_scale.tofile(os.path.join(output_loc,'voltage_scale.np'))
    
    return output_loc
    
def load_and_pack_timestamps(base_loc,sess,output_loc=None,):
    (subj,recording_loc,loc_name,modality,stim_name) = sess
    filename_filter = '{0}_{1}_{2}_{3}'.format(subj,loc_name,modality,stim_name)
    print(filename_filter)
    mat_file = [f for f in os.listdir(os.path.join(base_loc)) if f.endswith('.mat') and f.startswith(filename_filter)]
    mat_file_path = os.path.join(base_loc,mat_file[0])
    if not output_loc:
        output_loc = create_session_dir_if_needed(base_loc,filename_filter)
    # channel names are unique - get the underlying fieldname
    with h5py.File(mat_file_path,'r') as f:
        fieldname = f.keys()[0]
        fieldbase = fieldname.split('_')
        fieldname = fieldbase[0]
        for base in fieldbase[1:-1]:
            fieldname = fieldname+'_'
            fieldname = fieldname+base
        print('Field name base is {0}'.format(fieldname))
        try:
            t = f[fieldname+'_Ch19']['times'][0]
        except:
            for k in f.keys(): print(k)
            import pdb
            pdb.set_trace()
            raise Error()
        t.tofile(os.path.join(output_loc,'timestamps.np'))

def get_unit_details(loc):
    
    def get_cluster_waveform_all (kwik_model,cluster_id): 
        clusters = kwik_model.spike_clusters
        try:
            if not cluster_id in clusters:
                raise ValueError       
        except ValueError:
                print ("Exception: cluster_id (%d) not found !! " % cluster_id)
                return
        
        idx=np.argwhere(clusters==cluster_id)
        return kwik_model.all_waveforms[idx]
    
    # get the .kwik file and make KwikModel
    kwik_model = KwikModel(loc)
    kwik_model._open_kwik_if_needed()
    
    # all the stuff i want to extract
    output = {}
    output["n_samples_waveforms"] = kwik_model.n_samples_waveforms
    output["duration"] = kwik_model.duration
    output["num_channels"] = kwik_model.n_channels
    output["strong_threshold"] = kwik_model.metadata["threshold_strong_std_factor"]
    output["weak_threshold"] = kwik_model.metadata["threshold_weak_std_factor"]
    output["sample_rate"] = kwik_model.metadata["sample_rate"]
    output["high_pass_lo_f"] = kwik_model.metadata["filter_low"]
    output["high_pass_hi_f"] = kwik_model.metadata["filter_high_factor"]*output["sample_rate"]
    output["probe_name"] = kwik_model.metadata["prb_file"]
    output["data_type"] = kwik_model.metadata["dtype"]
    output["spikes_direction"] = kwik_model.metadata["detect_spikes"]
    output["klustakwik2_version"] = kwik_model.clustering_metadata["klustakwik2_version"]
    import pdb
    pdb.set_trace()
    ch_groups = kwik_model.channel_groups
    num_channel_groups = len(ch_groups)
    
    units = []
    
    for j,ch_grp in enumerate(ch_groups):
        print('Shank ', j+1,' of ',len(ch_groups))
        
        kwik_model.channel_group = ch_grp # go to that channel group
        kwik_model.clustering = 'main'
        
        cluster_quality = kwik_model.cluster_groups
        spike_time = kwik_model.spike_samples.astype(np.float64)/kwik_model.sample_rate
        spike_id = kwik_model.spike_ids
        cluster_id = kwik_model.spike_clusters
        cluster_ids = kwik_model.cluster_ids
        
        
        for i,clust_num in enumerate(cluster_ids):
            print('{0} of {1}'.format(i,cluster_ids.size))
            unit = {}
            unit["shank_no"] = ch_grp
            unit["cluster_id"] = clust_num
            unit["channels_in_shank"] = kwik_model.channels
            unit["manual_quality"] = cluster_quality[clust_num]
            # find the spike_ids corresponding to that cluster id
            spike_id_idx = np.argwhere(cluster_id==clust_num)
            spike_id_that_cluster = spike_id[spike_id_idx]
            unit["spike_ids"] = spike_id_that_cluster
            unit["spike_time"] = spike_time[spike_id_idx]
            
            
            waves = get_cluster_waveform_all(kwik_model,clust_num)
            mu_all = np.mean(waves[:,:,:],axis=0);
            std_all = np.std(waves[:,:,:],axis=0);
            unit["mean_waveform"] = np.mean(waves,axis=0)
            unit["std_waveform"] = np.std(waves,axis=0)
            
            unit["num_spikes"] = waves.shape[0]
            
            max_ind = np.unravel_index(np.argmin(mu_all),mu_all.shape)
            unit["loc_idx"] = max_ind
            max_ch = max_ind[1]
            
            unit["x_loc"] = kwik_model.channel_positions[max_ch,0]
            unit["y_loc"] = kwik_model.channel_positions[max_ch,1]
            
            units.append(unit)
    output["units"] = units
    
    kwik_model.close()
    return output
    
def get_unclustered_unit_details(loc):
    
    def get_cluster_waveform_all (kwik_model,cluster_id): 
        clusters = kwik_model.spike_clusters
        try:
            if not cluster_id in clusters:
                raise ValueError       
        except ValueError:
                print ("Exception: cluster_id (%d) not found !! " % cluster_id)
                return
        
        idx=np.argwhere(clusters==cluster_id)
        return kwik_model.all_waveforms[idx]
    
    # get the .kwik file and make KwikModel
    kwik_model = KwikModel(loc)
    kwik_model._open_kwik_if_needed()
    
    # all the stuff i want to extract
    output = {}
    output["n_samples_waveforms"] = kwik_model.n_samples_waveforms
    output["duration"] = kwik_model.duration
    output["num_channels"] = kwik_model.n_channels
    output["strong_threshold"] = kwik_model.metadata["threshold_strong_std_factor"]
    output["weak_threshold"] = kwik_model.metadata["threshold_weak_std_factor"]
    output["sample_rate"] = kwik_model.metadata["sample_rate"]
    output["high_pass_lo_f"] = kwik_model.metadata["filter_low"]
    output["high_pass_hi_f"] = kwik_model.metadata["filter_high_factor"]*output["sample_rate"]
    output["probe_name"] = kwik_model.metadata["prb_file"]
    output["data_type"] = kwik_model.metadata["dtype"]
    output["spikes_direction"] = kwik_model.metadata["detect_spikes"]
    ch_groups = kwik_model.channel_groups
    num_channel_groups = len(ch_groups)
    
    units = []
    
    for j,ch_grp in enumerate(ch_groups):
        print('Shank ', j+1,' of ',len(ch_groups))
        
        kwik_model.channel_group = ch_grp # go to that channel group
        
        spike_times = kwik_model.spike_times
        spike_id = kwik_model.spike_ids
        spike_channels = np.argmax(kwik_model.all_masks,axis=1)
        for i,chan_num in enumerate(kwik_model.channels):
            print('{0} of {1}'.format(i,np.max(kwik_model.channels)))
            unit = {}
            unit["shank_no"] = ch_grp
            unit["channel_no"] = chan_num
            spike_id_idx = np.argwhere(spike_channels==chan_num)
            spike_id_that_chan = spike_id[spike_id_idx]
            unit["spike_ids"] = spike_id_that_chan
            unit["spike_time"] = spike_times[spike_id_idx]
            unit["num_spikes"] = unit["spike_time"].shape[0]
            unit["x_loc"] = kwik_model.channel_positions[chan_num,0]
            unit["y_loc"] = kwik_model.channel_positions[chan_num,1]
            
            units.append(unit)
    output["units"] = units
    
    kwik_model.close()
    return output
    
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

def get_r_squared(data,fits):
    sse = np.sum(np.power(data-fits,2))
    sst = np.sum(np.power(data-data.mean(),2))
    return 1-sse/sst

def analyze_mua_by_channel_multistim(loc,type='raster',show_plot=True,min_z=5,z_zoom=1,interval=[-1,3]):
    # kwik_file
    kwik_file = [f for f in os.listdir(loc) if f.endswith('.kwik')]
    kwik_file = kwik_file[0]
    unit_details = get_unclustered_unit_details(os.path.join(loc,kwik_file))
    units = []
    # timestamps
    timestamps = np.fromfile(os.path.join(loc,'timestamps.np'))
    # get every other timestamp
    ratios = 2
    n_stim = 20
    fig, ax = plt.subplots(16,3,figsize=(20,20))
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
        if np.max(z_score)>max_z_lim: max_z_lim = np.max(z_score)
        offset = 0
        ax[15-i,1].plot([binned_time[0],binned_time[-1]],[offset,offset],'k--')
        ax[15-i,2].plot([binned_time[0],binned_time[-1]],[offset,offset],'k--')
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
                if np.max(total_z)>max_summed_z_lim: max_summed_z_lim = np.max(total_z)
                ax[15-i,1].plot(binned_time[:-1][that_time_bin][max_z_ind],offset+z_zoom*max_z,color=unit_color,marker='*')
                ax[15-i,2].plot(0.1*k+0.05,offset+total_z,color=unit_color,marker='x')
                max_z_for_unit.append(max_z)
                summed_z_for_unit.append(total_z)
                time_base.append(0.1*k+0.05)
            this_unit['mean_activity_ms'] = m_prestim
            this_unit['std_activity_ms'] = sd_prestim
            this_unit['max_z_ms'] = np.asarray(max_z_for_unit)
            this_unit['summed_z_ms'] = np.asarray(summed_z_for_unit)
            this_unit['stable_response_ms'] = np.mean(this_unit['max_z_ms'][-4:])/this_unit['max_z_ms'][0]
            this_unit['max_response_ms'] = np.max(max_z_for_unit)
            this_unit['first_pulse_ratio_ms'] = this_unit['max_z_ms'][1]/this_unit['max_z_ms'][0]
            this_unit['second_pulse_ratio_ms'] = this_unit['max_z_ms'][2]/this_unit['max_z_ms'][0]
            this_unit['channel_number_ms'] = i
            time_base = np.asarray(time_base)
            
            # look at off responses
            time_bin_off = np.bitwise_and(binned_time[:-1]>=2.0,binned_time[:-1]<=2.5)
            total_off_z = np.sum(z_score[time_bin_off])
            this_unit['off_response_ms'] = total_off_z
            
            # fit to exponential_func
            try:
                max_z = this_unit['max_z_ms']
                popt, pcov = curve_fit(exponential_func, time_base, max_z, p0=(max_z[0], 1e-6, 1),bounds=([max_z[0]*0.95,0,0],[max_z[0]*1.05,np.inf,np.inf]))
                xx = np.linspace(0, 2, 1000)
                yy = exponential_func(xx, *popt)
                ax[15-i,1].plot(xx,yy+offset,color=unit_color,linewidth=2)
                this_unit['max_z_t_ms'] = 1/popt[1]
                #ax[i,1].text(xx[-1],offset,'{0:.2f}'.format(1000/popt[1]),horizontalalignment='center',verticalalignment='bottom',fontsize=20)
                this_unit['max_z_t_quality_ms'] = get_r_squared(max_z,exponential_func(time_base, *popt))
            except RuntimeError:
                print('Runtime')
                this_unit['max_z_t_ms'] = np.nan
                this_unit['max_z_t_quality_ms'] = np.nan
            except ValueError: 
                this_unit['max_z_t_ms'] = np.nan
                this_unit['max_z_t_quality_ms'] = np.nan
            
            try:
                summed_z = this_unit['summed_z_ms']
                popt, pcov = curve_fit(exponential_func, time_base, summed_z, p0=(summed_z[0], 1e-6, 1),bounds=([summed_z[0]*0.95,0,0],[summed_z[0]*1.05,np.inf,np.inf]))
                xx = np.linspace(0, 2, 1000)
                yy = exponential_func(xx, *popt)
                ax[15-i,2].plot(xx,yy+offset,color=unit_color,linewidth=2)
                this_unit['total_z_t_ms'] = 1/popt[1]
                ax[15-i,2].text(xx[-1],2,'{0:.1f}'.format(1000/popt[1]),horizontalalignment='center',verticalalignment='bottom',fontsize=20)
                this_unit['total_z_t_quality_ms'] = get_r_squared(max_z,exponential_func(time_base, *popt))
            except RuntimeError:
                print('Runtime')
                this_unit['total_z_t_ms'] = np.nan
                this_unit['total_z_t_quality_ms'] = np.nan
            except ValueError: 
                this_unit['total_z_t_ms'] = np.nan
                this_unit['total_z_t_quality_ms'] = np.nan
            units.append(this_unit)
        else: ax[15-i,1].plot(binned_time[:-1],offset+z_zoom*z_score,color=(0.85,0.85,0.85))
        
        # get the ylim for ax[0]
        ylim = ax[i,0].get_ylim()
    
        for j in range(n_stim):
            stim1 = matplotlib.patches.Rectangle((j*0.1,0.0),0.05,ylim[1],color='blue',alpha=0.15)
            stim2 = matplotlib.patches.Rectangle((j*0.1,0.0),0.05,ylim[1],color='blue',alpha=0.15)
            ax[15-i,0].add_patch(stim1)
            ax[15-i,1].add_patch(stim2)
    fig.suptitle(loc)
    ax[0,0].set_title('raster',fontsize=20)
    ax[0,1].set_title('z-score',fontsize=20)
    ax[0,2].set_title('z-score summed',fontsize=20)
    if max_summed_z_lim==-np.inf: max_summed_z_lim=10
    try:
        for i in range(0,16):
            # ax[15-i,1].set_ylim(-2,np.ceil(1.1*max_z_lim))
            ax[15-i,1].set_ylim(-2,30)
            ax[15-i,1].set_yticklabels([])
            ax[15-i,1].set_xticklabels([])
            # ax[15-i,2].set_ylim(-2,np.ceil(1.1*max_summed_z_lim))
            ax[15-i,2].set_ylim(-2,40)
            ax[15-i,2].set_yticklabels([])
            ax[15-i,2].set_xticklabels([])
    except:
        import pdb;
        pdb.set_trace()
            
    ax[0,1].set_xticks([0,2])
    ax[0,1].set_xticklabels([0,2])    
    # ax[0,1].set_yticks([0,np.ceil(1.1*max_z_lim)])
    # ax[0,1].set_yticklabels([0,np.ceil(1.1*max_z_lim)])
    ax[0,1].set_yticks([0,30])
    ax[0,1].set_yticklabels([0,30])
    ax[0,2].set_xticks([0,2])
    ax[0,2].set_xticklabels([0,2])    
    # ax[0,2].set_yticks([0,np.ceil(1.1*max_summed_z_lim)])
    # ax[0,2].set_yticklabels([0,np.ceil(1.1*max_summed_z_lim)])
    ax[0,2].set_yticks([0,40])
    ax[0,2].set_yticklabels([0,40])
        
    if show_plot: plt.show()
    return fig,units
    
def analyze_mua_by_channel_multistim_old(loc,type='raster',show_plot=True,min_z=5,z_zoom=2,interval=[-1,3]):
    # kwik_file
    kwik_file = [f for f in os.listdir(loc) if f.endswith('.kwik')]
    kwik_file = kwik_file[0]
    unit_details = get_unclustered_unit_details(os.path.join(loc,kwik_file))
    units = []
    # timestamps
    timestamps = np.fromfile(os.path.join(loc,'timestamps.np'))
    # get every other timestamp
    ratios = 2
    n_stim = 20
    fig, ax = plt.subplots(16,3,figsize=(20,20),sharey=True)
    plt.subplots_adjust(left=0.125, bottom=0.1, right=0.9, top=0.9, wspace=0.2, hspace=0.2)

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
            ax[0].vlines(spike_that_trial,ymin=i*timestamps.size/ratios+j,ymax=i*timestamps.size/ratios+j+0.8,color=unit_color)
            
        m = np.mean(firing_rate,axis=0)/bin_interval
        sd = np.std(firing_rate,axis=0)/bin_interval
        m_prestim = np.mean(np.mean(firing_rate[:,binned_time[:-1]<0],axis=0))/bin_interval
        sd_prestim = np.mean(np.std(firing_rate[:,binned_time[:-1]<0],axis=0))/bin_interval
        z_score = np.divide(m-m_prestim,sd_prestim)
        offset = i*timestamps.size/ratios
        ax[1].plot([binned_time[0],binned_time[-1]],[offset,offset],'k--')
        ax[2].plot([binned_time[0],binned_time[-1]],[offset,offset],'k--')
        if np.max(z_score)>min_z:
            this_unit = {}
            ax[1].plot(binned_time[:-1],offset+z_zoom*z_score,color=unit_color,linewidth=3,alpha=0.5)
            max_z_for_unit = []
            summed_z_for_unit = []
            time_base = []
            # get sequential max_values
            for k in range(n_stim):
                that_time_bin = np.bitwise_and(binned_time[:-1]>=0.1*k,binned_time[:-1]<=0.1*(k+1))
                max_z = np.max(z_score[that_time_bin])
                max_z_ind = np.argmax(z_score[that_time_bin])
                total_z = np.sum(z_score[that_time_bin])
                ax[1].plot(binned_time[:-1][that_time_bin][max_z_ind],offset+z_zoom*max_z,color=unit_color,marker='*')
                ax[2].plot(0.1*k+0.05,offset+total_z,color=unit_color,marker='x')
                max_z_for_unit.append(max_z)
                summed_z_for_unit.append(total_z)
                time_base.append(0.1*k+0.05)
            this_unit['mean_activity_ms'] = m_prestim
            this_unit['std_activity_ms'] = sd_prestim
            this_unit['max_z_ms'] = np.asarray(max_z_for_unit)
            this_unit['summed_z_ms'] = np.asarray(summed_z_for_unit)
            this_unit['stable_response_ms'] = np.mean(this_unit['max_z_ms'][-4:])/this_unit['max_z_ms'][0]
            this_unit['max_response_ms'] = np.max(max_z_for_unit)
            this_unit['first_pulse_ratio_ms'] = this_unit['max_z_ms'][1]/this_unit['max_z_ms'][0]
            this_unit['second_pulse_ratio_ms'] = this_unit['max_z_ms'][2]/this_unit['max_z_ms'][0]
            this_unit['channel_number_ms'] = i
            time_base = np.asarray(time_base)
            
            # look at off responses
            time_bin_off = np.bitwise_and(binned_time[:-1]>=2.0,binned_time[:-1]<=2.5)
            total_off_z = np.sum(z_score[time_bin_off])
            this_unit['off_response_ms'] = total_off_z
            
            # fit to exponential_func
            try:
                max_z = this_unit['max_z_ms']
                popt, pcov = curve_fit(exponential_func, time_base, max_z, p0=(max_z[0], 1e-6, 1),bounds=([max_z[0]*0.95,0,0],[max_z[0]*1.05,np.inf,np.inf]))
                xx = np.linspace(0, 2, 1000)
                yy = exponential_func(xx, *popt)
                ax[1].plot(xx,yy+offset,color=unit_color,linewidth=2)
                this_unit['max_z_t_ms'] = 1/popt[1]
                ax[1].text(xx[-1],offset,'{0:.2f}'.format(1000/popt[1]),horizontalalignment='center',verticalalignment='bottom',fontsize=20)
                this_unit['max_z_t_quality_ms'] = get_r_squared(max_z,exponential_func(time_base, *popt))
            except RuntimeError:
                print('Runtime')
                this_unit['max_z_t_ms'] = np.nan
                this_unit['max_z_t_quality_ms'] = np.nan
            except ValueError: 
                this_unit['max_z_t_ms'] = np.nan
                this_unit['max_z_t_quality_ms'] = np.nan
            
            try:
                summed_z = this_unit['summed_z_ms']
                popt, pcov = curve_fit(exponential_func, time_base, summed_z, p0=(summed_z[0], 1e-6, 1),bounds=([summed_z[0]*0.95,0,0],[summed_z[0]*1.05,np.inf,np.inf]))
                xx = np.linspace(0, 2, 1000)
                yy = exponential_func(xx, *popt)
                ax[2].plot(xx,yy+offset,color=unit_color,linewidth=2)
                this_unit['total_z_t_ms'] = 1/popt[1]
                ax[2].text(xx[-1],offset,'{0:.2f}'.format(1000/popt[1]),horizontalalignment='center',verticalalignment='bottom',fontsize=20)
                this_unit['total_z_t_quality_ms'] = get_r_squared(max_z,exponential_func(time_base, *popt))
            except RuntimeError:
                print('Runtime')
                this_unit['total_z_t_ms'] = np.nan
                this_unit['total_z_t_quality_ms'] = np.nan
            except ValueError: 
                this_unit['total_z_t_ms'] = np.nan
                this_unit['total_z_t_quality_ms'] = np.nan
            units.append(this_unit)
        else: ax[1].plot(binned_time[:-1],offset+z_zoom*z_score,color=(0.85,0.85,0.85))
        
    # get the ylim for ax[0]
    ylim = ax[0].get_ylim()
    
    for i in range(n_stim):
        stim1 = matplotlib.patches.Rectangle((i*0.1,0.0),0.05,ylim[1],color='blue',alpha=0.15)
        stim2 = matplotlib.patches.Rectangle((i*0.1,0.0),0.05,ylim[1],color='blue',alpha=0.15)
        ax[0].add_patch(stim1)
        ax[1].add_patch(stim2)
    fig.suptitle(loc)
    ax[0].set_title('raster',fontsize=20)
    ax[1].set_title('z-score',fontsize=20)
    ax[2].set_title('z-score summed',fontsize=20)
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
    
def analyze_mua_by_channel_singlestim_old(loc,type='raster',show_plot=True,min_z=5,z_zoom=5,interval=[-1,1]):
    # kwik_file
    kwik_file = [f for f in os.listdir(loc) if f.endswith('.kwik')]
    kwik_file = kwik_file[0]
    unit_details = get_unclustered_unit_details(os.path.join(loc,kwik_file))
    units = []
    # timestamps
    timestamps = np.fromfile(os.path.join(loc,'timestamps.np'))
    ratios = 1
    n_stim = 1
    fig, ax = plt.subplots(1,2,figsize=(20,20),sharey=True)
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
            ax[0].vlines(spike_that_trial,ymin=i*timestamps.size/ratios+j,ymax=i*timestamps.size/ratios+j+0.8,color=unit_color)
            
        m = np.mean(firing_rate,axis=0)/bin_interval
        sd = np.std(firing_rate,axis=0)/bin_interval
        m_prestim = np.mean(np.mean(firing_rate[:,binned_time[:-1]<0],axis=0))/bin_interval
        sd_prestim = np.mean(np.std(firing_rate[:,binned_time[:-1]<0],axis=0))/bin_interval
        z_score = np.divide(m-m_prestim,sd_prestim)
        offset = i*timestamps.size/ratios
        ax[1].plot([binned_time[0],binned_time[-1]],[offset,offset],'k--')
        if np.max(z_score)>min_z:
            this_unit = {}
            ax[1].plot(binned_time[:-1],offset+z_zoom*z_score,color=unit_color,linewidth=3,alpha=0.5)
            max_z_for_unit = []
            summed_z_for_unit = []
            time_base = []
            # get sequential max_values
            for k in range(n_stim):
                that_time_bin = np.bitwise_and(binned_time[:-1]>=0.1*k,binned_time[:-1]<=0.1*(k+1))
                max_z = np.max(z_score[that_time_bin])
                max_z_ind = np.argmax(z_score[that_time_bin])
                total_z = np.sum(z_score[that_time_bin])
                ax[1].plot(binned_time[:-1][that_time_bin][max_z_ind],offset+z_zoom*max_z,color=unit_color,marker='*')
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
        else: ax[1].plot(binned_time[:-1],offset+z_zoom*z_score,color=(0.85,0.85,0.85))
        
    # get the ylim for ax[0]
    ylim = ax[0].get_ylim()
    ax[0].plot([0.1,0.1],ylim,'k--')
    ax[0].plot([0.2,0.2],ylim,'k--')
    ax[1].plot([0.1,0.1],ylim,'k--')
    ax[1].plot([0.2,0.2],ylim,'k--')
    
    
    for i in range(n_stim):
        stim1 = matplotlib.patches.Rectangle((i*0.1,0.0),0.05,ylim[1],color='blue',alpha=0.15)
        stim2 = matplotlib.patches.Rectangle((i*0.1,0.0),0.05,ylim[1],color='blue',alpha=0.15)
        ax[0].add_patch(stim1)
        ax[1].add_patch(stim2)
    fig.suptitle(loc)
    ax[0].set_title('raster')
    ax[1].set_title('z-score')
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

def exponential_func(x, a, b, c):
    return a*np.exp(-b*x)+c

def write_qsub_file(base,fold):
    with open(os.path.join(base,fold,'submit.qsub'),'w') as f:
        f.write('#!/bin/bash\n')
        f.write('#$ -N {0}\n'.format(fold))
        f.write('#$ -l h_rt=24:00:00\n')
        f.write('#$ -q long.q\n')
        f.write('#$ -wd /home/bsriram/data/{0}\n'.format(fold))
        f.write('#$ -j no\n')
        f.write('#$ -M balaji.sriram@biogen.com\n')
        f.write('#$ -m be\n')
        f.write('#$ -e error.log\n')
        f.write('#$ -o output.log\n')
        f.write('#$ -pe openmpi-fillup 12\n')
        f.write('\n')
        f.write('##########################################################################\n')
        f.write('# Start your script here #\n')
        f.write('##########################################################################\n')
        f.write('# Load the modules you need.\n')
        f.write('source /home/bsriram/miniconda3/bin/activate klusta\n')
        f.write('# Run some commands.\n')
        f.write('export LC_ALL=en_US.UTF-8\n')
        f.write('export LANG=en_US.UTF-8\n')
        f.write('export LANGUAGE=en_US.UTF-8\n')
        f.write('klusta --overwrite --detect-only params.prm\n')
        f.write('\n')
        f.write('rm -rf .spikedetekt\n')
        f.write('rm -rf .klustakwik2\n')
        f.write('# Exit successfully.\n')
        f.write('exit 0\n')
    return 0

def write_params_file(base,fold):
    with open(os.path.join(base,fold,'params.prm'),'w') as f:
        f.write("experiment_name = 'data'\n")
        f.write("prb_file = 'NNX_16Ch_Linear_SE.prb'\n")
        f.write("traces = dict(\n")
        f.write("    raw_data_files=[experiment_name + '.dat'],\n")
        f.write("    voltage_gain=0.523007705575,\n")
        f.write("    sample_rate=25000.0,\n")
        f.write("    n_channels=16,\n")
        f.write("    dtype='int16',\n")
        f.write(")\n")
        f.write("spikedetekt = dict(\n")
        f.write("    filter_low=500.,  # Low pass frequency (Hz)\n")
        f.write("    filter_high_factor=0.95 * .5,\n")
        f.write("    filter_butter_order=3,  # Order of Butterworth filter.\n")
        f.write("    filter_lfp_low=0,  # LFP filter low-pass frequency\n")
        f.write("    filter_lfp_high=300,  # LFP filter high-pass frequency\n")
        f.write("    chunk_size_seconds=1,\n")
        f.write("    chunk_overlap_seconds=.015,\n")
        f.write("    n_excerpts=50,\n")
        f.write("    excerpt_size_seconds=1,\n")
        f.write("    threshold_strong_std_factor=3.5,\n")
        f.write("    threshold_weak_std_factor=2.,\n")
        f.write("    detect_spikes='negative',\n")
        f.write("    connected_component_join_size=1,\n")
        f.write("    extract_s_before=40,\n")
        f.write("    extract_s_after=40,\n")
        f.write("    n_features_per_channel=3,  # Number of features per channel.\n")
        f.write("    pca_n_waveforms_max=10000,\n")
        f.write(")\n")
        f.write("klustakwik2 = dict(\n")
        f.write("    num_starting_clusters=100,\n")
        f.write("    max_iterations=1000,\n")
        f.write("    max_possible_clusters=100,\n")
        f.write(")\n")
    return 0

def plot_feature_across_areas(allData,feature,split_by_subj=True,sig_thresh=0.0025,plot_on=True,label='',filter_hi=None,filter_lo=None):
    if plot_on: fig,ax = plt.subplots(1,4,figsize=(20,5),sharey=True)
    locations = ['V1','LGN','S1','VPLM']
    Ts = []
    Ps = []
    Effects = []
    Vals_hom = []
    Vals_wt = []
    for i,location in enumerate(locations):
        hom = allData.loc[(allData['location']==location) & (allData['genotype']=='HOM'),feature]
        hom_sub = allData.loc[(allData['location']==location) & (allData['genotype']=='HOM'),'subject']
        if filter_hi:
            which = hom<filter_hi
            hom = hom[which]
            hom_sub = hom_sub[which]
        if filter_lo:
            which = hom>filter_lo
            hom = hom[which]
            hom_sub = hom_sub[which]
        wt = allData.loc[(allData['location']==location) & (allData['genotype']=='WT'),feature]
        wt_sub = allData.loc[(allData['location']==location) & (allData['genotype']=='WT'),'subject']
        if filter_hi:
            which = wt<filter_hi
            wt = wt[which]
            wt_sub = wt_sub[which]
        if filter_lo:
            which = wt>filter_lo
            wt = wt[which]
            wt_sub = wt_sub[which]
        hom_by_sub = []
        wt_by_sub = []
    
        # plot wild_type
        curr_color = (0.,0.,1.)
        for j,sub in enumerate(np.unique(wt_sub)):
            if plot_on: ax[i].scatter(0.9-0.25*(j+1)+0.1*np.random.rand(1,wt[wt_sub==sub].size),wt[wt_sub==sub],color=curr_color)
            wt_by_sub.append(np.nanmean(wt[wt_sub==sub]))
            curr_color=lighten_color(curr_color,amount=0.8)
        # plot homs
        curr_color = (1.,0.,0.)
        for j,sub in enumerate(np.unique(hom_sub)):
            if plot_on:ax[i].scatter(2.1+0.25*(j+1)+0.1*np.random.rand(1,hom[hom_sub==sub].size),hom[hom_sub==sub],color=curr_color)
            hom_by_sub.append(np.nanmean(hom[hom_sub==sub]))
            curr_color=lighten_color(curr_color,amount=0.8)
        if plot_on: ax[i].boxplot([wt,hom],showmeans=True)
        if plot_on: ax[i].set_title(location+' - '+feature)
        if plot_on: ax[i].set_xlim([-2,5])
        if plot_on: ax[i].set_xticks([0.75,2.25])
        if plot_on: ax[i].set_xticklabels(['WT','KO'])
        
        (t,p) = stats.ttest_ind(hom_by_sub,wt_by_sub)
        effect = np.mean(hom_by_sub)/np.mean(wt_by_sub)
        
        

        Ts.append(t)
        Ps.append(p)
        Effects.append(effect)
        Vals_hom.append(hom_by_sub)
        Vals_wt.append(wt_by_sub)
        
    ax[0].set_ylabel(label,fontsize=15)
    fig.suptitle(feature)
    return Ts,Ps,Effects,Vals_hom,Vals_wt


if __name__=='__main__':
    data = np.arange(1000)+1000000*(np.random.rand(1,1000)-0.5)
    fits = np.arange(1000)
    print(get_r_squared(data,fits))
    plt.scatter(data,fits)
    plt.show()
    # loc = r'C:\Users\bsriram\Desktop\Data\PGRN\Sessions\PGRN_HOM_1_Loc1_VisStim_B'
    # plot_mua_by_channel(loc)
    # Multi stim response
    # # All V1
    # timescale_i = []
    # ts_that = plot_all_units(r'C:\Users\bsriram\Desktop\Data\PGRN\Sessions\PGRN_WT_3_Loc2_VisStim_D',type='spikecount_100ms',plot_on = False)
    # timescale_i.extend(ts_that)
    # ts_that = plot_all_units(r'C:\Users\bsriram\Desktop\Data\PGRN\Sessions\PGRN_WT_4_Loc2_VisStim_D',type='spikecount_100ms',plot_on = False)
    # timescale_i.extend(ts_that)
    # ts_that = plot_all_units(r'C:\Users\bsriram\Desktop\Data\PGRN\Sessions\PGRN_WT_5_Loc2_VisStim_D',type='spikecount_100ms',plot_on = False)
    # timescale_i.extend(ts_that)
    # ts_that = plot_all_units(r'C:\Users\bsriram\Desktop\Data\PGRN\Sessions\PGRN_WT_6_Loc2_VisStim_D',type='spikecount_100ms',plot_on = False)
    # timescale_i.extend(ts_that)
    # ts_that = plot_all_units(r'C:\Users\bsriram\Desktop\Data\PGRN\Sessions\PGRN_WT_7_Loc2_VisStim_D',type='spikecount_100ms',plot_on = False)
    # timescale_i.extend(ts_that)
    # ts_that = plot_all_units(r'C:\Users\bsriram\Desktop\Data\PGRN\Sessions\PGRN_WT_8_Loc2_VisStim_D',type='spikecount_100ms',plot_on = False)
    # timescale_i.extend(ts_that)
    # ts_that = plot_all_units(r'C:\Users\bsriram\Desktop\Data\PGRN\Sessions\PGRN_WT_9_Loc4_VisStim_D',type='spikecount_100ms',plot_on = False)
    # timescale_i.extend(ts_that)
    # ts_that = plot_all_units(r'C:\Users\bsriram\Desktop\Data\PGRN\Sessions\PGRN_WT_10_Loc3_VisStim_D',type='spikecount_100ms',plot_on = False)
    # timescale_i.extend(ts_that)
    # ts_that = plot_all_units(r'C:\Users\bsriram\Desktop\Data\PGRN\Sessions\PGRN_WT_11_Loc3_VisStim_D',type='spikecount_100ms',plot_on = False)
    # timescale_i.extend(ts_that)
    
    # timescales = []
    # for t in timescale_i:
        # if not t: timescales.append(np.nan)
        # elif t<0: timescales.append(-1/t)
        # else: timescales.append(1/t)
    # timescales = np.array(timescales)

    # timescales_good = timescales[~np.isnan(timescales)]
    # timescales_good = timescales_good[timescales_good<3]
    # plt.hist(timescales_good,bins=20)
    # plt.show()
    
    # # All S1
    # timescale_i = []
    # ts_that = plot_all_units(r'C:\Users\bsriram\Desktop\Data\PGRN\Sessions\PGRN_WT_3_Loc1_whiskerstim_C',type='spikecount_100ms',plot_on = False)
    # timescale_i.extend(ts_that)
    # ts_that = plot_all_units(r'C:\Users\bsriram\Desktop\Data\PGRN\Sessions\PGRN_WT_4_Loc1_whiskerstim_C',type='spikecount_100ms',plot_on = False)
    # timescale_i.extend(ts_that)
    # ts_that = plot_all_units(r'C:\Users\bsriram\Desktop\Data\PGRN\Sessions\PGRN_WT_5_Loc1_whiskerstim_C',type='spikecount_100ms',plot_on = False)
    # timescale_i.extend(ts_that)
    # ts_that = plot_all_units(r'C:\Users\bsriram\Desktop\Data\PGRN\Sessions\PGRN_WT_6_Loc1_whiskerstim_C',type='spikecount_100ms',plot_on = False)
    # timescale_i.extend(ts_that)
    # ts_that = plot_all_units(r'C:\Users\bsriram\Desktop\Data\PGRN\Sessions\PGRN_WT_7_Loc1_whiskerstim_C',type='spikecount_100ms',plot_on = False)
    # timescale_i.extend(ts_that)
    # ts_that = plot_all_units(r'C:\Users\bsriram\Desktop\Data\PGRN\Sessions\PGRN_WT_8_Loc1_whiskerstim_C',type='spikecount_100ms',plot_on = False)
    # timescale_i.extend(ts_that)
    # ts_that = plot_all_units(r'C:\Users\bsriram\Desktop\Data\PGRN\Sessions\PGRN_WT_9_Loc1_whiskerstim_C',type='spikecount_100ms',plot_on = False)
    # timescale_i.extend(ts_that)
    # ts_that = plot_all_units(r'C:\Users\bsriram\Desktop\Data\PGRN\Sessions\PGRN_WT_10_Loc1_whiskerstim_C',type='spikecount_100ms',plot_on = False)
    # timescale_i.extend(ts_that)
    # ts_that = plot_all_units(r'C:\Users\bsriram\Desktop\Data\PGRN\Sessions\PGRN_WT_11_Loc1_whiskerstim_C',type='spikecount_100ms',plot_on = False)
    # timescale_i.extend(ts_that)
    
    # timescales = []
    # for t in timescale_i:
        # if not t: timescales.append(np.nan)
        # elif t<0: timescales.append(-1/t)
        # else: timescales.append(1/t)
    # timescales = np.array(timescales)

    # timescales_good = timescales[~np.isnan(timescales)]
    # timescales_good = timescales_good[timescales_good<3]
    # plt.hist(timescales_good,bins=20)
    # plt.show()
    
    # All dLG
    # timescale_i = []
    # ts_that = plot_all_units(r'C:\Users\bsriram\Desktop\Data\PGRN\Sessions\PGRN_WT_3_Loc3_VisStim_D',type='spikecount_100ms',plot_on = False)
    # timescale_i.extend(ts_that)
    # ts_that = plot_all_units(r'C:\Users\bsriram\Desktop\Data\PGRN\Sessions\PGRN_WT_4_Loc3_VisStim_D',type='spikecount_100ms',plot_on = False)
    # timescale_i.extend(ts_that)
    # ts_that = plot_all_units(r'C:\Users\bsriram\Desktop\Data\PGRN\Sessions\PGRN_WT_5_Loc3_VisStim_D',type='spikecount_100ms',plot_on = False)
    # timescale_i.extend(ts_that)
    # ts_that = plot_all_units(r'C:\Users\bsriram\Desktop\Data\PGRN\Sessions\PGRN_WT_6_Loc4_VisStim_D',type='spikecount_100ms',plot_on = False)
    # timescale_i.extend(ts_that)
    # ts_that = plot_all_units(r'C:\Users\bsriram\Desktop\Data\PGRN\Sessions\PGRN_WT_7_Loc3_VisStim_D',type='spikecount_100ms',plot_on = False)
    # timescale_i.extend(ts_that)
    # ts_that = plot_all_units(r'C:\Users\bsriram\Desktop\Data\PGRN\Sessions\PGRN_WT_8_Loc4_VisStim_D',type='spikecount_100ms',plot_on = False)
    # timescale_i.extend(ts_that)
    # ts_that = plot_all_units(r'C:\Users\bsriram\Desktop\Data\PGRN\Sessions\PGRN_WT_10_Loc4_VisStim_D',type='spikecount_100ms',plot_on = False)
    # timescale_i.extend(ts_that)
    # ts_that = plot_all_units(r'C:\Users\bsriram\Desktop\Data\PGRN\Sessions\PGRN_WT_10_Loc4_VisStim_D2',type='spikecount_100ms',plot_on = False)
    # timescale_i.extend(ts_that)
    # ts_that = plot_all_units(r'C:\Users\bsriram\Desktop\Data\PGRN\Sessions\PGRN_WT_11_Loc4_VisStim_D',type='spikecount_100ms',plot_on = False)
    # timescale_i.extend(ts_that)
    
    # timescales = []
    # for t in timescale_i:
        # if not t: timescales.append(np.nan)
        # elif t<0: timescales.append(-1/t)
        # else: timescales.append(1/t)
    # timescales = np.array(timescales)

    # timescales_good = timescales[~np.isnan(timescales)]
    # timescales_good = timescales_good[timescales_good<3]
    # plt.hist(timescales_good,bins=20)
    # plt.show()
    
    # All VPLM
    # timescale_i = []
    # ts_that = plot_all_units(r'C:\Users\bsriram\Desktop\Data\PGRN\Sessions\PGRN_WT_9_Loc2_whiskerstim_C',type='spikecount_100ms',plot_on = False)
    # timescale_i.extend(ts_that)
    # ts_that = plot_all_units(r'C:\Users\bsriram\Desktop\Data\PGRN\Sessions\PGRN_WT_9_Loc3_whiskerstim_C',type='spikecount_100ms',plot_on = False)
    # timescale_i.extend(ts_that)
    # ts_that = plot_all_units(r'C:\Users\bsriram\Desktop\Data\PGRN\Sessions\PGRN_WT_10_Loc2_whiskerstim_C',type='spikecount_100ms',plot_on = False)
    # timescale_i.extend(ts_that)
    # ts_that = plot_all_units(r'C:\Users\bsriram\Desktop\Data\PGRN\Sessions\PGRN_WT_11_Loc2_whiskerstim_C',type='spikecount_100ms',plot_on = False)
    # timescale_i.extend(ts_that)
    
    # timescales = []
    # for t in timescale_i:
        # if not t: timescales.append(np.nan)
        # elif t<0: timescales.append(-1/t)
        # else: timescales.append(1/t)
    # timescales = np.array(timescales)

    # timescales_good = timescales[~np.isnan(timescales)]
    # timescales_good = timescales_good[timescales_good<3]
    # plt.hist(timescales_good,bins=20)
    # plt.show()
    
    # All S1
    # plot_all_units(r'C:\Users\bsriram\Desktop\Data\PGRN\Sessions\PGRN_WT_3_Loc2_VisStim_D',type='spikecount_100ms')
    