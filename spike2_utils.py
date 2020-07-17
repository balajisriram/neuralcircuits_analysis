from __future__ import print_function
import pprint
import pickle
import sys
import os
import numpy as np
import h5py
import matplotlib.pyplot as plt
import matplotlib
from scipy.optimize import curve_fit
import matplotlib.colors as mc
import colorsys
import scipy.stats as stats
from tqdm import tqdm

ppr = pprint.PrettyPrinter(indent=2).pprint

def do_power_calculation(vals1,vals2,n_subsamples=[2,3,4,5,6,7,8,9,10,15,20,25,50,100],n_repeats=1000,alpha=0.05):
    vals1 = np.asarray(vals1).flatten()
    vals2 = np.asarray(vals2).flatten()
    power = []
    for n in tqdm(n_subsamples):
        succ = []
        for repeat in range(n_repeats):
            # sample n data points from vals1 nad vals2
            subsampled_v1 = np.random.choice(vals1,size=n,replace=True)
            subsampled_v2 = np.random.choice(vals2,size=n,replace=True)
            (t,p) = stats.ttest_ind(subsampled_v1,subsampled_v2)
            if p<alpha:
                succ.append(1)
            else:
                succ.append(0)
        succ = np.asarray(succ)
        power.append(succ.mean())
    return n_subsamples,power
    
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
    