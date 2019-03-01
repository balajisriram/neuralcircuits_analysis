from __future__ import print_function
from pypl2 import pl2_ad, pl2_spikes, pl2_events, pl2_info
import pprint
import pickle
import sys
import os
import numpy as np

ppr = pprint.PrettyPrinter(indent=2).pprint

def create_session_dir_if_needed(base,subj,date,output_name='output'):
    output_loc = os.path.join(base,output_name,subj+'_'+date)
    if not os.path.exists(output_loc):        
        os.makedirs(output_loc)
    return output_loc
    

def verify_names(filename,type='WB',contains='R1',n_req=32):
    spkinfo, evtinfo, adinfo = pl2_info(filename)
    
    signal_names = []
    n_samples = []
    
    if type=='WB':
        for n in range(len(adinfo)):
            if contains+'_'+type in adinfo[n].name and adinfo[n].n>0:
                signal_names.append(adinfo[n].name)
                if not n_samples:
                    n_samples = adinfo[n].n
                elif n_samples != adinfo[n].n:
                    ValueError('a channel does not have the same number of data points as another')
                else:
                    pass #all ok
        
        n_channels = len(signal_names)
        assert n_channels==n_req, 'obtained incommensurate number of channels. n_req=={0}, n_actual={1}'.format(n_req,n_channels)
        return signal_names,(n_samples,n_channels)
    elif type=='Running':
        for n in range(len(adinfo)):
            if contains in adinfo[n].name and adinfo[n].n>0:
                signal_names.append(adinfo[n].name)
                if not n_samples:n_samples = adinfo[n].n
                elif n_samples != adinfo[n].n:ValueError('a channel does not have the same number of data points as another')
                else:pass
        n_channels = len(signal_names)
        assert n_channels==n_req, 'obtained incommensurate number of channels. n_req=={0}, n_actual={1}'.format(n_req,n_channels)
        return signal_names,(n_samples,n_channels)
    elif type=='Events':
        for n in range(len(evtinfo)):
            if contains in evtinfo[n].name:
                signal_names.append(evtinfo[n].name)
                n_samples.append(evtinfo[n].n)
        return signal_names,(n_samples,len(n_samples))
        
        
def load_and_pack_WB(base_loc,sess,output_loc=None,n_recorded=32,dref = None,):
    """
        base_loc: folder of dates containing phys data
        sess: session details: (subj, date, session, rig)
        output_loc: folder where output lives. if None, standard value is created at base_loc
        n_recorded : number of channels recorded
        dref: digital referencing. 'ave','med',None
        
    """
    
    assert dref in ['ave','med'], 'only average and median referncing allowed'
    
    (subj,date,session,rig,running_channel,heastage,electrode) = sess
    pl2_file = [f for f in os.listdir(os.path.join(base_loc,date,session)) if f.endswith('.pl2')]
    pl2_file_path = os.path.join(base_loc,date,session,pl2_file[0])
    signal_names,details = verify_names(pl2_file_path, type='WB',contains=rig, n_req=n_recorded)
    
    if not output_loc:
        output_loc = create_session_dir_if_needed(base_loc,subj,date)
        
    n_samples,n_channels = details
    data_array = np.zeros([n_samples, n_channels], np.int16)
    max_value = np.NINF
    for i,signal_name in enumerate(signal_names):
        print('%d.' %i ,end='')
        sys.stdout.flush()
        temp_f,n,temp_ts,temp_frag_counts,temp = pl2_ad(pl2_file_path, signal_name)
        temp = np.asarray(temp)
        if np.max(np.abs(temp))>max_value:
            max_value = np.max(np.abs(temp))
    print("")
    print("maximum absolute value found = %s" % max_value)
    
    print("loading, normalizing and packing")
    for i,signal_name in enumerate(signal_names):
        print('%d.' %i ,end='')
        sys.stdout.flush()
        temp_f,n,temp_ts,temp_frag_counts,temp = pl2_ad(pl2_file_path, signal_name)
        temp = np.asarray(temp)
        temp = np.int16(temp/max_value*32767)
        data_array[:,i] = temp
    print("")
    wb_timestamp = temp_ts
    wb_timestamp = np.asarray(wb_timestamp)
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
        dat_filename = subj+'_'+date+'.dat'
    
    data_array.tofile(os.path.join(output_loc,dat_filename))
    voltage_scale.tofile(os.path.join(output_loc,'voltage_scale.np'))
    wb_timestamp.tofile(os.path.join(output_loc,'wb_timestamp.np'))
    
    return output_loc

def load_and_pack_running(base_loc,sess,output_loc=None):
    (subj,date,session,rig,running_channel,heastage,electrode) = sess
    pl2_file = [f for f in os.listdir(os.path.join(base_loc,date,session)) if f.endswith('.pl2')]
    pl2_file_path = os.path.join(base_loc,date,session,pl2_file[0])
    if running_channel=='RunningNotNamed':
        if rig=='R1': running_chan_details='AIF01'
        if rig=='R2': running_chan_details='AIF02'
    signal_names,details = verify_names(pl2_file_path, type='Running',contains=running_chan_details, n_req=1)
    
    if not output_loc:
        output_loc = create_session_dir_if_needed(base_loc,subj,date)
        
    running_f, n, running_ts, frag_count, running = pl2_ad(pl2_file_path, signal_names[0]) # only one
    running_f = np.asarray(running_f)
    running = np.asarray(running)
    running_ts = np.asarray(running_ts)
    running_f.tofile(os.path.join(output_loc,'running_samplerate.np'))
    running.tofile(os.path.join(output_loc,'running_speed.np'))
    running_ts.tofile(os.path.join(output_loc,'running_ts.np'))
    
def load_and_pack_events(base_loc,sess,output_loc=None):
    (subj,date,session,rig,running_channel,heastage,electrode) = sess
    pl2_file = [f for f in os.listdir(os.path.join(base_loc,date,session)) if f.endswith('.pl2')]
    pl2_file_path = os.path.join(base_loc,date,session,pl2_file[0])
    
    signal_names,details = verify_names(pl2_file_path, type='Events',contains=rig)
    
    if not output_loc:
        output_loc = create_session_dir_if_needed(base_loc,subj,date)
        
    event_data = {}
    for signal_name in signal_names:
        nums,ts,vals = pl2_events(pl2_file_path,signal_name)
        event_data[signal_name] = np.asarray(ts)
    
    with open(os.path.join(output_loc,'event_data.pickle'),'wb') as f:
        pickle.dump(event_data,f,protocol=pickle.HIGHEST_PROTOCOL)        