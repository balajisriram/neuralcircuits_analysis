from electrode_utils import get_headstage_map,get_electrode_map,get_graph_from_geometry
import os
import numpy as np

np64 = np.float64

def make_prb_file(headstage,electrode,save_path,mapping_loc=r'C:\Users\bsriram\Desktop\Code\neuralcircuits_analysis\mappings',probe_name=None):
    if not probe_name:
        probe_name = electrode+'_'+headstage+'.prb'
    
    hs_map = get_headstage_map(headstage,mapping_loc)
    elec_map = get_electrode_map(electrode,mapping_loc)
    
    with open(os.path.join(save_path,probe_name),'w') as f:
        f.write("# s[electrode_channel]->saved_raw_data_channel\n")
        f.write("s={\n")
        hs_map_order = hs_map['mapping']
        elec_map_order = elec_map['mapping']
        
        for hs_chan,elec_chan in zip(hs_map_order,elec_map_order):
            f.write("    {0}:{1},\n".format(elec_chan,hs_chan))
        f.write("  }\n")
        f.write("\n")
        f.write("\n")
        
        f.write("# now deal with channel groups\n")
        f.write("channel_groups = {\n")
        ch_grp_names = list(elec_map.keys())
        for key in ['mapping','s','name']:
            try:ch_grp_names.remove(key)
            except KeyError:pass
        for ch_grp_name in ch_grp_names:
            f.write("    {0}:".format(ch_grp_name))
            f.write("{\n")
            
            # prep the data for writing
            ch_grp = elec_map[ch_grp_name]
            geom = ch_grp['geometry']
            graph = ch_grp['graph']
            if not graph:graph = get_graph_from_geometry(geom)

            ## make the "channels" section
            f.write("        'channels':[")
            for ch in list(geom.keys()):f.write("s[{0}],".format(ch))
            f.write("],\n")
            
            ## make the "graph" section
            f.write("        'graph':[\n")
            for edge in graph:
                f.write("                (s[{0}],s[{1}]),".format(edge[0],edge[1]))
                f.write("\n")
            f.write("                ],\n")
            
            ## make the "geometry" section
            f.write("        'geometry':{\n")
            for ch in geom:
                f.write("                s[{0}]:({1},{2}),\n".format(ch,geom[ch][0],geom[ch][1]))
            f.write("                }\n")
            
            f.write("      },\n")
        f.write("}\n")
    return probe_name
    
def make_prm_file(session_location, sample_rate=40000., experiment_name=None, probe_file=None, voltage_gain=1.):

    dat_files = [f for f in os.listdir(session_location) if f.endswith('.dat')]
    if len(dat_files) !=1:
        RuntimeError('Too many or too few dat files. length={0}'.format(len(dat_files)))
    else:
        dat_file_name = dat_files[0]
    
    if not experiment_name:
        experiment_name, ext = os.path.splitext(dat_file_name)
             
    if not probe_file:
        prb_files = [f for f in os.listdir(session_location) if f.endswith('.prb')]
        if len(prb_files) !=1:
            RuntimeError('Too many or too few files. length={0}'.format(len(dat_files)))
        else:
            probe_file = prb_files[0]
    import math
    num_samp_before = math.ceil(sample_rate/1000)
    num_samp_after = num_samp_before
    with open(os.path.join(session_location, 'params.prm'),'w') as f:
        f.write("experiment_name = '{0}'\n".format(experiment_name))
        f.write("prb_file = '{0}'\n".format(probe_file))
        f.write("traces = dict(\n")
        f.write("    raw_data_files=[experiment_name + '.dat'],\n")
        f.write("    voltage_gain={0},\n".format(voltage_gain[0]))
        f.write("    sample_rate={0},\n".format(sample_rate))
        f.write("    n_channels=32,\n")
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
        f.write("    threshold_strong_std_factor=4.5,\n")
        f.write("    threshold_weak_std_factor=2.,\n")
        f.write("    detect_spikes='negative',\n")
        f.write("    connected_component_join_size=1,\n")
        f.write("    extract_s_before={0},\n".format(int(num_samp_before)))
        f.write("    extract_s_after={0},\n".format(int(num_samp_after)))
        f.write("    n_features_per_channel=3,  # Number of features per channel.\n")
        f.write("    pca_n_waveforms_max=10000,\n")
        f.write(")\n")
        f.write("klustakwik2 = dict(\n")
        f.write("    num_starting_clusters=100,\n")
        f.write("    max_iterations=1000,\n")
        f.write("    max_possible_clusters=100,\n")
        f.write(")\n")

def make_qsub_file(session_location, base_on_server='/home/bsriram/data_biogen/',):

    local_base,session_folder_name = os.path.split(session_location)
    server_wd = base_on_server+session_folder_name

    with open(os.path.join(session_location, 'submit.qsub'),mode='w') as f:
        f.write("#!/bin/bash\n")
        f.write("#$ -N {0}\n".format(session_folder_name+'_kwik'))
        f.write("#$ -l h_rt=24:00:00\n")
        f.write("#$ -q long.q\n")
        f.write("#$ -wd {0}\n".format(server_wd))
        f.write("#$ -j no\n")
        f.write("#$ -M balaji.sriram@biogen.com\n")
        f.write("#$ -m be\n")
        f.write("#$ -e error.log\n")
        f.write("#$ -o output.log\n")
        f.write("#$ -pe openmpi-fillup 12\n")
        f.write("##########################################################################\n")
        f.write("# Start your script here #\n")
        f.write("##########################################################################\n")
        f.write("# Load the modules you need.\n")
        f.write("export LC_ALL=en_US.utf8\n")
        f.write("export LANG=en_US.utf8\n")
        f.write("source /home/bsriram/miniconda3/bin/activate klusta\n")
        f.write("# Run some commands.\n")
        f.write("klusta --overwrite --detect-only params.prm\n")
        f.write("rm -rf .spikedetekt\n")
        f.write("klusta --cluster-only params.prm\n")
        f.write("# Exit successfully.\n")
        f.write("exit 0\n")

def make_qsub_file_1(session_location, base_on_server='/home/bsriram/data_biogen/',):

    local_base,session_folder_name = os.path.split(session_location)
    server_wd = base_on_server+session_folder_name
    with open(os.path.join(session_location, 'submit.qsub'),mode='w',newline='\n') as f:
        f.write(file_string)
    file_string = """#!/bin/bash
#$ -N {0}
#$ -l h_rt=24:00:00
#$ -q long.q
#$ -wd {1}
#$ -j no
#$ -M balaji.sriram@biogen.com
#$ -m be
#$ -e error.log
#$ -o output.log
#$ -pe openmpi-fillup 12
##########################################################################
# Start your script here #
##########################################################################
# Load the modules you need.
export LC_ALL=en_US.utf8
export LANG=en_US.utf8
source /home/bsriram/miniconda3/bin/activate klusta
# Run some commands.
klusta params.prm
# Exit successfully.
exit 0""".format(session_folder_name,server_wd)
    
def load_spike_and_trial_details(loc):
    spike_and_trial_details = {}
    print('Loading trial details')
    spike_and_trial_details['trial_records'] = load_trialrecs_to_dict(loc)
    print('Loading spike details')
    spike_and_trial_details['spike_records'] = load_spikerecs_to_dict(loc)
    return spike_and_trial_details
       
def load_trialrecs_to_dict(loc):
    file = [f for f in os.listdir(loc) if f.startswith('trialRecords')]
    if len(file) > 1 or len(file)==0:
        print(loc)
        error('too many or too few trial Records. how is this possible?')
    
    temp = scipy.io.loadmat(os.path.join(loc,file[0]))
    pdb.set_trace()
    tRs = temp['trialRecords'][0]
    numTrials = tRs.size
    
    trial_number = []
    refresh_rate = []
    step_name = []
    stim_manager_name = []
    trial_manager_name = []
    afc_grating_type = []
    trial_start_index = []
    trial_end_index = []
    
    pix_per_cycs = []
    driftfrequency = []
    orientation = []
    phase = []
    contrast = []
    max_duration = []
    radius = []
    annulus = []
    waveform = []
    radius_type = []
    location = []
    
    led_on = []
    led_intensity = []
    
    events = find_indices_for_all_trials(loc)
    events = get_channel_events(loc,events)
    for i,tR in enumerate(tRs):
        this_trial_number = np64(tR['trialNumber'][0][0])
        if this_trial_number in events['trial_number']:
            which_in_messages = [True if x==this_trial_number else False for x in events['trial_number']]
            this_start_index = [x for i,x in enumerate(events['start_index']) if which_in_messages[i]==True]
            this_end_index = [x for i,x in enumerate(events['end_index']) if which_in_messages[i]==True]
        else:
            continue
        
        trial_number.append(np64(tR['trialNumber'][0][0]))
        refresh_rate.append(np64(tR['refreshRate'][0][0]))
        
        step_name.append(tR['stepName'][0])
        stim_manager_name.append(tR['stimManagerClass'][0])
        trial_manager_name.append(tR['trialManagerClass'][0])
        afc_grating_type.append(tR['afcGratingType'][0])
        
        trial_start_index.append(this_start_index[0])
        trial_end_index.append(this_end_index[0])
        
        try: pix_per_cycs.append(np64(tR['stimulus'][0]['pixPerCyc'][0][0][0]))
        except: pix_per_cycs.append(np64(tR['stimulus'][0]['pixPerCycs'][0][0][0]))
        
        try: driftfrequency.append(np64(tR['stimulus'][0]['driftfrequency'][0][0][0]))
        except: driftfrequency.append(np64(tR['stimulus'][0]['driftfrequencies'][0][0][0]))
        
        try: orientation.append(np64(tR['stimulus'][0]['orientation'][0][0][0]))
        except: orientation.append(np64(tR['stimulus'][0]['orientations'][0][0][0]))
        
        try: phase.append(np64(tR['stimulus'][0]['phase'][0][0][0]))
        except: phase.append(np64(tR['stimulus'][0]['phases'][0][0][0]))
        
        try: contrast.append(np64(tR['stimulus'][0]['contrast'][0][0][0]))
        except: contrast.append(np64(tR['stimulus'][0]['contrasts'][0][0][0]))
        
        max_duration.append(np64(tR['stimulus'][0]['maxDuration'][0][0][0]))
        
        try: radius.append(np64(tR['stimulus'][0]['radius'][0][0][0]))
        except: radius.append(np64(tR['stimulus'][0]['radii'][0][0][0]))
        
        try: annulus.append(np64(tR['stimulus'][0]['annulus'][0][0][0]))
        except: annulus.append(np64(tR['stimulus'][0]['annuli'][0][0][0]))
        
        waveform.append(tR['stimulus'][0]['waveform'][0][0])
        radius_type.append(tR['stimulus'][0]['radiusType'][0][0])
        
        location.append(np64(tR['stimulus'][0]['location'][0][0]))
        
        led_on.append(np64(tR['LEDON'][0][0]))
        led_intensity.append(np64(tR['LEDIntensity'][0][0]))
        
    trial_records = dict([('trial_number',trial_number),\
                         ('refresh_rate',refresh_rate),\
                         ('step_name',step_name),\
                         ('stim_manager_name',stim_manager_name),\
                         ('trial_manager_name',trial_manager_name),\
                         ('afc_grating_type',afc_grating_type),\
                         ('trial_start_index',trial_start_index),\
                         ('trial_end_index',trial_end_index),\
                         ('pix_per_cycs',pix_per_cycs),\
                         ('driftfrequency',driftfrequency),\
                         ('orientation',orientation),\
                         ('phase',phase),\
                         ('contrast',contrast),\
                         ('max_duration',max_duration),\
                         ('radius',radius),\
                         ('annulus',annulus),\
                         ('waveform',waveform),\
                         ('radius_type',radius_type),\
                         ('location',location),\
                         ('led_on',led_on),\
                         ('led_intensity',led_intensity),\
                         ('events',events)])
                         
    return trial_records

def read_trial_records():
    pass