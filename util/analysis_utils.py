from electrode_utils import get_headstage_map,get_electrode_map,get_graph_from_geometry
import os

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
            f.write("                ]\n")
            
            ## make the "geometry" section
            f.write("        'geometry':{\n")
            for ch in geom:
                f.write("                s[{0}]:({1},{2}),\n".format(ch,geom[ch][0],geom[ch][1]))
            f.write("                }\n")
            
            f.write("      }\n")
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
        f.write("""
experiment_name = '{0}'
prb_file = '{1}'

traces = dict(
    raw_data_files=[experiment_name + '.dat'],
    voltage_gain={3},
    sample_rate={2},
    n_channels=32,
    dtype='int16',
)

spikedetekt = dict(
    filter_low=500.,  # Low pass frequency (Hz)
    filter_high_factor=0.95 * .5,
    filter_butter_order=3,  # Order of Butterworth filter.

    filter_lfp_low=0,  # LFP filter low-pass frequency
    filter_lfp_high=300,  # LFP filter high-pass frequency

    chunk_size_seconds=1,
    chunk_overlap_seconds=.015,

    n_excerpts=50,
    excerpt_size_seconds=1,
    threshold_strong_std_factor=4.5,
    threshold_weak_std_factor=2.,
    detect_spikes='both',

    connected_component_join_size=1,

    extract_s_before={4},
    extract_s_after={5},

    n_features_per_channel=3,  # Number of features per channel.
    pca_n_waveforms_max=10000,
)

klustakwik2 = dict(
    num_starting_clusters=100,
    max_iterations=1000,
    max_possible_clusters: 100,
)
        """.format(experiment_name,probe_file,sample_rate,voltage_gain,num_samp_before,num_samp_after))

def load_and_pack_dat(loc, rig='R1', n=32, outof=64):
    pass
        