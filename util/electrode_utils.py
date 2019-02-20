import xml.etree.ElementTree as ET
import json
import os
from scipy import spatial
from scipy.spatial.qhull import QhullError
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
from matplotlib.cm import get_cmap
from matplotlib.colors import Normalize
import numpy as np
import pprint
ppr = pprint.PrettyPrinter(indent=2).pprint


def get_headstage_map(name,loc=r'C:\Users\bsriram\Desktop\Code\neuralcircuits_analysis\mappings'):
    with open(os.path.join(loc,name+'.headstage_map'),'r') as f:
        headstage_map = json.load(f)
    return headstage_map

def get_electrode_map(name,loc=r'C:\Users\bsriram\Desktop\Code\neuralcircuits_analysis\mappings'):
    with open(os.path.join(loc,name+'.electrode_map'),'r') as f:
        electrode_map = json.load(f)
    return electrode_map

def get_impedance(xml_file):
    """
        xml file is numbered based on the channel 
    """
    
    tree = ET.parse(xml_file)
    impedance = {}
    phase = {}
    root = tree.getroot()
    for child in root:
        impedance[child.attrib['channel_number']] = child.attrib['magnitude']
        phase[child.attrib['channel_number']] = child.attrib['phase']
        
    return impedance,phase 

def get_graph_from_geometry(geometry):
    # let's transform the geometry into lists of channel names and coordinates
    chans,coords = zip(*[(ch,xy) for ch,xy in geometry.iteritems()])
    
    # we'll perform the triangulation
    try:
        tri = spatial.Delaunay(coords)
    except QhullError:
        chans,coords = list(chans),list(coords)
        x,y = zip(*coords)
        coords.append((max(x)+1,max(y)+1))
        tri = spatial.Delaunay(coords)
    
    # then build the list of edges from the triangulation
    indices, indptr = tri.vertex_neighbor_vertices
    edges = []
    for k in range(indices.shape[0]-1):
        for j in indptr[indices[k]:indices[k+1]]:
            try:
                edges.append((chans[k],chans[j]))
            except IndexError:
                pass
    return edges
    
def plot_edges_for_channel_group(edges,geometry,ax=None):
    if not ax:
        ax = plt.subplot()
    for edge in edges:
        e0 = edge[0]
        e1 = edge[1]
        ax.plot([geometry[e0][0],geometry[e1][0]],[geometry[e0][1],geometry[e1][1]],'k')    
        
    for p in geometry:
        ax.plot(geometry[p][0],geometry[p][1],'k.')
        ax.text(geometry[p][0],geometry[p][1],p)
    
def plot_edges(electrode_map):
    fig = plt.figure(figsize=(8.5,11))
    ax = plt.subplot()
    for ch_grp_num in electrode_map:
        if ch_grp_num in ['s','mapping','name']:
            continue
        else:
            # print(ch_grp_num)
            pass
        ch_grp = electrode_map[ch_grp_num]
        if not ch_grp['graph']:
            ch_grp['graph'] = get_graph_from_geometry(ch_grp['geometry'])
        plot_edges_for_channel_group(ch_grp['graph'],ch_grp['geometry'],ax=ax)
    ax.axis('tight')
    ax.axis('scaled')
    plt.show()
            
def plot_impedances(impedance, phase, electrode_map, headstage_map):
    fig = plt.figure()
    ax = plt.subplot()
    hs_mappings = np.squeeze(np.asarray(headstage_map['mapping']))
    electrode_mappings = np.squeeze(np.asarray(electrode_map['mapping']))
    geom = electrode_map['0']['geometry']
    cmap = get_cmap('bwr')
    cmap.set_over = (0.764,0.129,0.282)
    cmap.set_under = (0.2,0.2,0.6)
    cmap_norm = Normalize(vmin=-1,vmax=1)
    patches = []
    vals = []
    for ch in impedance:
        # ch is the channel on the recording
        #its location on the headstage is found by using the headstage_map
        imp = float(impedance[ch])
        loc_on_hs = np.squeeze(np.argwhere(hs_mappings==int(ch)))
        ch_on_electrode = electrode_mappings[loc_on_hs]
        #print(ch,':',loc_on_hs,':',ch_on_electrode,':(',geom[str(ch_on_electrode-1)][0],',',geom[str(ch_on_electrode-1)][1],')')
        (x,y) = (geom[str(ch_on_electrode)][0],geom[str(ch_on_electrode)][1] )
        col = cmap(cmap_norm(np.log10(imp/1000000.)))
        circ = Circle((x,y),7.5)#color=np.log10(imp/1000000.),alpha=0.5
        patches.append(circ)
        vals.append(cmap_norm(np.log10(imp/1000000.)))
    vals = np.array(vals)
    p = PatchCollection(patches)
    p.set_cmap(cmap)
    p.set_norm(cmap_norm)
    p.set_array(vals)
    ax.add_collection(p)
    cbar = fig.colorbar(p,ticks=[-1,-0.5,0,0.5,1])
    cbar.ax.set_yticklabels(['.1M','.3M','1M','3M','10M'])
    
    ax.axis('tight')
    ax.axis('scaled')
    plt.show()
    

        
    
