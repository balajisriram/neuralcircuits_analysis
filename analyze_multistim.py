from spike2_utils import analyze_mua_by_channel_multistim
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import pickle
import pandas as pd

allData = pd.DataFrame()

# V1
locs_V1_HOM = [r'C:\Data\PGRN_Old\PGRN_HOM_1_Loc1_VisStim_D',
        r'C:\Data\PGRN_Old\PGRN_HOM_2_Loc1_VisStim_D',
        r'C:\Data\PGRN_Old\PGRN_HOM_3_Loc1_VisStim_D',
        r'C:\Data\PGRN_Old\PGRN_HOM_4_Loc1_VisStim_D',
        r'C:\Data\PGRN_Old\PGRN_HOM_5_Loc1_VisStim_D',]
V1_HOM_subjects = ['HOM1','HOM2','HOM3','HOM4','HOM5']
V1_HOM_location = ['Loc1','Loc1','Loc1','Loc1','Loc1']
locs_V1_WT = [r'C:\Data\PGRN_Old\PGRN_WT_10_Loc3_VisStim_D',
        r'C:\Data\PGRN_Old\PGRN_WT_11_Loc3_VisStim_D',
        r'C:\Data\PGRN_Old\PGRN_WT_15_Loc1_VisStim_D',
        r'C:\Data\PGRN_Old\PGRN_WT_16_Loc1_VisStim_D',
        r'C:\Data\PGRN_Old\PGRN_WT_17_Loc1_VisStim_D',
        r'C:\Data\PGRN_Old\PGRN_WT_17_Loc2_VisStim_D',
        r'C:\Data\PGRN_Old\PGRN_WT_18_Loc1_VisStim_D',
        r'C:\Data\PGRN_Old\PGRN_WT_19_Loc1_VisStim_D',]
V1_WT_subjects = ['WT10','WT11','WT15','WT16','WT17','WT17','WT18','WT19']
V1_WT_location = ['Loc3','Loc3','Loc1','Loc1','Loc1','Loc2','Loc1','Loc1']
with PdfPages('V1_multistim.pdf') as pdf:
    f,ax = plt.subplots(1,1,figsize=(20,20))
    ax.text(0.5,0.5,'V1 Homozygous animals',horizontalalignment='center',verticalalignment='center',fontsize=48)
    pdf.savefig(f)
    plt.close()
    for loc,subj,penet in zip(locs_V1_HOM,V1_HOM_subjects,V1_HOM_location):
        f,units = analyze_mua_by_channel_multistim(loc,show_plot=False,min_z=2.5)
        curr_units = pd.DataFrame(units)
        curr_units['subject'] = subj
        curr_units['genotype'] = 'HOM'
        curr_units['location'] = 'V1'
        curr_units['penetration'] = penet
        allData = allData.append(curr_units,ignore_index=True,sort=False)
        pdf.savefig(f)
        plt.close()
    f,ax = plt.subplots(1,1,figsize=(20,20))
    ax.text(0.5,0.5,'V1 WT animals',horizontalalignment='center',verticalalignment='center',fontsize=48)
    pdf.savefig(f)
    plt.close()
    for loc,subj,penet in zip(locs_V1_WT,V1_WT_subjects,V1_WT_location):
        f,units = analyze_mua_by_channel_multistim(loc,show_plot=False,min_z=2.5)
        curr_units = pd.DataFrame(units)
        curr_units['subject'] = subj
        curr_units['genotype'] = 'WT'
        curr_units['location'] = 'V1'
        curr_units['penetration'] = penet
        allData = allData.append(curr_units,ignore_index=True,sort=False)
        pdf.savefig(f)
        plt.close()

# LGN
locs_LGN_HOM = [r'C:\Data\PGRN_Old\PGRN_HOM_1_Loc2_VisStim_D',
        r'C:\Data\PGRN_Old\PGRN_HOM_2_Loc2_VisStim_D',
        r'C:\Data\PGRN_Old\PGRN_HOM_3_Loc2_VisStim_D',
        r'C:\Data\PGRN_Old\PGRN_HOM_4_Loc2_VisStim_D',
        r'C:\Data\PGRN_Old\PGRN_HOM_5_Loc2_VisStim_D',
        r'C:\Data\PGRN_Old\PGRN_HOM_5_Loc3_VisStim_D',]
lgn_HOM_subjects = ['HOM1','HOM2','HOM3','HOM4','HOM5','HOM5']
lgn_HOM_location = ['Loc2','Loc2','Loc2','Loc2','Loc2','Loc3']
locs_LGN_WT = [r'C:\Data\PGRN_Old\PGRN_WT_10_Loc4_VisStim_D',
        r'C:\Data\PGRN_Old\PGRN_WT_10_Loc4_VisStim_D2',
        r'C:\Data\PGRN_Old\PGRN_WT_11_Loc4_VisStim_D',
        r'C:\Data\PGRN_Old\PGRN_WT_11_Loc5_VisStim_D',
        r'C:\Data\PGRN_Old\PGRN_WT_15_Loc3_VisStim_D',
        r'C:\Data\PGRN_Old\PGRN_WT_16_Loc2_VisStim_D',
        r'C:\Data\PGRN_Old\PGRN_WT_17_Loc3_VisStim_D',
        r'C:\Data\PGRN_Old\PGRN_WT_18_Loc2_VisStim_D',
        r'C:\Data\PGRN_Old\PGRN_WT_18_Loc3_VisStim_D',
        r'C:\Data\PGRN_Old\PGRN_WT_19_Loc2_VisStim_D',]
lgn_WT_subjects = ['WT10','WT10','WT11','WT11','WT15','WT16','WT17','WT18','WT18','WT19']
lgn_WT_location = ['Loc4','Loc4','Loc4','Loc5','Loc3','Loc2','Loc3','Loc2','Loc3','Loc2']
with PdfPages('LGN_multistim.pdf') as pdf:
    f,ax = plt.subplots(1,1,figsize=(20,20))
    ax.text(0.5,0.5,'LGN Homozygous animals',horizontalalignment='center',verticalalignment='center',fontsize=48)
    pdf.savefig(f)
    plt.close()
    for loc,subj,penet in zip(locs_LGN_HOM,lgn_HOM_subjects,lgn_HOM_location):
        f,units = analyze_mua_by_channel_multistim(loc,show_plot=False,min_z=2.5)
        curr_units = pd.DataFrame(units)
        curr_units['subject'] = subj
        curr_units['genotype'] = 'HOM'
        curr_units['location'] = 'LGN'
        curr_units['penetration'] = penet
        allData = allData.append(curr_units,ignore_index=True,sort=False)
        pdf.savefig(f)
        plt.close()
        
    f,ax = plt.subplots(1,1,figsize=(20,20))
    ax.text(0.5,0.5,'LGN WT animals',horizontalalignment='center',verticalalignment='center',fontsize=48)
    pdf.savefig(f)
    plt.close()
    wt_lgn_units = []
    for loc,subj,penet in zip(locs_LGN_WT,lgn_WT_subjects,lgn_WT_location):
        f,units = analyze_mua_by_channel_multistim(loc,show_plot=False,min_z=2.5)
        curr_units = pd.DataFrame(units)
        curr_units['subject'] = subj
        curr_units['genotype'] = 'WT'
        curr_units['location'] = 'LGN'
        curr_units['penetration'] = penet
        allData = allData.append(curr_units,ignore_index=True,sort=False)
        pdf.savefig(f)
        plt.close()

# S1
locs_S1_HOM = [r'C:\Data\PGRN_Old\PGRN_HOM_1_Loc3_WhiskStim_C',
        r'C:\Data\PGRN_Old\PGRN_HOM_2_Loc3_WhiskStim_C',
        r'C:\Data\PGRN_Old\PGRN_HOM_3_Loc3_WhiskStim_C',
        r'C:\Data\PGRN_Old\PGRN_HOM_4_Loc5_WhiskStim_C',
        r'C:\Data\PGRN_Old\PGRN_HOM_4_Loc5_WhiskStim_C1',
        r'C:\Data\PGRN_Old\PGRN_HOM_5_Loc5_WhiskStim_C',]
S1_HOM_subjects = ['HOM1','HOM2','HOM3','HOM4','HOM4','HOM5']
S1_HOM_location = ['Loc3','Loc3','Loc3','Loc5','Loc5','Loc5']
locs_S1_WT = [r'C:\Data\PGRN_Old\PGRN_WT_10_Loc1_whiskerstim_C',
        r'C:\Data\PGRN_Old\PGRN_WT_11_Loc1_whiskerstim_C',
        r'C:\Data\PGRN_Old\PGRN_WT_15_Loc4_WhiskStim_C',
        r'C:\Data\PGRN_Old\PGRN_WT_16_Loc3_WhiskStim_C',
        r'C:\Data\PGRN_Old\PGRN_WT_17_Loc4_WhiskStim_C',
        r'C:\Data\PGRN_Old\PGRN_WT_18_Loc4_WhiskStim_C',
        r'C:\Data\PGRN_Old\PGRN_WT_19_Loc3_WhiskStim_C',]
S1_WT_subjects = ['WT10','WT11','WT15','WT16','WT17','WT18','WT19']
S1_WT_location = ['Loc1','Loc1','Loc4','Loc3','Loc4','Loc4','Loc3']
with PdfPages('S1_multistim.pdf') as pdf:
    f,ax = plt.subplots(1,1,figsize=(20,20))
    ax.text(0.5,0.5,'S1 Homozygous animals',horizontalalignment='center',verticalalignment='center',fontsize=48)
    pdf.savefig(f)
    plt.close()
    for loc,subj,penet in zip(locs_S1_HOM,S1_HOM_subjects,S1_HOM_location):
        f,units = analyze_mua_by_channel_multistim(loc,show_plot=False,min_z=2.5)
        curr_units = pd.DataFrame(units)
        curr_units['subject'] = subj
        curr_units['genotype'] = 'HOM'
        curr_units['location'] = 'S1'
        curr_units['penetration'] = penet
        allData = allData.append(curr_units,ignore_index=True,sort=False)
        pdf.savefig(f)
        plt.close()
        
    f,ax = plt.subplots(1,1,figsize=(20,20))
    ax.text(0.5,0.5,'S1 WT animals',horizontalalignment='center',verticalalignment='center',fontsize=48)
    pdf.savefig(f)
    plt.close()
    for loc,subj,penet in zip(locs_S1_WT,S1_WT_subjects,S1_WT_location):
        f,units = analyze_mua_by_channel_multistim(loc,show_plot=False,min_z=2.5)
        curr_units = pd.DataFrame(units)
        curr_units['subject'] = subj
        curr_units['genotype'] = 'WT'
        curr_units['location'] = 'S1'
        curr_units['penetration'] = penet
        allData = allData.append(curr_units,ignore_index=True,sort=False)
        pdf.savefig(f)
        plt.close()
    
# VPLM
locs_VPLM_HOM = [r'C:\Data\PGRN_Old\PGRN_HOM_1_Loc4_WhiskStim_C',
        r'C:\Data\PGRN_Old\PGRN_HOM_2_Loc4_WhiskStim_C',
        r'C:\Data\PGRN_Old\PGRN_HOM_3_Loc4_WhiskStim_C',
        r'C:\Data\PGRN_Old\PGRN_HOM_3_Loc6_WhiskStim_C',
        r'C:\Data\PGRN_Old\PGRN_HOM_4_Loc3_WhiskStim_C',
        r'C:\Data\PGRN_Old\PGRN_HOM_4_Loc4_WhiskStim_C',
        r'C:\Data\PGRN_Old\PGRN_HOM_4_Loc6_WhiskStim_C',
        r'C:\Data\PGRN_Old\PGRN_HOM_5_Loc6_WhiskStim_C',]
vplm_HOM_subjects = ['HOM1','HOM2','HOM3','HOM3','HOM4','HOM4','HOM4','HOM5']
vplm_HOM_location = ['Loc4','Loc4','Loc4','Loc6','Loc3','Loc4','Loc6','Loc6']
locs_VPLM_WT = [r'C:\Data\PGRN_Old\PGRN_WT_10_Loc2_whiskerstim_C',
        r'C:\Data\PGRN_Old\PGRN_WT_11_Loc2_whiskerstim_C',
        r'C:\Data\PGRN_Old\PGRN_WT_15_Loc5_WhiskStim_C',
        r'C:\Data\PGRN_Old\PGRN_WT_16_Loc5_WhiskStim_C',
        r'C:\Data\PGRN_Old\PGRN_WT_17_Loc5_WhiskStim_C',
        r'C:\Data\PGRN_Old\PGRN_WT_18_Loc5_WhiskStim_C',
        r'C:\Data\PGRN_Old\PGRN_WT_19_Loc4_WhiskStim_C',]
vplm_WT_subjects = ['WT10','WT11','WT15','WT16','WT17','WT18','WT19']
vplm_WT_location = ['Loc2','Loc2','Loc5','Loc5','Loc5','Loc5','Loc4']
with PdfPages('VPLM_multistim.pdf') as pdf:
    f,ax = plt.subplots(1,1,figsize=(20,20))
    ax.text(0.5,0.5,'VPLM Homozygous animals',horizontalalignment='center',verticalalignment='center',fontsize=48)
    pdf.savefig(f)
    plt.close()
    for loc,subj,penet in zip(locs_VPLM_HOM,vplm_HOM_subjects,vplm_HOM_location):
        f,units = analyze_mua_by_channel_multistim(loc,show_plot=False,min_z=2.5)
        curr_units = pd.DataFrame(units)
        curr_units['subject'] = subj
        curr_units['genotype'] = 'HOM'
        curr_units['location'] = 'VPLM'
        curr_units['penetration'] = penet
        allData = allData.append(curr_units,ignore_index=True,sort=False)
        pdf.savefig(f)
        plt.close()
        
    f,ax = plt.subplots(1,1,figsize=(20,20))
    ax.text(0.5,0.5,'VPLM WT animals',horizontalalignment='center',verticalalignment='center',fontsize=48)
    pdf.savefig(f)
    plt.close()
    for loc,subj,penet in zip(locs_S1_WT,vplm_WT_subjects,vplm_WT_location):
        f,units = analyze_mua_by_channel_multistim(loc,show_plot=False,min_z=2.5)
        curr_units = pd.DataFrame(units)
        curr_units['subject'] = subj
        curr_units['genotype'] = 'WT'
        curr_units['location'] = 'VPLM'
        curr_units['penetration'] = penet
        allData = allData.append(curr_units,ignore_index=True,sort=False)
        pdf.savefig(f)
        plt.close()

allData.to_pickle('AllUnitData_MUA_multistim.pickle')