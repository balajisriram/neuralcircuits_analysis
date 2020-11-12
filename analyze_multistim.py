from util.whisker_analyses import analyze_mua_by_channel_multistim
from util.CohortDetails import tank_df as tank_details
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import pickle
import pandas as pd
# from tdt_utils import get_tdt_date, get_tdt_time, get_tdt_subject, get_tdt_genotype, has_ipsi_stim, get_repeat_number
import os
import shutil


def analyze_and_make_pdf(base_path=None, sessions = [], output_pdf='output.pdf',title='Unknown',data_df=None, chan_list=None):
    if not data_df:
        data_df=pd.DataFrame()
    with PdfPages(os.path.join(base_path,output_pdf)) as pdf:
        # title page
        f,ax=plt.subplots(1,1,figsize=(20,20))
        ax.text(0.5,0.5,title,horizontalalignment='center',verticalalignment='center',fontsize=48)
        pdf.savefig(f)
        
        for row in sessions.itertuples():
            # get_subject,get_date,get_time,get_genotype,get_repeat = getters
            details={}
            details['subject_id'] = row.subject
            details['tank_name'] = row.tank_name
            details['date'] = row.date
            details['time'] = row.time
            details['genotype'] = row.genotype
            details['repeat'] = row.repeat
            details['drug_state'] = row.drug_details
            if row.stim_lateralism=='bilateral':
                details['has_ipsi'] = True
            else:
                details['has_ipsi'] = False
            # if chan_list:
                # chan_that_session = chan_list[details['subject_id']]
                # print(details['subject_id'],':',chan_that_session)
            chan_that_session = range(16)
            details['stim_location'] = 'contra'
            # if base_path: session_path=os.path.join(base_path,row.tank_name)
            session_path=os.path.join(row.raw_data_path,row.tank_name)
            
            f,units = analyze_mua_by_channel_multistim(session_path,show_plot=False,min_z=0,stim_time='timestamps_R.np',common_details=details,chans=chan_that_session)
            f.suptitle(row.tank_name+' CONTRA stim',fontsize=20)
            curr_units=pd.DataFrame(units)
            data_df=data_df.append(curr_units,ignore_index=True,sort=False)
            pdf.savefig(f)
            plt.close()
            
            if row.stim_lateralism=='bilateral':
                details={}
                details['subject_id'] = row.subject
                details['tank_name'] = row.tank_name
                details['date'] = row.date
                details['time'] = row.time
                details['genotype'] = row.genotype
                details['repeat'] = row.repeat
                details['drug_state'] = row.drug_details
                details['has_ipsi'] = True
                details['stim_location'] = 'ipsi'
                
                f,units = analyze_mua_by_channel_multistim(session_path,show_plot=False,min_z=0,stim_time='timestamps_L.np',common_details=details)
                f.suptitle(row.tank_name+' IPSI stim',fontsize=20)
                curr_units=pd.DataFrame(units)
                data_df=data_df.append(curr_units,ignore_index=True,sort=False)
                
                pdf.savefig(f)
                plt.close()        
    return data_df
 
if __name__=='__main__':

    chan_list = {
    'PGRN_369': [ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 13],
    'PGRN_370': [ 4,  5,  6,  7,  8,  9, 12],
    'PGRN_371': [ 1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14],
    'PGRN_372': [ 0,  1,  4,  5,  6,  7,  9, 10], 
    'PGRN_374': [ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13],
    'PGRN_375': [ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 13, 14, 15],
    'PGRN_376': [ 3,  4,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15], 
    'PGRN_377': [ 0,  3,  4,  5,  6,  7,  8,  9, 13],
    'PGRN_378': [ 0,  1,  2,  3,  4,  8, 11, 12, 13, 14],
    'PGRN_379': [ 0,  1,  2,  3,  5,  6, 11],
    'PGRN_380': [ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10],
    'PGRN_381': [ 0,  1,  2,  3,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14],
    'PGRN_382': [ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10],
    'PGRN_383': [ 3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15],
    'PGRN_384': [ 0,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13],
    'PGRN_388': [ 0,  1,  2,  3,  4,  5,  6,  7, 10, 11, 14],
    'PGRN_389': [ 1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11],
    'PGRN_390': [ 0,  6,  7,  8,  9, 10, 11],
    'PGRN_391': [ 1,  2,  8,  9, 10, 13, 15],
    'PGRN_393': [9],
    'PGRN_396': [ 2,  3,  4,  5,  6,  7,  8,  9, 10, 11],
    'PGRN_400': [ 0,  1, 12],
    'PGRN_401': [ 5,  7, 12],}
    # tdt_tanks=tdt_tanks_all
    # getters=get_tdt_subject,get_tdt_date,get_tdt_time,get_tdt_genotype,get_repeat_number
    base_path='/home/bsriram/code/neuralcircuits_analysis'
    which_tanks = tank_details.loc[tank_details.tank_name=='PGRN_390_HET-200721-125503']
    data_df = analyze_and_make_pdf(base_path = base_path, sessions=which_tanks, output_pdf='trial.pdf',title='trial')
    data_df.to_pickle(os.path.join(base_path,'output.pickle'))

    
    
if False:
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