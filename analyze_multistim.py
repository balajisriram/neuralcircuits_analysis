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
    if os.name=='nt':
        base_path=r'C:\Users\bsriram\Desktop\Code\neuralcircuits_analysis\Results'
    else:
        base_path='/home/bsriram/code/neuralcircuits_analysis'
    which_tanks = tank_details.loc[tank_details.n_stim>0]
    data_df = analyze_and_make_pdf(base_path = base_path, sessions=which_tanks, output_pdf='trial.pdf',title='trial')
    data_df.to_pickle(os.path.join(base_path,'output.pickle')