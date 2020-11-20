from util.whisker_analyses import analyze_mua_by_channel_multistim,plot_autocorr_by_channel
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
        curr_row = 0
        for row in sessions.itertuples():
            print('{0} of {1}::{2}'.format(curr_row,len(sessions),row.tank_name))
            curr_row +=1
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
            session_path=os.path.join(row.analysis_data_path,row.tank_name)
            
            f,units = analyze_mua_by_channel_multistim(session_path,show_plot=False,min_z=0,stim_time='timestamps_R.np',common_details=details,chans=chan_that_session)
            if units:
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
                if units:
                    f.suptitle(row.tank_name+' IPSI stim',fontsize=20)
                    curr_units=pd.DataFrame(units)
                    data_df=data_df.append(curr_units,ignore_index=True,sort=False)
                    pdf.savefig(f)
                    plt.close()        
    return data_df

def get_autocorrelation_function(sessions,base_path):
    curr_row = 0
    all_units = []
    with PdfPages(os.path.join(base_path,'Autocorr.pdf')) as pdf:
        for row in sessions.itertuples():
            print('{0} of {1}::{2}'.format(curr_row,len(sessions),row.tank_name))
            curr_row +=1
            session_path=os.path.join(row.analysis_data_path,row.tank_name)
            f,units = plot_autocorr_by_channel(session_path,row)
            f.suptitle(row.tank_name,fontsize=20)
            pdf.savefig(f)
            plt.close()
            
            if not all_units:
                all_units = units
            else:
                all_units.extend(units)
    return all_units
        
        
if __name__=='__main__':
    # tdt_tanks=tdt_tanks_all
    # getters=get_tdt_subject,get_tdt_date,get_tdt_time,get_tdt_genotype,get_repeat_number
    
    if os.name=='nt':
        base_path=r'C:\Users\bsriram\Desktop\Code\neuralcircuits_analysis\Results'
    else:
        base_path='/home/bsriram/code/neuralcircuits_analysis'
    
    which_tanks = tank_details.loc[tank_details.n_stim==0]
    units = get_autocorrelation_function(which_tanks,base_path)
    df = pd.DataFrame(units)
    df.to_pickle(os.path.join(base_path,'Autocorr_df.pickle'))
    