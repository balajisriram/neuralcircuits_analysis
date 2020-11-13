from util.CohortDetails import tank_df as tank_details
import matplotlib.pyplot as plt
import pickle
import pandas as pd
import os
import shutil
import tdt


def analyze_raw(base_path=None, sessions = [], output_pdf='output.pdf',title='Unknown',data_df=None, chan_list=None):
    if not data_df:
        data_df=pd.DataFrame()
    curr_row = 0
    for row in sessions.itertuples():
        print('{0} of {1}::{2}'.format(curr_row,len(sessions),row.tank_name))
        curr_row +=1
        # get_subject,get_date,get_time,get_genotype,get_repeat = getters
        that_session={}
        session_path=os.path.join(row.raw_data_path,row.tank_name)
        
        
        f,units = analyze_mua_by_channel_multistim(session_path,show_plot=False,min_z=0,stim_time='timestamps_R.np',common_details=details,chans=chan_that_session)
        f.suptitle(row.tank_name+' CONTRA stim',fontsize=20)
        curr_units=pd.DataFrame(units)
        data_df=data_df.append(curr_units,ignore_index=True,sort=False)
        pdf.savefig(f)
        plt.close()
        
    return data_df
 
def assert_raw_exists(sessions):
    curr_row = 0
    for row in sessions.itertuples():
        # get_subject,get_date,get_time,get_genotype,get_repeat = getters
        that_session={}
        session_path=os.path.join(row.raw_data_path,row.tank_name)
        sev_files = [f for f in os.listdir(session_path) if f.endswith('sev')]
        print('{0} of {1}::{2} contains {3} sev files'.format(curr_row,len(sessions),row.tank_name,len(sev_files)))
        curr_row +=1
        
 
if __name__=='__main__':
    
    if os.name=='nt':
        base_path=r'C:\Users\bsriram\Desktop\Code\neuralcircuits_analysis\Results'
    else:
        base_path='/home/bsriram/code/neuralcircuits_analysis'
    
    which_tanks = tank_details.loc[tank_details.n_stim>0]
    assert_raw_exists(which_tanks)
    