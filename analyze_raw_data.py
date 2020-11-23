from util.CohortDetails import tank_df as tank_details
import matplotlib.pyplot as plt
import pickle
import pandas as pd
import os
import shutil
import tdt
import numpy as np

def analyze_raw(sessions):
    curr_row = 0
    all_units = []
    for row in sessions.itertuples():
        print('{0} of {1}::{2}'.format(curr_row,len(sessions),row.tank_name))
        curr_row +=1
        # get_subject,get_date,get_time,get_genotype,get_repeat = getters
        that_session={}
        session_path=os.path.join(row.raw_data_path,row.tank_name)
        temp = tdt.read_block(session_path)
        data = temp['streams']['BBLF']['data']
        data_std = np.std(data,axis=1)
        for channel in range(16):
            this_unit = {}
            this_unit['unit_and_tank_id'] = row.tank_name+'_'+str(channel)
            this_unit['raw_signal_std'] = data_std[channel]
            all_units.append(this_unit)
        
    return all_units
 
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
    units = analyze_raw(tank_details)
    data_df=pd.DataFrame(units)
    data_df.to_pickle(os.path.join(base_path,'rawData_std.pickle'))
    
