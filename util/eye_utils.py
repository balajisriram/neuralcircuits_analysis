import os
import time
import numpy as np
import ffmpy
import platform

def convert_video(video_file,save_path):
    if platform.system()=='Windows':
        executable_loc = r'C:\ffmpeg-20190204-075fd5b-win64-static\bin\ffmpeg.exe'
        assert os.path.exists(executable_loc),'unable to find the executable for ffmpeg'
    elif platform.system()=='Linux':
        executable_loc = 'ffmpeg'        
    new_filename = 'eye_trk.mp4'
    new_file_path = os.path.join(save_path,new_filename)
    if not os.path.isdir(save_path):
        os.makedirs(save_path)
    job = ffmpy.FFmpeg(executable=executable_loc,inputs={video_file:None},outputs={new_file_path:'-vf scale=320:-1 -crf 10'})
    job.run()
    
    return new_file_path
    
def make_eye_analysis_script(sess,save_loc,base_on_server='/home/bsriram/data_biogen/'):
    subj,date,sess_folder,rig,running,hs,elec = sess
    
    if not base_on_server: 
        wd=save_loc
    else:
        sep = '/'
        dat_path,trk_dir = os.path.split(save_loc)
        base,dat_folder = os.path.split(dat_path)
        wd = base_on_server+dat_folder+sep+trk_dir
        # wd = os.path.join(base_on_server,dat_folder,trk_dir)
        
    with open(os.path.join(save_loc,'track_eye.py'),'w') as f:
        f.write("""import deeplabcut
import tensorflow as tf
import os
import time
path_config_file = '/home/bsriram/code/eye-trk/examples/Eye-Tracking-BAS-2019-02-17/config.yaml'
os.environ["DLClight"]="True"
videos = ['{0}']
t_start = time.time()
print('analyzing')
deeplabcut.analyze_videos(path_config_file,videos)
print('creating labeled_video')
deeplabcut.create_labeled_video(path_config_file,videos)
print("It took ",time.time()-t_start)""".format(wd+sep+'eye_trk.mp4'))
    print('wrote into '+os.path.join(save_loc,'track_eye.py'))

def make_eye_trk_qsub_script(sess,save_loc,base_on_server='/home/bsriram/data_biogen/'):
    subj,date,sess_folder,rig,running,hs,elec = sess
    if not base_on_server: 
        wd=save_loc
        sep = '\\'
    else:
        sep = '/'
        dat_path,trk_dir = os.path.split(save_loc)
        base,dat_folder = os.path.split(dat_path)
        wd = base_on_server+dat_folder+sep+trk_dir
        # wd = os.path.join(base_on_server,dat_folder,trk_dir)
    
    with open(os.path.join(save_loc,'submit.qsub'),'w') as f:
        f.write("""#!/bin/bash
#$ -N {0}
#$ -l h_rt=24:00:00 
#$ -q long.q 
#$ -wd {1}
#$ -j no 
#$ -M balaji.sriram@biogen.com
#$ -m be 
#$ -e error.log 
#$ -o output.log 
#$ -pe openmpi-fillup 1

##########################################################################
# Start your script here #
##########################################################################
# Load the modules you need.
source /home/bsriram/miniconda3/bin/activate dlc_cpu
# Run some commands.
python track_eye.py

# Exit successfully.
exit 0
""".format(subj+'_'+date+'_eyetrk',wd))
    print('wrote into '+os.path.join(save_loc,'submit.qsub'))
        
def save_eye_track_video(sess,src_loc,sav_loc):
    subj,date,sess_folder,rig,running,hs,elec = sess
 
    # any other sibject in that session?
    if rig=='R1':avi_ending='_1.AVI'
    if rig=='R2':avi_ending='_2.AVI'
    avi_file = []
    # try searching for endswith(_1.AVI)
    avi_file = [f for f in os.listdir(src_loc) if f.endswith(avi_ending)]
    # if not found, then some dont have the _rig.AVI if only R1 is run
    if not avi_file and rig=='R1':
        avi_file = [f for f in os.listdir(src_loc) if f.endswith('1.AVI')]
        
    if len(avi_file)==0:
        print(subj,'\t:',sess_folder,'\t:','no eye tracking found. continuing')
        save_path = []
    elif len(avi_file)>1:
        print(subj,'\t:',sess_folder,'\t:','too many tracking found. ')
    else:
        video_file_path = os.path.join(src_loc,avi_file[0])
        save_path = os.path.join(sav_loc,'eye_tracking')
        print(subj,'\t:',sess_folder,'\t:','save '+video_file_path+' in '+save_path)
        convert_video(video_file_path,save_path)
        
    return save_path