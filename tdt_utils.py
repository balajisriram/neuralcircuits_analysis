import tdt
import os
import numpy as np

from util.whisker_analyses import analyze_mua_by_channel_multistim

def load_tdt_tank(dirname=r'C:\Data\PGRN\PGRN_155_WT-200401-154248'):
    temp = tdt.read_block(dirname)
    data = temp['streams']['BBLF']['data']
    ratio  = np.floor(32767/data.max())
    data = np.int16(data*ratio)
    return data.transpose(), ratio
    
def pack(data,path,filename='raw.dat',dref=None,sr = np.array(24414.0625), ratio=None):
    data.tofile(os.path.join(path,filename))
    ratio.tofile(os.path.join(path,'voltage_scale.np'))
    sr.tofile(os.path.join(path,'sample_rate.np'))

def pack_timestamps(path,):
    temp = tdt.read_block(path)
    onset = temp['epocs']['TRG1']['onset']
    offset = temp['epocs']['TRG1']['offset']
    timediff = np.median(offset-onset)
    onset.tofile(os.path.join(path,'timestamps.np'))
    timediff.tofile(os.path.join(path,'timediff.np'))
    
def pack_timestamps_dual_whisker_and_cams(path,):
    temp = tdt.read_block(path)
    
    try:
        onset_L = temp['epocs']['TRGL']['onset']
        offset_L = temp['epocs']['TRGL']['offset']
        timediff_L = np.median(offset_L-onset_L)
        onset_L.tofile(os.path.join(path,'timestamps_L.np'))
        timediff_L.tofile(os.path.join(path,'timediff_L.np'))
    except AttributeError:
        print('No left stimulus found')
    
    onset_R = temp['epocs']['TRGR']['onset']
    offset_R = temp['epocs']['TRGR']['offset']
    timediff_R = np.median(offset_R-onset_R)
    onset_R.tofile(os.path.join(path,'timestamps_R.np'))
    timediff_R.tofile(os.path.join(path,'timediff_R.np'))
    
    try:
        cam_onset = temp['epocs']['PtC0']['onset']
        cam_offset = temp['epocs']['PtC0']['offset']
        timediff_cam = np.median(cam_offset-cam_onset)
        cam_onset.tofile(os.path.join(path,'timestamps_WhiskerCam.np'))
        timediff_cam.tofile(os.path.join(path,'exposure_WhiskerCam.np'))
    except AttributeError:
        print('Whisker Camera not found')
        
    try:
        cam_onset = temp['epocs']['PtC1']['onset']
        cam_offset = temp['epocs']['PtC1']['offset']
        timediff_cam = np.median(cam_offset-cam_onset)
        cam_onset.tofile(os.path.join(path,'timestamps_PupilCam.np'))
        timediff_cam.tofile(os.path.join(path,'exposure_PupilCam.np'))
    except AttributeError:
        print('Pupil Camera not found')
    
awake_cohort1_tanks = ['PGRN_155_WT-200330-130817','PGRN_155_WT-200401-154248','PGRN_155_WT-200401-155139','PGRN_155_WT-200428-142415','PGRN_155_WT-200428-143352','PGRN_155_WT-200501-144905','PGRN_155_WT-200501-150247','PGRN_158_WT-200330-134401','PGRN_158_WT-200401-162227','PGRN_158_WT-200401-163127','PGRN_158_WT-200428-150301','PGRN_158_WT-200428-151153','PGRN_158_WT-200501-155657','PGRN_158_WT-200501-161040','PGRN_175_HET-200330-125152','PGRN_175_HET-200401-152118','PGRN_175_HET-200401-153033','PGRN_175_HET-200428-140532','PGRN_175_HET-200428-141408','PGRN_175_HET-200501-141903','PGRN_175_HET-200501-143232','PGRN_181_HET-200330-132613','PGRN_181_HET-200401-160230','PGRN_181_HET-200401-161132','PGRN_181_HET-200428-144511','PGRN_181_HET-200428-145342','PGRN_181_HET-200501-152447','PGRN_181_HET-200501-153924','PGRN_182_WT-200330-123844','PGRN_182_WT-200401-150003','PGRN_182_WT-200401-151015','PGRN_182_WT-200428-135320','PGRN_182_WT-200501-134811','PGRN_182_WT-200501-140216','PGRN_183_HOM-200401-164914','PGRN_183_HOM-200401-165806','PGRN_183_HOM-200428-133435','PGRN_183_HOM-200428-134317','PGRN_183_HOM-200501-131710','PGRN_183_HOM-200501-133051','PGRN_184_HOM-200330-121157','PGRN_184_HOM-200330-121906','PGRN_184_HOM-200401-144100','PGRN_184_HOM-200401-170826','PGRN_184_HOM-200401-171613','PGRN_184_HOM-200428-131207','PGRN_184_HOM-200428-132202','PGRN_184_HOM-200501-123637','PGRN_184_HOM-200501-125015','PGRN_HT_175_LMRS2-200303-144850','PGRN_WT_158_RS-200302-151344']

awake_cohort2_tanks = ['PGRN_369_WT-200713-140327','PGRN_370_HOM-200708-131420','PGRN_371_WT-200713-134124','PGRN_372_HET-200708-165326','PGRN_374_HOM-200713-102856','PGRN_375_HOM-200713-105801','PGRN_376_WT-200708-124614','PGRN_377_HET-200708-155016','PGRN_378_WT-200708-150655','PGRN_379_HET-200713-111941','PGRN_380_WT-200708-161417','PGRN_381_HET-200713-130234','PGRN_382_WT-200708-163354','PGRN_383_HOM-200708-143833','PGRN_384_WT-200706-164212','PGRN_384_WT-200715-123104','PGRN_388_HOM-200713-123705','PGRN_389_WT-200706-161140','PGRN_389_WT-200715-114822','PGRN_390_HET-200706-154900','PGRN_390_HET-200715-112446','PGRN_391_HOM-200706-163047','PGRN_391_HOM-200715-120933','PGRN_393_HOM-200713-145953','PGRN_396_HOM-200713-142749','PGRN_400_HET-200708-133938','PGRN_401_HET-200708-140920',]

awake_cohort2_tanks2=['PGRN_370_HOM-200721-140639','PGRN_370_HOM-200721-141333','PGRN_376_WT-200721-134456','PGRN_376_WT-200721-135218','PGRN_377_HET-200721-143537','PGRN_377_HET-200721-144337','PGRN_378_WT-200721-145141','PGRN_378_WT-200721-145931','PGRN_384_WT-200721-132153','PGRN_384_WT-200721-132927','PGRN_389_WT-200721-122505','PGRN_389_WT-200721-123728','PGRN_390_HET-200721-124605','PGRN_390_HET-200721-125503','PGRN_391_HOM-200721-130354','PGRN_391_HOM-200721-131154','PGRN_400_HET-200721-142125','PGRN_400_HET-200721-142857','PGRN_401_HET-200721-150644','PGRN_401_HET-200721-151431',]

awake_cohort2_tanks3=['PGRN_369_WT-200723-132704','PGRN_369_WT-200723-133516','PGRN_371_WT-200723-131157','PGRN_371_WT-200723-131853','PGRN_372_HET-200723-110143','PGRN_372_HET-200723-111033','PGRN_374_HOM-200723-112530','PGRN_374_HOM-200723-113210','PGRN_375_HOM-200723-114142','PGRN_375_HOM-200723-114840','PGRN_379_HET-200723-115754','PGRN_379_HET-200723-120738','PGRN_380_WT-200723-102746','PGRN_380_WT-200723-103602','PGRN_381_HET-200723-125642','PGRN_381_HET-200723-130329','PGRN_382_WT-200723-104510','PGRN_382_WT-200723-105226','PGRN_383_HOM-200723-123322','PGRN_383_HOM-200723-124229','PGRN_388_HOM-200723-121521','PGRN_388_HOM-200723-122209','PGRN_393_HOM-200723-140955','PGRN_393_HOM-200723-141642','PGRN_396_HOM-200723-134711','PGRN_396_HOM-200723-135345',]

if __name__=='__main__':
    base_path = r'C:\Users\bsriram\Desktop\Data\PGRN_more'
    tanks = awake_cohort2_tanks3
        
    for tank in tanks:
        tank_path = os.path.join(base_path,tank)
        data,ratio = load_tdt_tank(tank_path)
        pack(data,tank_path,ratio=ratio)
        pack_timestamps_dual_whisker_and_cams(tank_path)
        print('done for {0}'.format(tank))
        print()
        print()
        #analyze_mua_by_channel_multistim(tank_path,type='raster',show_plot=True,min_z=2,z_zoom=2,interval=[-1,3])