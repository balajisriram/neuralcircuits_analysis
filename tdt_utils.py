import tdt
import os
import numpy as np
import os.path as op
import shutil
import datetime

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
    
    try:
        onset_R = temp['epocs']['TRGR']['onset']
        offset_R = temp['epocs']['TRGR']['offset']
        timediff_R = np.median(offset_R-onset_R)
        onset_R.tofile(os.path.join(path,'timestamps_R.np'))
        timediff_R.tofile(os.path.join(path,'timediff_R.np'))
    except AttributeError:
        print('No right stimulus found')
        
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

def copy_analyses(src=r'C:\Users\bsriram\Desktop\Data\PGRN_Coh2',des=r'C:\Users\bsriram\Desktop\Data\PGRN_Coh2_ephys_analysis',tanks=None, types=['.np','.prb','.prm','.kwik','.kwx',]):
    opj=os.path.join
    if not tanks: 
        tanks = [f for f in os.listdir(src) if os.path.isdir(opj(src,f))]
    for tank in tanks:
        print(tank)
        src_folder = opj(src,tank)
        des_folder = opj(des,tank)
        if not os.path.exists(des_folder): os.mkdir(des_folder)
        for type in types:
            files = [f for f in os.listdir(src_folder) if f.endswith(type)]
            for file in files:
                shutil.copyfile(opj(src_folder,file),opj(des_folder,file))
                print('    ',file,'...DONE')
        print()
        print()
        print()
    
def copy_prm_prb(src_folder=r'\\camhpcisixcifs.biogen.com\dept\electrophysiology\invivo_mouse\PGRN\PGRN_AwakeCohort2\DataAnalyzed\PGRN_369_WT-200713-140327',des=r'C:\Users\bsriram\Desktop\Data\PGRN_Coh2',tanks=None,types=['.prm','.prb']):
    opj = os.path.join
    if not tanks: tanks = [f for f in os.listdir(des) if os.path.isdir(opj(des,f))]
    for tank in tanks:
        print(tank)
        des_folder = opj(des,tank)
        for type in types:
            files = [f for f in os.listdir(src_folder) if f.endswith(type)]
            for file in files:
                shutil.copyfile(opj(src_folder,file),opj(des_folder,file))
                print('    ',file,'...DONE')
        print()
        print()
        print()

def get_subjects_in_tank_list(tanks):
    subjects = np.asarray([t.split('-')[0] for t in tanks])
    return np.unique(subjects)
    
def get_tank_sequencing(tanks):
    orig_list = tanks
    running_list = orig_list
    
    subjects = get_subjects_in_tank_list(tanks)
    pre_tank = []
    post_tank1 = []
    post_tank2 = []
    for subject in subjects:
        sess_time = sorted([datetime.datetime.strptime('{0}-{1}'.format(t.split('-')[1],t.split('-')[2]),'%y%m%d-%H%M%S') for t in tanks if subject in t])
        pre_tank.append('{0}-{1}'.format(subject,datetime.datetime.strftime(sess_time[0],'%y%m%d-%H%M%S')))
        post_tank1.append('{0}-{1}'.format(subject,datetime.datetime.strftime(sess_time[1],'%y%m%d-%H%M%S')))
        post_tank2.append('{0}-{1}'.format(subject,datetime.datetime.strftime(sess_time[2],'%y%m%d-%H%M%S')))
    return pre_tank,post_tank1,post_tank2


awake_cohort1_tanks = ['PGRN_155_WT-200330-130817','PGRN_155_WT-200401-154248','PGRN_155_WT-200401-155139','PGRN_155_WT-200428-142415','PGRN_155_WT-200428-143352','PGRN_155_WT-200501-144905','PGRN_155_WT-200501-150247','PGRN_158_WT-200330-134401','PGRN_158_WT-200401-162227','PGRN_158_WT-200401-163127','PGRN_158_WT-200428-150301','PGRN_158_WT-200428-151153','PGRN_158_WT-200501-155657','PGRN_158_WT-200501-161040','PGRN_175_HET-200330-125152','PGRN_175_HET-200401-152118','PGRN_175_HET-200401-153033','PGRN_175_HET-200428-140532','PGRN_175_HET-200428-141408','PGRN_175_HET-200501-141903','PGRN_175_HET-200501-143232','PGRN_181_HET-200330-132613','PGRN_181_HET-200401-160230','PGRN_181_HET-200401-161132','PGRN_181_HET-200428-144511','PGRN_181_HET-200428-145342','PGRN_181_HET-200501-152447','PGRN_181_HET-200501-153924','PGRN_182_WT-200330-123844','PGRN_182_WT-200401-150003','PGRN_182_WT-200401-151015','PGRN_182_WT-200428-135320','PGRN_182_WT-200501-134811','PGRN_182_WT-200501-140216','PGRN_183_HOM-200401-164914','PGRN_183_HOM-200401-165806','PGRN_183_HOM-200428-133435','PGRN_183_HOM-200428-134317','PGRN_183_HOM-200501-131710','PGRN_183_HOM-200501-133051','PGRN_184_HOM-200330-121157','PGRN_184_HOM-200330-121906','PGRN_184_HOM-200401-144100','PGRN_184_HOM-200401-170826','PGRN_184_HOM-200401-171613','PGRN_184_HOM-200428-131207','PGRN_184_HOM-200428-132202','PGRN_184_HOM-200501-123637','PGRN_184_HOM-200501-125015','PGRN_HT_175_LMRS2-200303-144850','PGRN_WT_158_RS-200302-151344']

awake_cohort2_tanks = ['PGRN_369_WT-200713-140327','PGRN_370_HOM-200708-131420','PGRN_371_WT-200713-134124','PGRN_372_HET-200708-165326','PGRN_374_HOM-200713-102856','PGRN_375_HOM-200713-105801','PGRN_376_WT-200708-124614','PGRN_377_HET-200708-155016','PGRN_378_WT-200708-150655','PGRN_379_HET-200713-111941','PGRN_380_WT-200708-161417','PGRN_381_HET-200713-130234','PGRN_382_WT-200708-163354','PGRN_383_HOM-200708-143833','PGRN_384_WT-200706-164212','PGRN_384_WT-200715-123104','PGRN_388_HOM-200713-123705','PGRN_389_WT-200706-161140','PGRN_389_WT-200715-114822','PGRN_390_HET-200706-154900','PGRN_390_HET-200715-112446','PGRN_391_HOM-200706-163047','PGRN_391_HOM-200715-120933','PGRN_393_HOM-200713-145953','PGRN_396_HOM-200713-142749','PGRN_400_HET-200708-133938','PGRN_401_HET-200708-140920',]

awake_cohort2_tanks2=['PGRN_370_HOM-200721-140639','PGRN_370_HOM-200721-141333','PGRN_376_WT-200721-134456','PGRN_376_WT-200721-135218','PGRN_377_HET-200721-143537','PGRN_377_HET-200721-144337','PGRN_378_WT-200721-145141','PGRN_378_WT-200721-145931','PGRN_384_WT-200721-132153','PGRN_384_WT-200721-132927','PGRN_389_WT-200721-122505','PGRN_389_WT-200721-123728','PGRN_390_HET-200721-124605','PGRN_390_HET-200721-125503','PGRN_391_HOM-200721-130354','PGRN_391_HOM-200721-131154','PGRN_400_HET-200721-142125','PGRN_400_HET-200721-142857','PGRN_401_HET-200721-150644','PGRN_401_HET-200721-151431',]

awake_cohort2_tanks3=['PGRN_369_WT-200723-132704','PGRN_369_WT-200723-133516','PGRN_371_WT-200723-131157','PGRN_371_WT-200723-131853','PGRN_372_HET-200723-110143','PGRN_372_HET-200723-111033','PGRN_374_HOM-200723-112530','PGRN_374_HOM-200723-113210','PGRN_375_HOM-200723-114142','PGRN_375_HOM-200723-114840','PGRN_379_HET-200723-115754','PGRN_379_HET-200723-120738','PGRN_380_WT-200723-102746','PGRN_380_WT-200723-103602','PGRN_381_HET-200723-125642','PGRN_381_HET-200723-130329','PGRN_382_WT-200723-104510','PGRN_382_WT-200723-105226','PGRN_383_HOM-200723-123322','PGRN_383_HOM-200723-124229','PGRN_388_HOM-200723-121521','PGRN_388_HOM-200723-122209','PGRN_393_HOM-200723-140955','PGRN_393_HOM-200723-141642','PGRN_396_HOM-200723-134711','PGRN_396_HOM-200723-135345',]

awake_cohort2_tanks4=['PGRN_369_WT-200820-114036','PGRN_369_WT-200820-114919','PGRN_370_HOM-200818-123855','PGRN_370_HOM-200818-124535','PGRN_371_WT-200820-112534','PGRN_371_WT-200820-113158','PGRN_372_HET-200818-142954','PGRN_372_HET-200818-143637','PGRN_374_HOM-200818-144524','PGRN_374_HOM-200818-145210','PGRN_375_HOM-200818-150015','PGRN_375_HOM-200818-150645','PGRN_376_WT-200818-122154','PGRN_376_WT-200818-122812','PGRN_377_HET-200818-133837','PGRN_377_HET-200818-134502','PGRN_378_WT-200818-132329','PGRN_378_WT-200818-133007','PGRN_379_HET-200818-151436','PGRN_379_HET-200818-152104','PGRN_380_WT-200818-135417','PGRN_380_WT-200818-140105','PGRN_381_HET-200820-111010','PGRN_381_HET-200820-111634','PGRN_382_WT-200818-140937','PGRN_382_WT-200818-141607','PGRN_384_WT-200818-114324','PGRN_384_WT-200818-114953','PGRN_388_HOM-200818-153425','PGRN_388_HOM-200818-154112','PGRN_389_WT-200818-110750','PGRN_389_WT-200818-111418','PGRN_390_HET-200818-105318','PGRN_390_HET-200818-110039','PGRN_391_HOM-200818-112909','PGRN_391_HOM-200818-113530','PGRN_393_HOM-200820-121424','PGRN_393_HOM-200820-122041','PGRN_396_HOM-200820-120015','PGRN_396_HOM-200820-120631','PGRN_400_HET-200818-125348','PGRN_400_HET-200818-130007','PGRN_401_HET-200818-130825','PGRN_401_HET-200818-131443',]

awake_cohort2_tanks5=['PGRN_369_WT-200925-170556','PGRN_369_WT-200925-171323','PGRN_369_WT-200925-172542','PGRN_371_WT-200925-162655','PGRN_371_WT-200925-163541','PGRN_371_WT-200925-164756','PGRN_372_HET-200925-114651','PGRN_372_HET-200925-115922','PGRN_372_HET-200925-121226','PGRN_374_HOM-200925-123008','PGRN_374_HOM-200925-123753','PGRN_374_HOM-200925-125104','PGRN_375_HOM-200925-131008','PGRN_375_HOM-200925-131704','PGRN_375_HOM-200925-133107','PGRN_376_WT-200923-164027','PGRN_376_WT-200923-164823','PGRN_376_WT-200923-170108','PGRN_377_HET-200924-155214','PGRN_377_HET-200924-155954','PGRN_377_HET-200924-161222','PGRN_378_WT-200924-151011','PGRN_378_WT-200924-152045','PGRN_378_WT-200924-153304','PGRN_379_HET-200925-143015','PGRN_379_HET-200925-143852','PGRN_379_HET-200925-145125','PGRN_380_WT-200924-162939','PGRN_380_WT-200924-163722','PGRN_380_WT-200924-165826','PGRN_381_HET-200925-154922','PGRN_381_HET-200925-155742','PGRN_381_HET-200925-160959','PGRN_382_WT-200924-172207','PGRN_382_WT-200924-173146','PGRN_382_WT-200924-174422','PGRN_384_WT-200923-160011','PGRN_384_WT-200923-160859','PGRN_384_WT-200923-162127','PGRN_388_HOM-200925-151039','PGRN_388_HOM-200925-151821','PGRN_388_HOM-200925-153040','PGRN_389_WT-200923-143845','PGRN_389_WT-200923-145259','PGRN_389_WT-200923-150634','PGRN_391_HOM-200923-152335','PGRN_391_HOM-200923-153141','PGRN_391_HOM-200923-154439','PGRN_393_HOM-200925-181822','PGRN_393_HOM-200925-182617','PGRN_393_HOM-200925-183840','PGRN_396_HOM-200925-174223','PGRN_396_HOM-200925-174956','PGRN_396_HOM-200925-180218','PGRN_400_HET-200924-132658','PGRN_400_HET-200924-133908','PGRN_400_HET-200924-135146','PGRN_401_HET-200924-143215','PGRN_401_HET-200924-144141','PGRN_401_HET-200924-145359',]

awake_cohort2_pre_drug = ['PGRN_369_WT-200925-170556', 'PGRN_371_WT-200925-162655', 'PGRN_372_HET-200925-114651', 'PGRN_374_HOM-200925-123008', 'PGRN_375_HOM-200925-131008', 'PGRN_376_WT-200923-164027', 'PGRN_377_HET-200924-155214', 'PGRN_378_WT-200924-151011', 'PGRN_379_HET-200925-143015', 'PGRN_380_WT-200924-162939', 'PGRN_381_HET-200925-154922', 'PGRN_382_WT-200924-172207', 'PGRN_384_WT-200923-160011', 'PGRN_388_HOM-200925-151039', 'PGRN_389_WT-200923-143845', 'PGRN_391_HOM-200923-152335', 'PGRN_393_HOM-200925-181822', 'PGRN_396_HOM-200925-174223', 'PGRN_400_HET-200924-132658', 'PGRN_401_HET-200924-143215']

awake_cohort2_post_drug = ['PGRN_369_WT-200925-171323', 'PGRN_371_WT-200925-163541', 'PGRN_372_HET-200925-115922', 'PGRN_374_HOM-200925-123753', 'PGRN_375_HOM-200925-131704', 'PGRN_376_WT-200923-164823', 'PGRN_377_HET-200924-155954', 'PGRN_378_WT-200924-152045', 'PGRN_379_HET-200925-143852', 'PGRN_380_WT-200924-163722', 'PGRN_381_HET-200925-155742', 'PGRN_382_WT-200924-173146', 'PGRN_384_WT-200923-160859', 'PGRN_388_HOM-200925-151821', 'PGRN_389_WT-200923-145259', 'PGRN_391_HOM-200923-153141', 'PGRN_393_HOM-200925-182617', 'PGRN_396_HOM-200925-174956', 'PGRN_400_HET-200924-133908', 'PGRN_401_HET-200924-144141']

awake_cohort2_post_drug2 = ['PGRN_369_WT-200925-172542', 'PGRN_371_WT-200925-164756', 'PGRN_372_HET-200925-121226', 'PGRN_374_HOM-200925-125104', 'PGRN_375_HOM-200925-133107', 'PGRN_376_WT-200923-170108', 'PGRN_377_HET-200924-161222', 'PGRN_378_WT-200924-153304', 'PGRN_379_HET-200925-145125', 'PGRN_380_WT-200924-165826', 'PGRN_381_HET-200925-160959', 'PGRN_382_WT-200924-174422', 'PGRN_384_WT-200923-162127', 'PGRN_388_HOM-200925-153040', 'PGRN_389_WT-200923-150634', 'PGRN_391_HOM-200923-154439', 'PGRN_393_HOM-200925-183840', 'PGRN_396_HOM-200925-180218', 'PGRN_400_HET-200924-135146', 'PGRN_401_HET-200924-145359']

if __name__=='__main__':
    base_path = r'C:\Users\bsriram\Desktop\Data\PGRN_Coh2'
    tanks = awake_cohort2_tanks5
    pre_tank,post_tank1,post_tank2 = get_tank_sequencing(tanks)
    print(pre_tank)
    print(post_tank1)
    print(post_tank2)
    # for tank in post_tank1:
        # tank_path = os.path.join(base_path,tank)
        # data,ratio = load_tdt_tank(tank_path)
        # pack(data,tank_path,ratio=ratio)
        # pack_timestamps_dual_whisker_and_cams(tank_path)
        # print('done for {0}'.format(tank))
        # print()
        # print()
        # analyze_mua_by_channel_multistim(tank_path,type='raster',show_plot=True,min_z=2,z_zoom=2,interval=[-1,3])