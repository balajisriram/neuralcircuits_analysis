import pandas as pd
import numpy as np
import os

ls = os.listdir
isdir = os.path.isdir
isfile = os.path.isfile
opj = os.path.join

def get_tdt_subject(tank):
    splits = tank.split('_')
    return(splits[0]+'_'+splits[1])
    
def get_tdt_genotype(tank):
    splits = tank.split('_')
    return splits[2].split('-')[0]

def get_tdt_date(tank):
    splits = tank.split('_')
    return splits[2].split('-')[1]
    
def get_tdt_time(tank):
    splits = tank.split('_')
    return splits[2].split('-')[2]

def has_ipsi_stim(path):
    return op.exists(op.join(path,'timestamps_L.np'))
    
def get_repeat_number(tank):
    if '200708' in tank or '200713' in tank or '200706' in tank or '200715' in tank:
        return 'Repeat1'
    elif '200723' in tank or '200721' in tank:
        return 'Repeat2'
    elif '200820' in tank or '200818' in tank:
        return 'Repeat3'
    elif '200923' in tank or '200924' in tank or '200925' in tank:
        return 'Repeat4'
    else:
        return 'Unknown'

def get_n_stim(tank,base_path=r'\\camhpcisixcifs.biogen.com\dept\electrophysiology\invivo_mouse\PGRN\PGRN_AwakeCohort2\DataAnalyzed'):
    try:
        # timestamps
        timestamps = np.fromfile(opj(base_path,tank,'timestamps_R.np'))
        istimi = np.mean(np.diff(timestamps))
        
        # is it single stim or multi_stim
        if np.abs(istimi-14)<1 or np.abs(istimi-7)<1:
            n_stim=20
        elif np.abs(istimi-10)<1 or np.abs(istimi-5)<1:
            n_stim=1
        else:
            print('Interstim interval is {0} for tank {1}. Setting n_stim to None'.format(istimi,tank))
            n_stim=None
    except Exception as ex:
        print(ex)
        print('tank::',tank)
        n_stim=None
    return n_stim
    
def get_stim_lateralism(tank,base_path=r'\\camhpcisixcifs.biogen.com\dept\electrophysiology\invivo_mouse\PGRN\PGRN_AwakeCohort2\DataAnalyzed'):
    has_l = False
    has_r = False
    try:
        if isfile(opj(base_path,tank,'timestamps_R.np')): has_r=True
    except Exception as ex:
        print(ex)
        
    try:
        if isfile(opj(base_path,tank,'timestamps_L.np')): has_l=True
    except Exception as ex:
        print(ex)
        
    if has_l and has_r : return 'bilateral'
    elif has_l and not has_r: return 'lateral_left'
    elif has_r and not has_l: return 'lateral_right'
    else: return 'no_stim'
    
awake_cohort1_tanks = ['PGRN_155_WT-200330-130817','PGRN_155_WT-200401-154248','PGRN_155_WT-200401-155139','PGRN_155_WT-200428-142415','PGRN_155_WT-200428-143352','PGRN_155_WT-200501-144905','PGRN_155_WT-200501-150247','PGRN_158_WT-200330-134401','PGRN_158_WT-200401-162227','PGRN_158_WT-200401-163127','PGRN_158_WT-200428-150301','PGRN_158_WT-200428-151153','PGRN_158_WT-200501-155657','PGRN_158_WT-200501-161040','PGRN_175_HET-200330-125152','PGRN_175_HET-200401-152118','PGRN_175_HET-200401-153033','PGRN_175_HET-200428-140532','PGRN_175_HET-200428-141408','PGRN_175_HET-200501-141903','PGRN_175_HET-200501-143232','PGRN_181_HET-200330-132613','PGRN_181_HET-200401-160230','PGRN_181_HET-200401-161132','PGRN_181_HET-200428-144511','PGRN_181_HET-200428-145342','PGRN_181_HET-200501-152447','PGRN_181_HET-200501-153924','PGRN_182_WT-200330-123844','PGRN_182_WT-200401-150003','PGRN_182_WT-200401-151015','PGRN_182_WT-200428-135320','PGRN_182_WT-200501-134811','PGRN_182_WT-200501-140216','PGRN_183_HOM-200401-164914','PGRN_183_HOM-200401-165806','PGRN_183_HOM-200428-133435','PGRN_183_HOM-200428-134317','PGRN_183_HOM-200501-131710','PGRN_183_HOM-200501-133051','PGRN_184_HOM-200330-121157','PGRN_184_HOM-200330-121906','PGRN_184_HOM-200401-144100','PGRN_184_HOM-200401-170826','PGRN_184_HOM-200401-171613','PGRN_184_HOM-200428-131207','PGRN_184_HOM-200428-132202','PGRN_184_HOM-200501-123637','PGRN_184_HOM-200501-125015','PGRN_HT_175_LMRS2-200303-144850','PGRN_WT_158_RS-200302-151344']

awake_cohort2_repeat1 = ['PGRN_369_WT-200713-140327','PGRN_370_HOM-200708-131420','PGRN_371_WT-200713-134124','PGRN_372_HET-200708-165326','PGRN_374_HOM-200713-102856','PGRN_375_HOM-200713-105801','PGRN_376_WT-200708-124614','PGRN_377_HET-200708-155016','PGRN_378_WT-200708-150655','PGRN_379_HET-200713-111941','PGRN_380_WT-200708-161417','PGRN_381_HET-200713-130234','PGRN_382_WT-200708-163354','PGRN_383_HOM-200708-143833','PGRN_384_WT-200706-164212','PGRN_384_WT-200715-123104','PGRN_388_HOM-200713-123705','PGRN_389_WT-200706-161140','PGRN_389_WT-200715-114822','PGRN_390_HET-200706-154900','PGRN_390_HET-200715-112446','PGRN_391_HOM-200706-163047','PGRN_391_HOM-200715-120933','PGRN_393_HOM-200713-145953','PGRN_396_HOM-200713-142749','PGRN_400_HET-200708-133938','PGRN_401_HET-200708-140920',]

awake_cohort2_repeat2=['PGRN_370_HOM-200721-140639','PGRN_370_HOM-200721-141333','PGRN_376_WT-200721-134456','PGRN_376_WT-200721-135218','PGRN_377_HET-200721-143537','PGRN_377_HET-200721-144337','PGRN_378_WT-200721-145141','PGRN_378_WT-200721-145931','PGRN_384_WT-200721-132153','PGRN_384_WT-200721-132927','PGRN_389_WT-200721-122505','PGRN_389_WT-200721-123728','PGRN_390_HET-200721-124605','PGRN_390_HET-200721-125503','PGRN_391_HOM-200721-130354','PGRN_391_HOM-200721-131154','PGRN_400_HET-200721-142125','PGRN_400_HET-200721-142857','PGRN_401_HET-200721-150644','PGRN_401_HET-200721-151431',]

awake_cohort2_repeat2_2=['PGRN_369_WT-200723-132704','PGRN_369_WT-200723-133516','PGRN_371_WT-200723-131157','PGRN_371_WT-200723-131853','PGRN_372_HET-200723-110143','PGRN_372_HET-200723-111033','PGRN_374_HOM-200723-112530','PGRN_374_HOM-200723-113210','PGRN_375_HOM-200723-114142','PGRN_375_HOM-200723-114840','PGRN_379_HET-200723-115754','PGRN_379_HET-200723-120738','PGRN_380_WT-200723-102746','PGRN_380_WT-200723-103602','PGRN_381_HET-200723-125642','PGRN_381_HET-200723-130329','PGRN_382_WT-200723-104510','PGRN_382_WT-200723-105226','PGRN_383_HOM-200723-123322','PGRN_383_HOM-200723-124229','PGRN_388_HOM-200723-121521','PGRN_388_HOM-200723-122209','PGRN_393_HOM-200723-140955','PGRN_393_HOM-200723-141642','PGRN_396_HOM-200723-134711','PGRN_396_HOM-200723-135345',]

awake_cohort2_repeat3=['PGRN_369_WT-200820-114036','PGRN_369_WT-200820-114919','PGRN_370_HOM-200818-123855','PGRN_370_HOM-200818-124535','PGRN_371_WT-200820-112534','PGRN_371_WT-200820-113158','PGRN_372_HET-200818-142954','PGRN_372_HET-200818-143637','PGRN_374_HOM-200818-144524','PGRN_374_HOM-200818-145210','PGRN_375_HOM-200818-150015','PGRN_375_HOM-200818-150645','PGRN_376_WT-200818-122154','PGRN_376_WT-200818-122812','PGRN_377_HET-200818-133837','PGRN_377_HET-200818-134502','PGRN_378_WT-200818-132329','PGRN_378_WT-200818-133007','PGRN_379_HET-200818-151436','PGRN_379_HET-200818-152104','PGRN_380_WT-200818-135417','PGRN_380_WT-200818-140105','PGRN_381_HET-200820-111010','PGRN_381_HET-200820-111634','PGRN_382_WT-200818-140937','PGRN_382_WT-200818-141607','PGRN_384_WT-200818-114324','PGRN_384_WT-200818-114953','PGRN_388_HOM-200818-153425','PGRN_388_HOM-200818-154112','PGRN_389_WT-200818-110750','PGRN_389_WT-200818-111418','PGRN_390_HET-200818-105318','PGRN_390_HET-200818-110039','PGRN_391_HOM-200818-112909','PGRN_391_HOM-200818-113530','PGRN_393_HOM-200820-121424','PGRN_393_HOM-200820-122041','PGRN_396_HOM-200820-120015','PGRN_396_HOM-200820-120631','PGRN_400_HET-200818-125348','PGRN_400_HET-200818-130007','PGRN_401_HET-200818-130825','PGRN_401_HET-200818-131443',]

awake_cohort2_tanks5=['PGRN_369_WT-200925-170556','PGRN_369_WT-200925-171323','PGRN_369_WT-200925-172542','PGRN_371_WT-200925-162655','PGRN_371_WT-200925-163541','PGRN_371_WT-200925-164756','PGRN_372_HET-200925-114651','PGRN_372_HET-200925-115922','PGRN_372_HET-200925-121226','PGRN_374_HOM-200925-123008','PGRN_374_HOM-200925-123753','PGRN_374_HOM-200925-125104','PGRN_375_HOM-200925-131008','PGRN_375_HOM-200925-131704','PGRN_375_HOM-200925-133107','PGRN_376_WT-200923-164027','PGRN_376_WT-200923-164823','PGRN_376_WT-200923-170108','PGRN_377_HET-200924-155214','PGRN_377_HET-200924-155954','PGRN_377_HET-200924-161222','PGRN_378_WT-200924-151011','PGRN_378_WT-200924-152045','PGRN_378_WT-200924-153304','PGRN_379_HET-200925-143015','PGRN_379_HET-200925-143852','PGRN_379_HET-200925-145125','PGRN_380_WT-200924-162939','PGRN_380_WT-200924-163722','PGRN_380_WT-200924-165826','PGRN_381_HET-200925-154922','PGRN_381_HET-200925-155742','PGRN_381_HET-200925-160959','PGRN_382_WT-200924-172207','PGRN_382_WT-200924-173146','PGRN_382_WT-200924-174422','PGRN_384_WT-200923-160011','PGRN_384_WT-200923-160859','PGRN_384_WT-200923-162127','PGRN_388_HOM-200925-151039','PGRN_388_HOM-200925-151821','PGRN_388_HOM-200925-153040','PGRN_389_WT-200923-143845','PGRN_389_WT-200923-145259','PGRN_389_WT-200923-150634','PGRN_391_HOM-200923-152335','PGRN_391_HOM-200923-153141','PGRN_391_HOM-200923-154439','PGRN_393_HOM-200925-181822','PGRN_393_HOM-200925-182617','PGRN_393_HOM-200925-183840','PGRN_396_HOM-200925-174223','PGRN_396_HOM-200925-174956','PGRN_396_HOM-200925-180218','PGRN_400_HET-200924-132658','PGRN_400_HET-200924-133908','PGRN_400_HET-200924-135146','PGRN_401_HET-200924-143215','PGRN_401_HET-200924-144141','PGRN_401_HET-200924-145359',]

awake_cohort2_pre_drug = ['PGRN_369_WT-200925-170556', 'PGRN_371_WT-200925-162655', 'PGRN_372_HET-200925-114651', 'PGRN_374_HOM-200925-123008', 'PGRN_375_HOM-200925-131008', 'PGRN_376_WT-200923-164027', 'PGRN_377_HET-200924-155214', 'PGRN_378_WT-200924-151011', 'PGRN_379_HET-200925-143015', 'PGRN_380_WT-200924-162939', 'PGRN_381_HET-200925-154922', 'PGRN_382_WT-200924-172207', 'PGRN_384_WT-200923-160011', 'PGRN_388_HOM-200925-151039', 'PGRN_389_WT-200923-143845', 'PGRN_391_HOM-200923-152335', 'PGRN_393_HOM-200925-181822', 'PGRN_396_HOM-200925-174223', 'PGRN_400_HET-200924-132658', 'PGRN_401_HET-200924-143215']

awake_cohort2_post_drug = ['PGRN_369_WT-200925-171323', 'PGRN_371_WT-200925-163541', 'PGRN_372_HET-200925-115922', 'PGRN_374_HOM-200925-123753', 'PGRN_375_HOM-200925-131704', 'PGRN_376_WT-200923-164823', 'PGRN_377_HET-200924-155954', 'PGRN_378_WT-200924-152045', 'PGRN_379_HET-200925-143852', 'PGRN_380_WT-200924-163722', 'PGRN_381_HET-200925-155742', 'PGRN_382_WT-200924-173146', 'PGRN_384_WT-200923-160859', 'PGRN_388_HOM-200925-151821', 'PGRN_389_WT-200923-145259', 'PGRN_391_HOM-200923-153141', 'PGRN_393_HOM-200925-182617', 'PGRN_396_HOM-200925-174956', 'PGRN_400_HET-200924-133908', 'PGRN_401_HET-200924-144141']

awake_cohort2_post_drug2 = ['PGRN_369_WT-200925-172542', 'PGRN_371_WT-200925-164756', 'PGRN_372_HET-200925-121226', 'PGRN_374_HOM-200925-125104', 'PGRN_375_HOM-200925-133107', 'PGRN_376_WT-200923-170108', 'PGRN_377_HET-200924-161222', 'PGRN_378_WT-200924-153304', 'PGRN_379_HET-200925-145125', 'PGRN_380_WT-200924-165826', 'PGRN_381_HET-200925-160959', 'PGRN_382_WT-200924-174422', 'PGRN_384_WT-200923-162127', 'PGRN_388_HOM-200925-153040', 'PGRN_389_WT-200923-150634', 'PGRN_391_HOM-200923-154439', 'PGRN_393_HOM-200925-183840', 'PGRN_396_HOM-200925-180218', 'PGRN_400_HET-200924-135146', 'PGRN_401_HET-200924-145359']


tank_details = []

# for tank in awake_cohort1_tanks:
    # this_tank = {}
    # this_tank['tank_name'] = tank
    # this_tank['subject'] = get_tdt_subject(tank)
    # this_tank['genotype'] = get_tdt_genotype(tank)
    # this_tank['date'] = get_tdt_date(tank)
    # this_tank['time'] = get_tdt_time(tank)
    # this_tank['cohort'] = 'Cohort1'
    # this_tank['repeat'] = 'N/A'
    # this_tank['raw_data_path'] = r'\\camhpcisixcifs.biogen.com\dept\electrophysiology\invivo_mouse\PGRN\PGRN_AwakeCohort1'
    # this_tank['drug_details'] = 'N/A'
    # tank_details.append(this_tank)

for tank in awake_cohort2_repeat1:
    this_tank = {}
    this_tank['tank_name'] = tank
    this_tank['subject'] = get_tdt_subject(tank)
    this_tank['genotype'] = get_tdt_genotype(tank)
    this_tank['date'] = get_tdt_date(tank)
    this_tank['time'] = get_tdt_time(tank)
    this_tank['cohort'] = 'Cohort2'
    this_tank['repeat'] = get_repeat_number(tank)
    this_tank['raw_data_path'] = r'\\camhpcisixcifs.biogen.com\dept\electrophysiology\invivo_mouse\PGRN\PGRN_AwakeCohort2\DataAnalyzed'
    this_tank['n_stim'] = get_n_stim(tank,base_path=this_tank['raw_data_path'])
    this_tank['stim_lateralism'] = get_stim_lateralism(tank,base_path=this_tank['raw_data_path'])
    this_tank['drug_details'] = 'N/A'
    assert isdir(opj(this_tank['raw_data_path'],this_tank['tank_name'])),'Couldnt find the raw data dir {0}'.format(opj(this_tank['raw_data_path'],this_tank['tank_name']))
    tank_details.append(this_tank)

for tank in awake_cohort2_repeat2:
    this_tank = {}
    this_tank['tank_name'] = tank
    this_tank['subject'] = get_tdt_subject(tank)
    this_tank['genotype'] = get_tdt_genotype(tank)
    this_tank['date'] = get_tdt_date(tank)
    this_tank['time'] = get_tdt_time(tank)
    this_tank['cohort'] = 'Cohort2'
    this_tank['repeat'] = get_repeat_number(tank)
    this_tank['raw_data_path'] = r'\\camhpcisixcifs.biogen.com\dept\electrophysiology\invivo_mouse\PGRN\PGRN_AwakeCohort2\DataAnalyzed'
    this_tank['n_stim'] = get_n_stim(tank,base_path=this_tank['raw_data_path'])
    this_tank['stim_lateralism'] = get_stim_lateralism(tank,base_path=this_tank['raw_data_path'])
    this_tank['drug_details'] = 'N/A'
    assert isdir(opj(this_tank['raw_data_path'],this_tank['tank_name'])),'Couldnt find the raw data dir {0}'.format(opj(this_tank['raw_data_path'],this_tank['tank_name']))
    tank_details.append(this_tank)
    
for tank in awake_cohort2_repeat2_2:
    this_tank = {}
    this_tank['tank_name'] = tank
    this_tank['subject'] = get_tdt_subject(tank)
    this_tank['genotype'] = get_tdt_genotype(tank)
    this_tank['date'] = get_tdt_date(tank)
    this_tank['time'] = get_tdt_time(tank)
    this_tank['cohort'] = 'Cohort2'
    this_tank['repeat'] = get_repeat_number(tank)
    this_tank['raw_data_path'] = r'\\camhpcisixcifs.biogen.com\dept\electrophysiology\invivo_mouse\PGRN\PGRN_AwakeCohort2\DataAnalyzed'
    this_tank['n_stim'] = get_n_stim(tank,base_path=this_tank['raw_data_path'])
    this_tank['stim_lateralism'] = get_stim_lateralism(tank,base_path=this_tank['raw_data_path'])
    this_tank['drug_details'] = 'N/A'
    assert isdir(opj(this_tank['raw_data_path'],this_tank['tank_name'])),'Couldnt find the raw data dir {0}'.format(opj(this_tank['raw_data_path'],this_tank['tank_name']))
    tank_details.append(this_tank)
    
for tank in awake_cohort2_repeat3:
    this_tank = {}
    this_tank['tank_name'] = tank
    this_tank['subject'] = get_tdt_subject(tank)
    this_tank['genotype'] = get_tdt_genotype(tank)
    this_tank['date'] = get_tdt_date(tank)
    this_tank['time'] = get_tdt_time(tank)
    this_tank['cohort'] = 'Cohort2'
    this_tank['repeat'] = get_repeat_number(tank)
    this_tank['raw_data_path'] = r'\\camhpcisixcifs.biogen.com\dept\electrophysiology\invivo_mouse\PGRN\PGRN_AwakeCohort2\DataAnalyzed'
    this_tank['n_stim'] = get_n_stim(tank,base_path=this_tank['raw_data_path'])
    this_tank['stim_lateralism'] = get_stim_lateralism(tank,base_path=this_tank['raw_data_path'])
    this_tank['drug_details'] = 'N/A'
    assert isdir(opj(this_tank['raw_data_path'],this_tank['tank_name'])),'Couldnt find the raw data dir {0}'.format(opj(this_tank['raw_data_path'],this_tank['tank_name']))
    tank_details.append(this_tank)
    
for tank in awake_cohort2_pre_drug:
    this_tank = {}
    this_tank['tank_name'] = tank
    this_tank['subject'] = get_tdt_subject(tank)
    this_tank['genotype'] = get_tdt_genotype(tank)
    this_tank['date'] = get_tdt_date(tank)
    this_tank['time'] = get_tdt_time(tank)
    this_tank['cohort'] = 'Cohort2'
    this_tank['repeat'] = get_repeat_number(tank)
    this_tank['raw_data_path'] = r'C:\Users\bsriram\Desktop\Data\PGRN_Coh2'
    this_tank['n_stim'] = get_n_stim(tank,base_path=this_tank['raw_data_path'])
    this_tank['stim_lateralism'] = get_stim_lateralism(tank,base_path=this_tank['raw_data_path'])
    this_tank['drug_details'] = 'None'
    assert isdir(opj(this_tank['raw_data_path'],this_tank['tank_name'])),'Couldnt find the raw data dir {0}'.format(opj(this_tank['raw_data_path'],this_tank['tank_name']))
    tank_details.append(this_tank)
    
for tank in awake_cohort2_post_drug:
    this_tank = {}
    this_tank['tank_name'] = tank
    this_tank['subject'] = get_tdt_subject(tank)
    this_tank['genotype'] = get_tdt_genotype(tank)
    this_tank['date'] = get_tdt_date(tank)
    this_tank['time'] = get_tdt_time(tank)
    this_tank['cohort'] = 'Cohort2'
    this_tank['repeat'] = get_repeat_number(tank)
    this_tank['raw_data_path'] = r'C:\Users\bsriram\Desktop\Data\PGRN_Coh2'
    this_tank['n_stim'] = get_n_stim(tank,base_path=this_tank['raw_data_path'])
    this_tank['stim_lateralism'] = get_stim_lateralism(tank,base_path=this_tank['raw_data_path'])
    this_tank['drug_details'] = 'EarlyDiazepam'
    assert isdir(opj(this_tank['raw_data_path'],this_tank['tank_name'])),'Couldnt find the raw data dir {0}'.format(opj(this_tank['raw_data_path'],this_tank['tank_name']))
    tank_details.append(this_tank)
    
for tank in awake_cohort2_post_drug2:
    this_tank = {}
    this_tank['tank_name'] = tank
    this_tank['subject'] = get_tdt_subject(tank)
    this_tank['genotype'] = get_tdt_genotype(tank)
    this_tank['date'] = get_tdt_date(tank)
    this_tank['time'] = get_tdt_time(tank)
    this_tank['cohort'] = 'Cohort2'
    this_tank['repeat'] = get_repeat_number(tank)
    this_tank['raw_data_path'] = r'C:\Users\bsriram\Desktop\Data\PGRN_Coh2'
    this_tank['n_stim'] = get_n_stim(tank,base_path=this_tank['raw_data_path'])
    this_tank['stim_lateralism'] = get_stim_lateralism(tank,base_path=this_tank['raw_data_path'])
    this_tank['drug_details'] = 'LateDiazepam'
    assert isdir(opj(this_tank['raw_data_path'],this_tank['tank_name'])),'Couldnt find the raw data dir {0}'.format(opj(this_tank['raw_data_path'],this_tank['tank_name']))
    tank_details.append(this_tank)

tank_df = pd.DataFrame(tank_details)
breakpoint()