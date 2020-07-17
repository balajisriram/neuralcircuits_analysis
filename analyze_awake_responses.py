from spike2_utils import analyze_mua_by_channel_multistim
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import pickle
import pandas as pd
import os
import numpy as np
from scipy.stats import ttest_ind as ttest

allData = pd.DataFrame()
base_path = 'C:\Data\PGRN'
tanks = ['PGRN_155_WT-200330-130817','PGRN_155_WT-200401-154248','PGRN_155_WT-200401-155139','PGRN_155_WT-200428-142415','PGRN_155_WT-200428-143352','PGRN_155_WT-200501-144905','PGRN_155_WT-200501-150247','PGRN_158_WT-200330-134401','PGRN_158_WT-200401-162227','PGRN_158_WT-200401-163127','PGRN_158_WT-200428-150301','PGRN_158_WT-200428-151153','PGRN_158_WT-200501-155657','PGRN_158_WT-200501-161040','PGRN_175_HET-200330-125152','PGRN_175_HET-200401-152118','PGRN_175_HET-200401-153033','PGRN_175_HET-200428-140532','PGRN_175_HET-200428-141408','PGRN_175_HET-200501-141903','PGRN_175_HET-200501-143232','PGRN_181_HET-200330-132613','PGRN_181_HET-200401-160230','PGRN_181_HET-200401-161132','PGRN_181_HET-200428-144511','PGRN_181_HET-200428-145342','PGRN_181_HET-200501-152447','PGRN_181_HET-200501-153924','PGRN_182_WT-200330-123844','PGRN_182_WT-200401-150003','PGRN_182_WT-200401-151015','PGRN_182_WT-200428-135320','PGRN_182_WT-200501-134811','PGRN_182_WT-200501-140216','PGRN_183_HOM-200401-164914','PGRN_183_HOM-200401-165806','PGRN_183_HOM-200428-133435','PGRN_183_HOM-200428-134317','PGRN_183_HOM-200501-131710','PGRN_183_HOM-200501-133051','PGRN_184_HOM-200330-121157','PGRN_184_HOM-200330-121906','PGRN_184_HOM-200401-144100','PGRN_184_HOM-200401-170826','PGRN_184_HOM-200401-171613','PGRN_184_HOM-200428-131207','PGRN_184_HOM-200428-132202','PGRN_184_HOM-200501-123637','PGRN_184_HOM-200501-125015','PGRN_HT_175_LMRS2-200303-144850','PGRN_WT_158_RS-200302-151344']
dates = ['Acclimation','March','March','April','April','KX','KX','Acclimation','March','March','April','April','KX','KX','Acclimation','March','March','April','April','KX','KX','Acclimation','March','March','April','April','KX','KX','Acclimation','March','March','April','April','KX','March','March','April','April','KX','KX','Acclimation','Acclimation','March','March','March','April','April','KX','KX','Prior','Prior']
subj_ID = ['PGRN_155','PGRN_155','PGRN_155','PGRN_155','PGRN_155','PGRN_155','PGRN_155','PGRN_158','PGRN_158','PGRN_158','PGRN_158','PGRN_158','PGRN_158','PGRN_158','PGRN_175','PGRN_175','PGRN_175','PGRN_175','PGRN_175','PGRN_175','PGRN_175','PGRN_181','PGRN_181','PGRN_181','PGRN_181','PGRN_181','PGRN_181','PGRN_181','PGRN_182','PGRN_182','PGRN_182','PGRN_182','PGRN_182','PGRN_182','PGRN_183','PGRN_183','PGRN_183','PGRN_183','PGRN_183','PGRN_183','PGRN_184','PGRN_184','PGRN_184','PGRN_184','PGRN_184','PGRN_184','PGRN_184','PGRN_184','PGRN_184','PGRN_175','PGRN_158']
genotype = ['WT','WT','WT','WT','WT','WT','WT','WT','WT','WT','WT','WT','WT','WT','HET','HET','HET','HET','HET','HET','HET','HET','HET','HET','HET','HET','HET','HET','WT','WT','WT','WT','WT','WT','HOM','HOM','HOM','HOM','HOM','HOM','HOM','HOM','HOM','HOM','HOM','HOM','HOM','HOM','HOM','HET','WT']

assert len(tanks)==len(subj_ID), 'needs to be the same length'
assert len(tanks)==len(genotype), 'needs to be the same length'

all_locs = [os.path.join(base_path,t) for t in tanks]

# with PdfPages('All_responses_awake_new.pdf') as pdf:
    # f,ax = plt.subplots(1,1,figsize=(20,20))
    # ax.text(0.5,0.5,'All Responses',horizontalalignment='center',verticalalignment='center',fontsize=48)
    # pdf.savefig(f)
    # plt.close()
    # for loc,subj,gt in zip(all_locs,subj_ID,genotype):
        # print(loc)
        # f,units = analyze_mua_by_channel_multistim(loc,show_plot=False,min_z=1,z_zoom=10)
        # curr_units = pd.DataFrame(units)
        # curr_units['subject'] = subj
        # curr_units['genotype'] = gt
        # curr_units['location'] = 'VPLM'
        # allData = allData.append(curr_units,ignore_index=True,sort=False)
        # pdf.savefig(f)
        # plt.close()
    # f,ax = plt.subplots(1,1,figsize=(20,20))
    
# allData.to_pickle('AllUnitData_MUA_multistim_awake.pickle')

# for tank,sub,gen in zip(tanks,subj_ID,genotype):
    # print('tank::{0} - subj::{1} - genotype::{2}'.format(tank,sub,gen))
    
    
t_155 = [
(303,np.nan,289,np.nan,1071,1101,1021,np.nan,994,995,894,780,735,826,np.nan,790),
(108,np.nan,np.nan,1577,1344,1254,1275,np.nan,1014,1167,1490,np.nan,1368,1097,1291,1105),
(249,608,575,795,764,780,772,np.nan,820,1001,777,np.nan,734,674,658,653),]

t_158 = [
(np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,675,np.nan,np.nan,909,875,734,781,788,793,868),
(np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,1078,np.nan,1043,np.nan,981,np.nan,850,908,903,946),
]

t_175 = [
(212,113,229,174,np.nan,np.nan,np.nan,956,859,914,967,882,881,1066,1191,1047),
(280,400,416,155,391,859,np.nan,1092,1074,868,967,979,770,938,1150,1037),
(390,837,545,346,684,864,854,829,778,800,906,911,863,883,1062,981)
]

t_181 = [
(np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,953,1136,1632,1412,1124,1208,np.nan,np.nan,),
(np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,281,320,108,743,993,570,944,880,np.nan,np.nan,),
(np.nan,np.nan,np.nan,np.nan,126,829,102,965,105,102,1697,np.nan,np.nan,1144,1596,1142),
]

t_183 = [
(np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,108,132,170,139,123,142,163,632,171,163,),
(np.nan,np.nan,np.nan,np.nan,245,200,np.nan,np.nan,np.nan,291,np.nan,np.nan,194,533,161,284),
]

t_184 = [
(719,867,757,785,741,821,860,np.nan,738,841,np.nan,np.nan,869,918,np.nan,np.nan,),
]

# PLOT RESPOSES
x = [np.nan]
y = [np.nan]
z = [np.nan]
# import pdb; pdb.set_trace()

# f,ax = plt.subplots(1,figsize=(10,10))
# plt.scatter(np.array(t_155[0]),np.array(t_155[1]),color=[1,0,0],marker='.',s=50);x.extend(t_155[0]);y.extend(t_155[1])
# plt.scatter(np.array(t_155[0]),np.array(t_155[2]),color=[0.75,0,0],marker='.',s=50);x.extend(t_155[0]);y.extend(t_155[2])
# plt.scatter(np.array(t_155[2]),np.array(t_155[1]),color=[0.5,0,0],marker='.',s=50);x.extend(t_155[2]);y.extend(t_155[1])

# plt.scatter(np.array(t_158[0]),np.array(t_158[1]),color=[0,1,0],marker='.',s=50);x.extend(t_158[0]);y.extend(t_158[1])

# plt.scatter(np.array(t_175[0]),np.array(t_175[1]),color=[0,0,1],marker='.',s=50);x.extend(t_175[0]);y.extend(t_175[1])
# plt.scatter(np.array(t_175[0]),np.array(t_175[2]),color=[0,0,0.75],marker='.',s=50);x.extend(t_175[0]);y.extend(t_175[2])
# plt.scatter(np.array(t_175[2]),np.array(t_175[1]),color=[0,0,0.5],marker='.',s=50);x.extend(t_175[2]);y.extend(t_175[1])

# plt.scatter(np.array(t_181[0]),np.array(t_181[1]),color=[1,0,1],marker='.',s=50);x.extend(t_181[0]);y.extend(t_181[1])
# plt.scatter(np.array(t_181[0]),np.array(t_181[2]),color=[0.75,0,0.75],marker='.',s=50);x.extend(t_181[0]);y.extend(t_181[2])
# plt.scatter(np.array(t_181[2]),np.array(t_181[1]),color=[0.5,0,0.5],marker='.',s=50);x.extend(t_181[2]);y.extend(t_181[1])

# plt.scatter(np.array(t_183[0]),np.array(t_183[1]),color=[1,1,0],marker='.',s=50);x.extend(t_183[0]);y.extend(t_183[1])
# plt.plot([0,1200],[0,1200],'k--')
# ax.set_xlabel('Day 1 timecourse (ms)')
# ax.set_ylabel('Day 2 timecourse (ms)')
# ax.set_title('Day-to-day variability')

# assert len(x)==len(y),'lengths are different'
# x = np.array(x)
# y = np.array(y)
# bad = np.bitwise_or(np.isnan(x),np.isnan(y))
# temp = np.corrcoef(x[~bad],y[~bad])
# print(temp)


# plt.show()

x.extend(t_155[1]);x.extend(t_158[0]);x = np.array(x)
y.extend(t_175[1]);y.extend(t_181[1]);y = np.array(y)
z.extend(t_183[1]);z.extend(t_184[0]);z = np.array(z)

print('x:',np.nanmean(x),':',np.nanstd(x))
print('y:',np.nanmean(y),':',np.nanstd(y))
print('z:',np.nanmean(z),':',np.nanstd(z))

st,p = ttest(x,y,nan_policy='omit')
print('x vs y:',p)
st,p = ttest(x,z,nan_policy='omit')
print('x vs z:',p)
st,p = ttest(z,y,nan_policy='omit')
print('y vs z:',p)



# 