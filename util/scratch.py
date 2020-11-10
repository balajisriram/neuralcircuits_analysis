import os

ls = os.listdir
isdir = os.path.isdir
opj = os.path.join

base = r'\\camhpcisixcifs.biogen.com\dept\electrophysiology\invivo_mouse\PGRN\PGRN_AwakeCohort2\DataAnalyzed'
folders = [f for f in ls(base) if isdir(opj(base,f))]

for folder in folders:
    dat_files = [f for f in ls(opj(base,folder)) if 'dat' in f]
    if len(dat_files)==0:
        print(folder)