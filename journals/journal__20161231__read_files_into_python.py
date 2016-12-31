import ipdb
ipdb.set_trace()

import scipy.io as sio

datDir    = '/SPENCEdata/Research/database/RENU2/processed--from_Dave_Hysell/'

#Add filenames for E data to list
temp_list = ['VLF_data.mat','plasma_freq.mat','HF_data.mat','E-fieldcalcs.mat','20160325-R2-convection.mat']
temp_list = [datDir + s for s in temp_list]


vlfdat    = sio.loadmat(temp_list[0])
pfreqdat  = sio.loadmat(temp_list[1])
hfdat     = sio.loadmat(temp_list[2])
efdat     = sio.loadmat(temp_list[3])
convdat   = sio.loadmat(temp_list[4])

print("Loaded VLFdat, PFreqdat, HFdat, EFdat, convdat")

