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

#Here are your keys for these guys
print('vlfdat has : ',vlfdat.keys())
print('pfreqdat has : ',pfreqdat.keys())
print('hfdat has : ',hfdat.keys())
print('efdat has : ',efdat.keys())
print('convdat has : ',convdat.keys())

l_vlfdat = [len(v) for v in vlfdat.values()]
l_pfreqdat = [len(v) for v in pfreqdat.values()]
l_hfdat = [len(v) for v in hfdat.values()]
l_efdat = [len(v) for v in efdat.values()]
l_convdat = [len(v) for v in convdat.values()]

print '*****VLFDAT*****'
for key in vlfdat:
    #    print key ,' : ', type(vlfdat[key])
    print '{:15} : {:25}'.format(key , type(vlfdat[key]), len(vlfdat[key]))
print ' '

print '*****PFREQDAT*****'
for key in pfreqdat:
    print '{:15} : {:25} ({:7} elements)'.format(key , type(pfreqdat[key]), len(pfreqdat[key]))
print ' '

print '*****HFDAT*****'
for key in hfdat:
    print '{:15} : {:25} ({:7} elements)'.format(key , type(hfdat[key]), len(hfdat[key]))
print ' '

print '*****EFDAT*****'
for key in efdat:
    print '{:15} : {:25} ({:7} elements)'.format(key , type(efdat[key]), len(efdat[key]))
print ' '

print '*****CONVDAT*****'
for key in convdat:
    print '{:15} : {:25} ({:7} elements)'.format(key , type(convdat[key]), len(convdat[key]))
print ' '


#type(vlfdat['freqs34'])
