import ipdb
ipdb.set_trace()

import scipy.io as sio

datDir    = '/SPENCEdata/Research/database/RENU2/'

#Add filenames for E data to list
temp_list = ['velMAG.mat']
temp_list = [datDir + s for s in temp_list]


velmagdat    = sio.loadmat(temp_list[0])

print("Loaded velmagdat")

print('velmagdat has : ',velmagdat.keys())

l_velmagdat = [len(v) for v in velmagdat.values()]

print '*****VELMAGDAT*****'
for key in velmagdat:
    #    print key ,' : ', type(velmagdat[key])
    print '{:15} : {:25}'.format(key , type(velmagdat[key]), len(velmagdat[key]))

    print ' '
