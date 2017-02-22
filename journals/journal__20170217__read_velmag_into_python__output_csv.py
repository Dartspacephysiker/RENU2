# import ipdb
# ipdb.set_trace()

import numpy as np
import scipy.io as sio
import hatch_utils as hu

datDir = '/SPENCEdata/Research/database/RENU2/'
outfil = datDir + 'velMAG.csv'
hdrStr = 'thistGPS,vel1,vel2,vel3'
fmtStr = '%3d,%12.5f,%12.5f,%12.5f'

# Add filenames for E data to list
temp_list = ['velMAG.mat']
temp_list = [datDir + s for s in temp_list]


velmagdat = sio.loadmat(temp_list[0])

print("Loaded velmagdat")

print('velmagdat has : ', velmagdat.keys())

l_velmagdat = [len(v) for v in velmagdat.values()]

print '*****VELMAGDAT*****'
hu.print_dict(velmagdat)

alldat = np.concatenate((velmagdat['thistGPS'], velmagdat['velMAG']), axis=1)
np.savetxt(outfil, alldat, fmt='%3d,%12.5f,%12.5f,%12.5f',
           delimiter=',', newline='\n', header=hdrStr)
