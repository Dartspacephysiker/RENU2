import ipdb
import numpy as np
# import scipy.io as sio
import hatch_utils as hu

datDir = '/SPENCEdata/Research/database/RENU2/'
infil = 'velMAG.mat'
outfil = datDir + 'velMAG.csv'
hdrStr = 'thistGPS,vel1,vel2,vel3'
fmtStr = '%3d,%12.5f,%12.5f,%12.5f'

# Add filenames for E data to list
# temp_list = [infil]
# temp_list = [datDir + s for s in temp_list]


# velmagdat = sio.loadmat(temp_list[0])

# print("Loaded velmagdat")

# print('velmagdat has : ', velmagdat.keys())

# l_velmagdat = [len(v) for v in velmagdat.values()]

# print '*****VELMAGDAT*****'
# hu.print_dict(velmagdat)

velmagdat = hu.read_mat_file(infil, datDir)

alldat = np.concatenate((velmagdat['thistGPS'], velmagdat['velMAG']), axis=1)

print("Saving to {}".format(outfil))
np.savetxt(outfil, alldat, fmt='%3d,%12.5f,%12.5f,%12.5f',
           delimiter=',', newline='\n', header=hdrStr)
