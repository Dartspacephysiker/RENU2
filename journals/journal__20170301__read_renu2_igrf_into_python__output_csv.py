import ipdb
import numpy as np
# import scipy.io as sio
import hatch_utils as hu

datdir = '/SPENCEdata/Research/database/RENU2/IGRF/'
fpref = 'RENU2_IGRF'
infil = fpref + '.mat'
outfil = fpref + '.csv'
outfil2 = fpref + '2.csv'
# hdrstr = 'thistGPS,vel1,vel2,vel3'
# fmtstr = '%3d,%12.5f,%12.5f,%12.5f'
fmtstrsing = '%.8E'
fmtstrsing2 = '%d'

igrf = hu.read_mat_file(infil, datdir)

sizes = []
arrkeys = []
for key in igrf:
    if type(igrf[key]) is np.ndarray:
        sizes.append(igrf[key].size)
        arrkeys.append(key)
pickit = max(sizes)
pickit2 = min(sizes)
arrkeys.sort()

arrkeys2 = arrkeys

hdrstr = ''
fmtstr = ''
ind = 0
keepkeys = []
while not hdrstr:
    if igrf[arrkeys[ind]].size == pickit:
        hdrstr = arrkeys.pop(ind)
        fmtstr += fmtstrsing
        keepkeys.append(hdrstr)
        break
    else:
        ind += 1

nkeepers = 1
for key in arrkeys:
    if igrf[key].size == pickit:
        hdrstr = hdrstr + ',{!s}'.format(key)
        fmtstr += ',' + fmtstrsing
        keepkeys.append(key)
        nkeepers += 1
print('Nkeepers ({}): {}'.format(nkeepers, hdrstr))

# Now the smaller stuff
hdrstr2 = ''
fmtstr2 = ''
ind = 0
keepkeys2 = []
exit = 0
while not hdrstr2:
    if igrf[arrkeys2[ind]].size == pickit2:
        try:
            float(igrf[arrkeys2[ind]].astype('int'))
            hdrstr2 = arrkeys2.pop(ind)
            fmtstr2 += fmtstrsing2
            keepkeys2.append(hdrstr2)
            print("keepkeys2:", keepkeys2[-1])
            break
        except ValueError:
            print("Not an int: {}".format(igrf[arrkeys2[ind]][0]))
            ind += 1
    else:
        ind += 1


nkeepers2 = 1
for key in arrkeys2:
    if igrf[key].size == pickit2:
        try:
            float(igrf[key].astype('int'))
            hdrstr2 = hdrstr2 + ',{!s}'.format(key)
            fmtstr2 += ',' + fmtstrsing2
            keepkeys2.append(key)
            nkeepers2 += 1
        except:
            print("Not an int: {}".format(igrf[key][0]))

print('Nkeepers2 ({}): {}'.format(nkeepers2, hdrstr2))


final = np.empty((nkeepers, pickit))
final2 = np.empty((nkeepers2, pickit2))

for i in enumerate(keepkeys):
    final[i[0], 0:] = (igrf[i[1]])[0:]

for i in enumerate(keepkeys2):
    print(i)
    print(igrf[i[1]], type(igrf[i[1]].astype('int')))
    final2[i[0], 0] = (igrf[i[1]])[0]

print("Saving to {}".format(outfil))

np.savetxt(datdir + outfil, final.transpose(), fmt=fmtstr,
           delimiter=',', newline='\n', header=hdrstr)

print("Saving to {}".format(outfil2))

np.savetxt(datdir + outfil2, final2.transpose(), fmt=fmtstr2,
           delimiter=',', newline='\n', header=hdrstr2)
