# 2017/03/03
# def journal__20170303__read_conversions_from_csv():
# import pandas as pd
import csv
import hatch_utils as hu
# import numpy as np
import re
import scipy.io as sio
# import sys

indir = '/SPENCEdata/Research/database/RENU2/csv/'
infil = 'RENU2_VEL_FAC_CAR.csv'
outfil = re.sub('\.csv$', '.mat', infil)

data = []
f = open(indir + infil)
reader = csv.reader(f)
headers = next(reader)
column = {}
for h in headers:
    column[h] = []

print("Headers: {!s}".format(headers))

for row in reader:
    for h, v in zip(headers, row):
        column[h].append(v)

hu.print_dict(column)

print("Saving {} rows of {} entries apiece ({} total) ...".format(
    len(column), len(column[headers[0]]),
    len(column) * len(column[headers[0]])))

print("output: {!s}".format(indir + outfil))
sio.savemat(indir + outfil, column)

# EX. 1: The following works OK, but the above should be an improvement

# with open(indir + infil) as f:
#     reader = csv.reader(f)
#     for row in reader:
#         rowData = [float(elem) for elem in row]
#         data.append(rowData)

# matrix = np.array(data)
# sio.savemat(indir + outfil, {'csvmatrix': matrix})

# my_data = np.genfromtxt(indir + infil, delimiter=',', names=True)
# my_data = np.loadtxt(indir + infil, delimiter=',', skiprows=1)
# sio.savemat(indir + outfil, mdict, appendmat=True, format='5',
#             long_field_names=False, do_compression=False, oned_as='row')
