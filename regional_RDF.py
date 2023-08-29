# python script to calculate RDF in different regions
# Shenli, 06/07/18
#-----running command-------#
# python regional_RDF.py filename zlow zhigh outname
# e.g. python regional_RDF.py final.sci 55 100 final

import numpy as np
import sys
from math import sqrt, pi
import matplotlib.pyplot as plt


def rdf_func(data, type1, type2, volume):
    deltar = 0.2
    rdf = []
    data1 = data[np.logical_or.reduce([data[:, 0] == type1])]
    data2 = data[np.logical_or.reduce([data[:, 0] == type2])]
    N1, N2 = data1.shape[0], data2.shape[0]
    for i in range(1, 50):
        r = deltar * i
        n = 0
        for item in data1:
            x0, y0, z0 = item[1], item[2], item[3]
            for item2 in data2:
                x1, y1, z1 = item2[1], item2[2], item2[3]
                if x1 == x0 and y1 == y0 and z1 == z0:
                    continue
                else:
                    distance = sqrt((x1 - x0)**2 + (y1 - y0)
                                    ** 2 + (z1 - z0)**2)
                    if distance >= r - deltar and distance < r:
                        n += 1
        value = n / (N1 * 4 * pi * r**2 * deltar * float(N2 / volume))
        rdf.append([r, value])
    return np.asarray(rdf)


def select_data(data, zlow, zhigh, s):
    newdata = []
    for item in data:
        if item[3] >= zlow and item[3] <= zhigh:
            newdata.append(item)
    volume = s * (zhigh - zlow)
    return np.asarray(newdata), volume


def convert2xyz(filename):
    data = []
    i = 0
    with open(filename) as f:
        for line in f:
            if not line in ['\n', '\r\n'] and len(line.split()) > 1:
                if line.split()[1] == '{}':
                    if line.split()[0] == '14':
                        index = 1
                    if line.split()[0] == '17':
                        index = 2
                    x = float(line.split()[3].strip('{'))
                    y = float(line.split()[4])
                    z = float(line.split()[5].strip('}'))
                    data.append([index, x, y, z])
            if "{a b c A B G} {0 0 0}" in line:
                lattice = line.split('} ')[2].strip('{')
                lx = float(lattice.split()[0])
                ly = float(lattice.split()[1])
                lz = float(lattice.split()[2])
    # print lx, ly, lz
    outdata = []
    for item in data:
        outdata.append([item[0], item[1] * lx, item[2] * ly, item[3] * lz])
    outdata = np.asarray(outdata)
    return outdata, lx * ly

#------main-------------#
filename = sys.argv[1]
zlow, zhigh = float(sys.argv[2]), float(sys.argv[3])
outname = sys.argv[4]
data, s = convert2xyz(filename)
selected, volume = select_data(data, zlow, zhigh, s)
rdf11 = rdf_func(selected, 1, 1, volume)
np.savetxt(outname + 'SiSi_rdf_%.0f_%.0f.txt' % (zlow, zhigh), rdf11)
#rdf12 = rdf_func(selected, 1, 2, volume)
#np.savetxt(outname + 'SiCl_rdf_%.0f_%.0f.txt' % (zlow, zhigh), rdf12)
plt.plot(rdf11[:, 0], rdf11[:, 1], label="SiSi")
#plt.plot(rdf12[:, 0], rdf12[:, 1], label="SiCl")
plt.legend(loc='best')
plt.xlabel('r diantce ($\AA$)')
plt.ylabel('g(r)')
plt.savefig(outname + '_%.0f_%.0f.png' % (zlow, zhigh))
