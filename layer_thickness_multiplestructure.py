# python script to characterize the layer thickness of amorphous region
# Shenli, 05/03/18

import numpy as np
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt
import sys
import os

number = int(sys.argv[1])


def density_profile(zmin, zmax, deltaz, data):
    size = int((zmax - zmin) / deltaz) + 1
    Sinumber, Clnumber = [], []
    for i in range(0, size):
        z = zmin + deltaz * i
        Sicount, Clcount = 0, 0
        for item in data:
            if z - float(deltaz / 2) <= item[3] < z + float(deltaz / 2):
                if item[0] == 1:
                    Sicount += 1
                if item[0] == 2:
                    Clcount += 1
        Sidensity = float(Sicount / deltaz)
        Cldensity = float(Clcount / deltaz)
        Sinumber.append([z - zmin, Sidensity])
        Clnumber.append([z - zmin, Cldensity])
    Sinumber = np.asarray(Sinumber)
    Clnumber = np.asarray(Clnumber)
    return Sinumber, Clnumber


def FWHM(data):
    x = data[:, 0]
    spline = UnivariateSpline(x, data[:, 1] - np.max(data[:, 1]) / 2, s=0)
    r1, r2 = spline.roots()
    return r2 - r1

zmin, zmax = float(sys.argv[2]), float(sys.argv[3])
deltaz = 2.7

outfile = sys.argv[4]
try:
    os.remove(outfile)
except OSError:
    pass
f = open(outfile, 'a')
f.write("#structure No. width(Angstrom) \n")
for i in range(1, number + 1):
    try:
        filename = '5_%d.3_%d.xyz' % (i, i)
    except IOError:
        continue
    data = np.genfromtxt(filename, skip_header=2)
# zmin, zmax = 54, 154
    Sinumber, Clnumber = density_profile(zmin, zmax, deltaz, data)
    width = FWHM(Clnumber)
    f.write('%d %.4f \n' % (i + 1, width))
f.close()
