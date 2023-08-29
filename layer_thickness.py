# python script to characterize the layer thickness of amorphous region
# Shenli, 05/03/18

import numpy as np
import matplotlib.pyplot as plt
import sys

filename = raw_input("Please enter the sci file to transfer:")
data = np.genfromtxt(filename, skip_header=2)


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

zmin, zmax = float(sys.argv[1]), float(sys.argv[2])
#zmin, zmax = 54, 154
deltaz = 2.7
Sinumber, Clnumber = density_profile(zmin, zmax, deltaz, data)
plt.plot(Sinumber[:, 0], Sinumber[:, 1], label='Si')
plt.plot(Clnumber[:, 0], Clnumber[:, 1], label='Cl')
plt.xlabel('z distance ($\AA$)')
plt.ylabel('number density (#/$\AA$)')
plt.legend()
plt.savefig("final_density_profile.png")
