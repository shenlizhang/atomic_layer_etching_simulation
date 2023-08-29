# python script to delete etched products
# Shenli, 05/02/18

import numpy as np
from scipy.spatial import cKDTree
import sys
import os


# filename = "final.xyz"  # can be edited
filename = sys.argv[1]
output = sys.argv[2]
data = np.genfromtxt(filename, skip_header=2)
number = data.shape[0]
# number = 50
# add an attribute to data to mark its cluster
# 0: keep; 1: to be determined; 2:visited; 3-N: cluster numbers need to be deleted
# initial the data to be determined
data = np.concatenate(
    (np.arange(number).reshape((number, 1)), data), axis=1)
data = np.concatenate((data, np.ones((number, 1), dtype=int)), axis=1)
# take periodic boundary into consideration, for x and y direction, limit=3 A
# lx=(low, high,dx)


def periodic_replicate(data, lx, ly, limit):
    newdata = []
    for item in data:
        if item[2] <= lx[0] + limit:
            newdata.append(
                [item[0], item[1], item[2] + lx[2], item[3], item[4], item[5]])
        if item[2] >= lx[1] - limit:
            newdata.append(
                [item[0], item[1], item[2] - lx[2], item[3], item[4], item[5]])
        if item[3] <= ly[0] + limit:
            newdata.append([item[0], item[1], item[2],
                            item[3] + ly[2], item[4], item[5]])
        if item[3] >= ly[1] - limit:
            newdata.append([item[0], item[1], item[2],
                            item[3] - ly[2], item[4], item[5]])
    newdata = np.asarray(newdata)
    newdata = np.concatenate((data, newdata), axis=0)
    return newdata

newdata = periodic_replicate(data, [0, 38.3155, 38.3155], [
                             0, 38.3155, 38.3155], 3)
# get nearest neighbor list for each atom


def nearest_neighbors(x):
    Sidata = newdata[np.logical_or.reduce([newdata[:, 1] == 1])]
    Cldata = newdata[np.logical_or.reduce([newdata[:, 1] == 2])]
    tree1 = cKDTree(Sidata[:, 2:5])
    tree2 = cKDTree(Cldata[:, 2:5])
    ddSi, iiSi = tree1.query(x, k=[2, 3, 4, 5], distance_upper_bound=3)
    ddCl, iiCl = tree2.query(x, k=[1, 2, 3, 4], distance_upper_bound=2.6)
    iitotal = []
    for item in iiSi:
        try:
            iitotal.append(Sidata[item][0])
        except:
            pass
    for item in iiCl:
        try:
            iitotal.append(Cldata[item][0])
        except:
            pass
    return np.array(iitotal)


def mark_atoms(atom, data):
    keep, visited = True, False
    neighbor_list = np.unique(nearest_neighbors(atom[2:5]))
    neighbor_list = neighbor_list.astype(int)
    for neighbor in neighbor_list:
        if data[neighbor][5] == 0:
            break
        if data[neighbor][5] >= 3:
            keep, visited = False, False
            break
        atype = data[neighbor][1]
        if data[neighbor][5] == 2:
            continue
        if atype == 2:
            keep = False
        if atype == 1:
            visited = True
    return keep, visited, neighbor_list


def keep_function(atom, data, neighbor_list):  # where keep==True
    data[int(atom[0])][5] = 0
    for neighbor in neighbor_list:
        data[neighbor][5] = 0
    for i in range(0, data.shape[0]):
        if data[i][5] == 2:
            data[i][5] = 0
    return data


# if keep == False and visited == False:
def delcluster_function(s, atom, data, neighbor_list):
    data[int(atom[0])][5] = s
    for neighbor in neighbor_list:
        data[neighbor][5] = s
    for i in range(0, data.shape[0]):
        if data[i][5] == 2:
            data[i][5] = s
    s += 1
    return s, data


def determine(atom, data, neighbor_list):
    # print "visited=True for atom", atom[0]
    data[int(atom[0])][5] = 2
    # print data
    nextlist = []
    for neighbor in neighbor_list:
        if data[neighbor][5] == 2:
            continue
        else:
            data[neighbor][5] = 2
            nextlist.append(data[neighbor])
    return nextlist[0], data


def main(data):
    s = 3
    l = data[np.logical_or.reduce([data[:, 1] == 1])]
    for i in range(0, l.shape[0]):
        atom = l[i]
        keep, visited, neighbor_list = mark_atoms(atom, data)
        if keep == True:
            data = keep_function(atom, data, neighbor_list)

        if keep == False and visited == False:
            s, data = delcluster_function(s, atom, data, neighbor_list)

        if keep == False and visited == True:
            atom, data = determine(atom, data, neighbor_list)

    for i in range(0, data.shape[0]):
        if data[i][1] == 2 and data[i][5] == 1:
            data[i][5] = s + 1
    return data

final = main(data)
#-------counting how many Si and Cl will be removed----------------#
nSi, nCl = 0, 0
for item in final:
    if item[5] >= 3:
        if item[1]==1:
            nSi+=1
        if item[1]==2:
            nCl+=1
print "%d Si atoms and %d Cl atoms will be deleted." %(nSi, nCl)

#output = 'final_cluster.lammpstrj'
number = data.shape[0]
try:
    os.remove(output)
except OSError:
    pass
f2 = open(output, 'a')
f2.write("ITEM: TIMESTEP \n")
f2.write("0 \n ITEM: NUMBER OF ATOMS \n %d \n" % number)
f2.write("ITEM: BOX BOUNDS pp pp pp \n")
f2.write("%e %e \n" % (0, 38.3155))
f2.write("%e %e \n" % (0, 38.3155))
f2.write("%e %e \n" % (0, 156.698082))
f2.write("ITEM: ATOMS id type x y z mol \n")
for item in final:
    f2.write("%d %d %.6f %.6f %.6f %d \n" %
             (item[0], item[1], item[2], item[3], item[4], item[5]))
f2.close()
