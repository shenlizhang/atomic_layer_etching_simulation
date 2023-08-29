# python script to delete etched products
# Shenli, 05/02/18
#-----running command-------#
# python delete_etched_products_custom.py

import numpy as np
from scipy.spatial import cKDTree
import sys
import os


class basic_functions(object):

    def __init__(self):
        return

    def convert2xyz(self, filename, frame):
        atom_num = int(np.genfromtxt(filename, skip_header=3, max_rows=1))
        data = np.genfromtxt(filename, skip_header=frame *
                             (atom_num + 9) + 9, max_rows=atom_num, usecols=(1, 2, 3, 4))
        alx = np.genfromtxt(filename, skip_header=frame *
                            (atom_num + 9) + 5, max_rows=1)
        aly = np.genfromtxt(filename, skip_header=frame *
                            (atom_num + 9) + 6, max_rows=1)
        alz = np.genfromtxt(filename, skip_header=frame *
                            (atom_num + 9) + 7, max_rows=1)
        lx = alx[1] - alx[0]
        ly = aly[1] - aly[0]
        lz = alz[1] - alz[0]
        # print lx, ly, lz
        return data, lx, ly, lz

    def periodic_replicate(self, data, lx, ly, limit):
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

    # get nearest neighbor list for each atom

    def nearest_neighbors(self, x, newdata):
        Gedata = newdata[np.logical_or.reduce([newdata[:, 1] == 1])]
        Cldata = newdata[np.logical_or.reduce([newdata[:, 1] == 2])]
        tree1 = cKDTree(Gedata[:, 2:5])
        tree2 = cKDTree(Cldata[:, 2:5])
        ddGe, iiGe = tree1.query(x, k=[2, 3, 4, 5], distance_upper_bound=2.85)
        ddCl, iiCl = tree2.query(x, k=[1, 2, 3, 4], distance_upper_bound=2.65)
        iitotal = []
        for item in iiGe:
            try:
                iitotal.append(Gedata[item][0])
            except:
                pass
        for item in iiCl:
            try:
                iitotal.append(Cldata[item][0])
            except:
                pass
        return np.array(iitotal)

    def mark_atoms(self, atom, data, newdata):
        keep, visited = True, False
        neighbor_list = np.unique(self.nearest_neighbors(atom[2:5], newdata))
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

    def keep_function(self, atom, data, neighbor_list):  # where keep==True
        data[int(atom[0])][5] = 0
        for neighbor in neighbor_list:
            data[neighbor][5] = 0
        for i in range(0, data.shape[0]):
            if data[i][5] == 2:
                data[i][5] = 0
        return data

    # if keep == False and visited == False:
    def delcluster_function(self, s, atom, data, neighbor_list):
        data[int(atom[0])][5] = s
        for neighbor in neighbor_list:
            data[neighbor][5] = s
        for i in range(0, data.shape[0]):
            if data[i][5] == 2:
                data[i][5] = s
        s += 1
        return s, data

    def determine(self, atom, data, neighbor_list):
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

    def clustering(self, data, n):
        GeCl, GeCl2, GeCl3, GeCl4, GenClx = 0, 0, 0, 0, 0
        for i in range(3, n + 1):
            cluster = data[np.logical_or.reduce([data[:, 5].astype(int) == i])]
            Ge, Cl = 0, 0
            for item in cluster:
                if int(item[1]) == 1:
                    Ge += 1
                if int(item[1]) == 2:
                    Cl += 1
            if Ge == 1 and Cl == 1:
                GeCl += 1
            if Ge == 1 and Cl == 2:
                GeCl2 += 1
            if Ge == 1 and Cl == 3:
                GeCl3 += 1
            if Ge == 1 and Cl == 4:
                GeCl4 += 1
            if Ge > 1:
                GenClx += 1
        a = np.array([GeCl, GeCl2, GeCl3, GeCl4, GenClx])
        total = np.sum(a)
        if total == 0:
            percent = np.zeros(5)
        else:
            percent = a / float(total)
        return percent, total

    def main(self, data, newdata):
        s = 3
        l = data[np.logical_or.reduce([data[:, 1] == 1])]
        for i in range(0, l.shape[0]):
            atom = l[i]
            keep, visited, neighbor_list = self.mark_atoms(atom, data, newdata)
            if keep == True:
                data = self.keep_function(atom, data, neighbor_list)

            if keep == False and visited == False:
                s, data = self.delcluster_function(
                    s, atom, data, neighbor_list)

            if keep == False and visited == True:
                atom, data = self.determine(atom, data, neighbor_list)

        for i in range(0, data.shape[0]):
            if data[i][1] == 2 and data[i][5] == 1:
                data[i][5] = s + 1
        return data, s + 1

    def write_configuration(self, output, lx, ly, lz, data):
        n = 0
        for item in data:
            if item[5] <= 2:
                n = n + 1
        try:
            os.remove(output)
        except OSError:
            pass
        f2 = open(output, 'a')
        f2.write("ITEM: TIMESTEP \n")
        f2.write("0 \nITEM: NUMBER OF ATOMS\n%d\n" % n)
        f2.write("ITEM: BOX BOUNDS pp pp pp \n")
        f2.write("%e %e \n" % (0, lx))
        f2.write("%e %e \n" % (0, ly))
        f2.write("%e %e \n" % (0, lz))
        f2.write("ITEM: ATOMS id type x y z\n")
        n = 0
        for item in data:
            if item[5] <= 2:
                f2.write("%d %d %.6f %.6f %.6f \n" %
                         (n, item[1], item[2], item[3], item[4]))
            n += 1
        f2.close()
        return

    def write_results(self, i, output1, output2, data, s):
        totGe, totCl = [], []
        totpercent = []
        nGe, nCl = 0, 0
        for item in data:
            if item[5] >= 3:
                if item[1] == 1:
                    nGe += 1
                if item[1] == 2:
                    nCl += 1

        totGe.append(nGe)
        totCl.append(nCl)
        f = open(output1, 'a')
        f.write(" %d %d %d \n" % (i, nGe, nCl))
        f.close()
        percent, totalC = self.clustering(data, s)
        totpercent.append(percent)
        f2 = open(output2, 'a')
        f2.write("%d %.3f %.3f %.3f %.3f %.3f %d \n" %
                 (i, percent[0], percent[1], percent[2], percent[3], percent[4], totalC))
        f2.close()
        return
