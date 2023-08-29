# python script to control Ar deposition and removel cycle
# Shenli,created 03/13/2018

import os
import sys
import subprocess
import re
import time
import shutil
import etched_product_cyclic
import numpy as np

#----------lammps-python loop-----------------#
output1 = "cluster_deletion.txt"
output2 = "cluster_type.txt"

for i in range(1, 200):  # define loop cycle number n here

    ###----------------------------- submit Ar deposition, removal and equilibration LAMMPs file-------------###
    '''in Ar_deposit lammps input file, it reads previously generated configuration,
       named "NVT_equilibrated_${i}.xyz", for example,
       and write out a configuration named "Ar_deposited_${i}.xyz" in the end.
       "i" is the variable name in lammps script.
       we submit the lammps script using batch file run.sh
       in the batch file, it contains a submission command:
       lmp -var i number -in lammps.input
       "number" should be changed to the value wanted'''

    # here it creates a new batch script for each cycle i'''
    os.system('cp run.sh run_%d.sh' % i)

    '''change the lammps variable "i" value number here using the follwing line.
    Notice here the read in configuration number is i-1 instead of i, as it is
    from previous cycle.
    '''
    os.system('sed -i "14s/\${i}/%d/g" run_%d.sh' % (i - 1, i))
    os.system('sed -i "14s/\${i}/%d/g" run_%d.sh' % (i - 1, i))
    # 14 is the line number in run_Ar_deposit_%d.sh, should be changed

    '''the following two lines are to get the submitted LAMMPS job id
    sbatch -N 1 -n 32 run.sh is the command to submit the simulation jobs'''
    submit_info = os.popen("sbatch -N 1 -n 32 run.sh").readline()
    job_id = submit_info.split()[-1]

    '''the following line checks if the current simulation has finished every 3600 second.
    if the job is still running, this procedure will simply wait and the following command
    lines won't be exected.'''
    while job_id in os.popen("squeue | grep " + job_id).readline():
        time.sleep(60)  # 3600 second could be changed


###---------------------use python script to delete etched products---------------------###
    etch1 = etched_product_cyclic.basic_functions()
    filename = "after_bombard_treated_%d.custom" % (i - 1)
    data, lx, ly, lz = etch1.convert2xyz(filename, 1)
    number = data.shape[0]
    data = np.concatenate(
        (np.arange(number).reshape((number, 1)), data), axis=1)
    data = np.concatenate((data, np.ones((number, 1), dtype=int)), axis=1)
    newdata = etch1.periodic_replicate(data, [0, lx, lx], [
        0, ly, ly], 3)
    final, s = etch1.main(data, newdata)
    outfile = "after_deposit_800_treated_%d.custom" % i
    etch1.write_configuration(outfile, lx, ly, lz, final)
    etch1.write_results(i, outpu1, output2, final, s)
