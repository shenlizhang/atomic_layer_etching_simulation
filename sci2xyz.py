# python script to transfer .sci file to .xyz file
# To run the script, open a command window and go to the directory
# containing both this script and sci files, and then type "python sci2xyz.py"
# the program will ask you to input the filename, just type final.sci for example (no "" is needed)


import os

filename = raw_input("Please enter the sci file to transfer:")
#filename = 'final.sci'
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
print lx, ly, lz
atom_num = len(data)

output = filename.split('.')[0] + ".xyz"
try:
    os.remove(output)
except OSError:
    pass
f = open(output, 'a')
f.write('%d \n' % atom_num)
f.write('Atoms. Timestep:0 \n')
for item in data:
    f.write('%d %f %f %f \n' %
            (item[0], item[1] * lx, item[2] * ly, item[3] * lz))
f.close()
