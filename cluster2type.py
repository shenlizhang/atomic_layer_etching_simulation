#python script to change cluster output to xyz format data
import numpy as np

filename="final_cluster.txt"

output="final_cluster.xyz"
f=open(output,'a')
f.write('6211 \n Atoms. Timestep:0 \n')
data=np.genfromtxt(filename, skip_header=2)

for item in data:
	if item[5]==0:
		atype==1
	if 
	f.write('%d %.6f %.6f %.6f \n' %(item[5]))