import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import makepdb

def writeColumn(data, buffer_lines=2):
    col_width = max(len(word) for row in data for word in row) + buffer_lines
    out_list = []
    for row in data:
        out_list.append("".join(word.ljust(col_width) for word in row))
    return out_list

n = 720
 
golden_angle = np.pi * (3 - np.sqrt(5))
theta = golden_angle * np.arange(n)
z = np.linspace(1 - 1.0 / n, 1.0 / n - 1, n)
radius = np.sqrt(1 - z * z)
 
points = np.zeros((n, 3))
points[:,0] = radius * np.cos(theta)
points[:,1] = radius * np.sin(theta)
points[:,2] = z

a = points

b = np.asarray(a)
append_list = []
for entry in b:
    if entry[0] >= -0:
        append_list.append(entry.tolist())



new_list = []
for entry in append_list:
    temp = entry[:]
    temp[0] = -temp[0]
    new_list.append(temp)

np.asarray(new_list)
append_list = np.asarray(append_list)
append_list = np.concatenate((append_list,new_list), axis=0)




 
n = 120
 
radius = np.sqrt(np.arange(n) / float(n))
 
golden_angle = np.pi * (3 - np.sqrt(5))
theta = golden_angle * np.arange(n)
 
npoints = np.zeros((n, 3))
npoints[:,1] = np.cos(theta)
npoints[:,2] = np.sin(theta)
npoints *= radius.reshape((n, 1))
temp_list = []
for entry in npoints:
    if entry[1] > 0 and entry[2] > 0:
        temp_list.append(entry.tolist())

new_list = []        
for entry in temp_list:
    temp = entry[:]
    temp[1] = -temp[1]
    new_list.append(temp)
    temp = entry[:]
    temp[2] = -temp[2]
    new_list.append(temp)
    temp = entry[:]
    temp[1] = -temp[1]
    temp[2] = -temp[2]
    new_list.append(temp)


new_list.append([0,0,0])
temp_list = np.asarray(temp_list)
temp_list = np.concatenate((temp_list, new_list), axis=0)




append_list = np.concatenate((append_list, temp_list), axis=0)
#append_list = temp_list

append_list = append_list + 3
append_list = append_list*20
col_data = append_list.tolist()
for i, entry in enumerate(col_data):
    for j, number in enumerate(entry):
        print number
        col_data[i][j] = "{0:.3f}".format(number)
print col_data
output = writeColumn(col_data)
print output

with open('temp.txt', 'w') as f:
    for line in output:
        f.write('N  '+line+'\n')
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(append_list[:,0], append_list[:,1], append_list[:,2])
#ax.scatter(points[:,0], points[:,1], points[:,2])

ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
plt.show()

pdb = makepdb.PdbMaker(numAtoms=len(col_data))
#pdb.placement['x'][0] = [str(x) for x in append_list[:,0]]
pdb.placement['x'][0] = [x[0] for x in col_data]

pdb.placement['y'][0] = [x[1] for x in col_data]
pdb.placement['z'][0] = [x[2] for x in col_data]
print pdb.placement['y'][0]
pdb.calc_spacing()
print pdb.numAtoms
print 'spacing', pdb.spacing
pdb.makePDB()




