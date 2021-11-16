import numpy as np

data=np.loadtxt('grid.dat')
rows=data.shape[0]
data=data[np.argsort(data[:,0])]

with open('sorted_grid.dat', 'w') as f:
    for i in range(rows):
        f.write(str(data[i,0]) + '\t' + str(data[i,1]) + '\n')
