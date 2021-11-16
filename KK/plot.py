import numpy as np
import scipy as sp
import matplotlib
matplotlib.use('Agg')
import sys
import matplotlib.pyplot as plt

data = np.loadtxt(sys.argv[1])

en1 = data[:,0]
imeps1 = data[:,1]
reeps1 = data[:,2]
#en2 = data[:,3]
#imeps2 = data[:,4]
#reeps2 = data[:,5]

fig,ax=plt.subplots()


ax.plot(en1,imeps1,color='red', linestyle='solid', linewidth = 1.5,label='imeps1')
ax.plot(en1,reeps1,color='blue', linestyle='solid', linewidth = 1.5,label='reeps1')
#ax.plot(en2,imeps2,color='green', linestyle='--', linewidth = 1.5,label='imeps2')
#ax.plot(en2,reeps2,color='black', linestyle='--', linewidth = 1.5,label='reeps2')


#plt.show()
ax.legend()
plt.show()
plt.savefig('plot.pdf')
