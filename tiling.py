from p2ptrans import tiling as t
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt


Acell = np.array([[1,0,0],[0,1,0],[0,0,1]])
Bcell = np.array([[1,0.3,0],[0.3,1,0],[0.6,0.3,0.3]]).T*2

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.quiver(np.zeros(3), np.zeros(3), np.zeros(3), Bcell[0,:], Bcell[1,:], Bcell[2,:])
ax.set_xlim([0, 10])
ax.set_ylim([0, 10])
ax.set_zlim([0, 10])

ASC = t.parallelepiped(Acell,5,5,5)
BSC = t.parallelepiped(Bcell,5,5,5)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(ASC[0,:], ASC[1,:], ASC[2,:])
maxXAxis = np.max([ASC[0,:].max(),ASC[1,:].max(),ASC[2,:].max()])
minXAxis = np.min([ASC[0,:].min(),ASC[1,:].min(),ASC[2,:].min()])
ax.set_xlim([minXAxis, maxXAxis])
ax.set_ylim([minXAxis, maxXAxis])
ax.set_zlim([minXAxis, maxXAxis])

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(BSC[0,:], BSC[1,:], BSC[2,:])
maxXAxis = np.max([BSC[0,:].max(),BSC[1,:].max(),BSC[2,:].max()])
minXAxis = np.min([BSC[0,:].min(),BSC[1,:].min(),BSC[2,:].min()])
ax.set_xlim([minXAxis, maxXAxis])
ax.set_ylim([minXAxis, maxXAxis])
ax.set_zlim([minXAxis, maxXAxis])


plt.show()
