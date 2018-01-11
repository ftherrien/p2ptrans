from p2ptrans import transform as tr
import numpy as np
import numpy.linalg as la
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from p2ptrans import tiling as t
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

# Setting the unit cells of A and B
Acell = np.array([[1,0,0],[0,1,0],[0,0,1]])
Bcell = np.array([[1,0,0],[0,1,0],[0,0,1]])

# # Plotting the cell vectors of A and B
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.quiver(np.zeros(3), np.zeros(3), np.zeros(3), Bcell[0,:], Bcell[1,:], Bcell[2,:])
# ax.quiver(np.zeros(3), np.zeros(3), np.zeros(3), Acell[0,:], Acell[1,:], Acell[2,:])
# ax.set_xlim([0, 10])
# ax.set_ylim([0, 10])
# ax.set_zlim([0, 10])

ASC = t.sphere(Acell,2)
BSC = t.sphere(Bcell,2)

# # Plot gamma points of each A cell
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.scatter(ASC[0,:], ASC[1,:], ASC[2,:])
# maxXAxis = np.max([ASC[0,:].max(),ASC[1,:].max(),ASC[2,:].max()])
# minXAxis = np.min([ASC[0,:].min(),ASC[1,:].min(),ASC[2,:].min()])
# ax.set_xlim([minXAxis, maxXAxis])
# ax.set_ylim([minXAxis, maxXAxis])
# ax.set_zlim([minXAxis, maxXAxis])

# # Plot gamma points of each B cell
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.scatter(BSC[0,:], BSC[1,:], BSC[2,:])
# maxXAxis = np.max([BSC[0,:].max(),BSC[1,:].max(),BSC[2,:].max()])
# minXAxis = np.min([BSC[0,:].min(),BSC[1,:].min(),BSC[2,:].min()])
# ax.set_xlim([minXAxis, maxXAxis])
# ax.set_ylim([minXAxis, maxXAxis])
# ax.set_zlim([minXAxis, maxXAxis])

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.set_xlim3d(-1,1)
# ax.set_ylim3d(-1,1)
# ax.set_zlim3d(-1,1)

# # Create a random Apos
# #Apos = np.random.random((3,n))
# Apos = np.zeros((3,n))
# Apos[0,:] = np.linspace(-0.5,0.5,n)


# # Transform Apos to get Bpos
# tetha = 2*np.pi*np.random.random()
# vec = np.random.random((3,1))-0.5
# u = np.random.random((3,1))-0.5
# u = u/la.norm(u)
# Bpos = np.zeros((3,n))
# Bpos[0,:] = np.linspace(-1,1,n)
# Bpos = np.asfortranarray(Bpos)

# tr.trans(Bpos,tetha,u,vec)

# #Bpos = Bpos + np.random.random((3,n))*0.05
# # <><><><><><><><><><><><><><><><><>

# Adds atoms to A and B
atom_Apos = np.array([[0,0,0],
                      [0.5,0.5,0.5]])
atom_Bpos = np.array([[0,0,0],
                      [0.5,0.5,0.5]])
atoms = np.array([1,1]) # One of each kind


Apos = np.array([[],[],[]])
Bpos = np.array([[],[],[]])
for i in range(np.shape(atom_Apos)[0]):
    Apos = np.concatenate((Apos, ASC + Acell.dot(atom_Apos[i:i+1,:].T).dot(np.ones((1,np.shape(ASC)[1])))), axis=1)
    Bpos = np.concatenate((Bpos, ASC + Bcell.dot(atom_Bpos[i:i+1,:].T).dot(np.ones((1,np.shape(BSC)[1])))), axis=1)

Apos = np.random.random((3,300))

tetha = 2*np.pi*np.random.random()
vec = np.random.random((3,1))-0.5
u = np.random.random((3,1))-0.5
u = u/la.norm(u)
Bpos = np.asfortranarray(np.array(Apos))
print("1",Bpos)
tr.trans(Bpos,tetha,u,vec)
print("2",Bpos)

# ax.scatter(Apos.T[:,0],Apos.T[:,1],Apos.T[:,2])
# ax.scatter(Bpos.T[:,0],Bpos.T[:,1],Bpos.T[:,2])

Apos = np.asfortranarray(Apos)
print(Apos)
print(Bpos)
mapMat, dmin = tr.fastmapping(Apos, Bpos,np.array([1]),10000, 0.1, 0.1, 0.5)
#mapMat = mapMat-1 # Fortran index to python index

fig = plt.figure()

ax = fig.add_subplot(111, projection='3d')
ax.set_xlim3d(-1,1)
ax.set_ylim3d(-1,1)
ax.set_zlim3d(-1,1)
ax.scatter(Apos.T[:,0],Apos.T[:,1],Apos.T[:,2])
ax.scatter(Bpos.T[:,0],Bpos.T[:,1],Bpos.T[:,2])

Bpos = np.asanyarray(Bpos)
Apos = np.asanyarray(Apos)

# Directional plot
dirA = np.array([1,0,0])

print(dmin)
print(mapMat)

plt.show()



