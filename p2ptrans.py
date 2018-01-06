from p2ptrans import transform as tr
import numpy as np
import numpy.linalg as la
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlim3d(-1,1)
ax.set_ylim3d(-1,1)
ax.set_zlim3d(-1,1)

n = 300

# Create a random Apos
#Apos = np.random.random((3,n))
Apos = np.zeros((3,n))
Apos[0,:] = np.linspace(-0.5,0.5,n)


# Transform Apos to get Bpos
tetha = 2*np.pi*np.random.random()
vec = np.random.random((3,1))-0.5
u = np.random.random((3,1))-0.5
u = u/la.norm(u)
Bpos = np.zeros((3,n))
Bpos[0,:] = np.linspace(-1,1,n)
Bpos = np.asfortranarray(Bpos)

tr.trans(Bpos,tetha,u,vec)

#Bpos = Bpos + np.random.random((3,n))*0.05
# <><><><><><><><><><><><><><><><><>

ax.scatter(Apos.T[:,0],Apos.T[:,1],Apos.T[:,2])
ax.scatter(Bpos.T[:,0],Bpos.T[:,1],Bpos.T[:,2])

Bpos = np.asfortranarray(Bpos)
Apos = np.asfortranarray(Apos)
cBpos = np.asfortranarray(Bpos)
cApos = np.asfortranarray(Apos)
cmapMat, cdmin = tr.fastmapping(cApos, cBpos,50, 0.1, 0.01)
mapMat, dmin = tr.fastmapping(Apos, Bpos, 50, 0.1, 0.01)
print("max diff",max(sum((cApos - Apos)**2,0)))
print("max diff",max(sum((cBpos - Bpos)**2,0)))
print(cdmin)
print(cmapMat)
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

print(dmin)
print(mapMat)

plt.show()



