from p2ptrans import transform as tr
import numpy as np
import numpy.linalg as la
# from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from p2ptrans import tiling as t
import time

random = False

# Setting the unit cells of A and B
Acell = np.array([[1,0],[0,1]]).T
#Bcell = np.array([[-1/2*1.1,1/2*1.1],[1,1]]).T
#Bcell = np.array([[-1/2,1/2],[1,1]]).T 
Bcell = np.array([[0.7,0],[0.3,0.5]]).T 


# Plotting the cell vectors of A and B
fig = plt.figure()
ax = fig.add_subplot(111)
ax.quiver(np.ones(3), np.ones(3), Bcell[0,:], Bcell[1,:])
ax.quiver(-np.ones(3), -np.ones(3), Acell[0,:], Acell[1,:])
ax.set_xlim([-5, 5])
ax.set_ylim([-5, 5])

ASC = t.circle(Acell,250)
BSC = t.circle(Bcell,250)

# Plot gamma points of each A cell
fig = plt.figure()
ax = fig.add_subplot(111)
ax.scatter(ASC[0,:], ASC[1,:])
maxXAxis = np.max([ASC[0,:].max(),ASC[1,:].max()])
minXAxis = np.min([ASC[0,:].min(),ASC[1,:].min()])
ax.set_xlim([minXAxis-1, maxXAxis+1])
ax.set_ylim([minXAxis-1, maxXAxis+1])
ax.set_aspect('equal')

# Plot gamma points of each B cell
fig = plt.figure()
ax = fig.add_subplot(111)
ax.scatter(BSC[0,:], BSC[1,:])
maxXAxis = np.max([BSC[0,:].max(),BSC[1,:].max()])
minXAxis = np.min([BSC[0,:].min(),BSC[1,:].min()])
ax.set_xlim([minXAxis-1, maxXAxis+1])
ax.set_ylim([minXAxis-1, maxXAxis+1])
ax.set_aspect('equal')


if random:
    # Create a random Apos and B with random small displacement
    atoms = np.array([1]) # One atom
    n = 30
    Apos = np.concatenate([np.random.random((2,n))*3, np.zeros((1,n))]) 
    
    # Transform Apos to get Bpos
    tetha = 2*np.pi*np.random.random()
    vec = np.random.random((3,1))-0.5
    vec[2] = 0
    u = np.array([0,0,1])
    Bpos = np.asfortranarray(np.array(Apos))

    tr.trans(Bpos,tetha,u,vec)

    randDisp = np.concatenate([np.random.random((2,n)), np.zeros((1,np.shape(Apos)[1]))])
    
    Bpos = Bpos + randDisp*0

    # Bpos = Bpos[:,:int(n/2)]
    
else:
    # Adds atoms to A and B
    atom_Apos = np.array([[0,0]])
    atom_Bpos = np.array([[0,0]])
    atoms = np.array([1]) # One atom

    Apos = np.array([[],[]])
    Bpos = np.array([[],[]])
    for i in range(np.shape(atom_Apos)[0]):
        Apos = np.concatenate((Apos, ASC + Acell.dot(atom_Apos[i:i+1,:].T).dot(np.ones((1,np.shape(ASC)[1])))), axis=1)
    for i in range(np.shape(atom_Bpos)[0]):
        Bpos = np.concatenate((Bpos, BSC + Bcell.dot(atom_Bpos[i:i+1,:].T).dot(np.ones((1,np.shape(BSC)[1])))), axis=1)

        
    Apos = np.concatenate([Apos,np.zeros((1,np.shape(Apos)[1]))]) # Adds the third dimension«
    Bpos = np.concatenate([Bpos,np.zeros((1,np.shape(Bpos)[1]))])
    
Apos = np.asfortranarray(Apos)
Bpos = np.asfortranarray(Bpos) 

t_time = time.time()
# mapMat, dmin = tr.fastmapping(Apos, Bpos, atoms,10000, 0.01, 0.00001, 5) # For dist 2 (not necessarily optimal)
frac = 0.5
nb = int(np.shape(Bpos)[1]*frac)
Acell_tmp = np.identity(3)
Acell_tmp[:2,:2] = Acell
saveBpos = np.array(Bpos)
Apos_map, Bpos, Bposst, n_map, tmat, dmin = tr.fastoptimization(Apos, Bpos, frac, Acell_tmp, la.inv(Acell_tmp), atoms, 250, 50, 3, 5, 5e-6, 5e-6) # For dist 6  
t_time = time.time() - t_time
Bpos = np.asanyarray(Bpos)
Apos = np.asanyarray(Apos)

Apos = Apos[:2,:]
Bpos = Bpos[:,:n_map]
Bposst = Bposst[:,:n_map]
Apos_map = Apos_map[:,:n_map]

#mapMat = mapMat-1 # Fortran index to python index)

# Plotting the Apos and Bpos overlayed
fig = plt.figure()
ax = fig.add_subplot(111)
#ax.scatter(Apos.T[:,0],Apos.T[:,1])
ax.scatter(Apos[:,:int(frac*125)].T[:,0],Apos[:,:int(frac*125)].T[:,1], c="C0")
ax.scatter(Apos[:,int(frac*125):].T[:,0],Apos[:,int(frac*125):].T[:,1], alpha=0.2, c="C0")
ax.scatter(Bpos.T[:,0],Bpos.T[:,1], alpha=0.5, c="C1")
maxXAxis = np.max([Apos.max(), Bpos.max()]) + 1
ax.set_xlim([-maxXAxis, maxXAxis])
ax.set_ylim([-maxXAxis, maxXAxis])
ax.set_aspect('equal')

# Mapping lattice
#Apos_map = Apos[:,mapMat[:np.shape(Bpos)[1]]]
disps = Apos_map - Bpos

#fig = plt.figure()
#ax = fig.add_subplot(111)
ax.quiver(Bpos.T[:,0], Bpos.T[:,1], disps.T[:,0], disps.T[:,1], scale_units='xy', scale=1)
maxXAxis = np.max([Apos.max(), Bpos.max()]) + 1
ax.set_xlim([-maxXAxis, maxXAxis])
ax.set_ylim([-maxXAxis, maxXAxis])
ax.set_aspect('equal')

# Plotting the Apos and Bposst overlayed
fig = plt.figure()
ax = fig.add_subplot(111)
ax.scatter(Apos[:,:int(frac*125)].T[:,0],Apos[:,:int(frac*125)].T[:,1], c="C0")
ax.scatter(Apos[:,int(frac*125):].T[:,0],Apos[:,int(frac*125):].T[:,1], alpha=0.2, c="C0")
ax.scatter(Bposst.T[:,0],Bposst.T[:,1], alpha=0.5, c="C1")
maxXAxis = np.max([Apos.max(), Bposst.max()]) + 1
ax.set_xlim([-maxXAxis, maxXAxis])
ax.set_ylim([-maxXAxis, maxXAxis])
ax.set_aspect('equal')

# Mapping lattice
#Apos_map = Apos[:,mapMat[:np.shape(Bposst)[1]]]
disps = Apos_map - Bposst

#fig = plt.figure()
#ax = fig.add_subplot(111)
ax.quiver(Bposst.T[:,0], Bposst.T[:,1],disps.T[:,0], disps.T[:,1], scale_units='xy', scale=1)
maxXAxis = np.max([Apos.max(), Bposst.max()]) + 1
ax.set_xlim([-maxXAxis, maxXAxis])
ax.set_ylim([-maxXAxis, maxXAxis])
ax.set_aspect('equal')


print(dmin)
print(sum(la.norm(disps,axis=0)))
#print(mapMat)
#print("Expected Order:", all(mapMat == np.arange(len(mapMat))))
print("Mapping time:", t_time)

plt.show()
