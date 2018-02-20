from mpl_toolkits.mplot3d import Axes3D

vec = np.array([0,0,0])

tr.center(saveBpos)

# tmat = np.identity(3)

# tmat[:,0] = tmat[:,0] * 1.5
# tmat[:,1] = tmat[:,1] * 1.1


# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.scatter(tmat.T[:,0],tmat.T[:,1],tmat.T[:,2])

# tmat = np.asfortranarray(tmat)
# tr.trans(tmat,np.pi/6,np.array([0,0,1]),np.zeros(3))
# tmat = np.asanyarray(tmat)

# ax.scatter(tmat.T[:,0],tmat.T[:,1],tmat.T[:,2],"r")

print("Initial tmat", tmat)

dist, stretch, rot_mat = tr.test_ds(saveBpos,saveBpos,frac,atoms,tmat,vec)

print("Stretching only:", rot_mat)
print("dist not in line:",dist, stretch)
print("dist not in line sure:", np.sum(np.sqrt(np.sum((saveBpos - rot_mat.dot(saveBpos))**2,axis=0))))


# rot_mat[1,:] = -rot_mat[1,:]

# Plotting the cell vectors of A and B
fig = plt.figure()
ax = fig.add_subplot(111)
ax.quiver(np.ones(3), np.ones(3), rot_mat[0,:2], rot_mat[1,:2])
ax.quiver(-np.ones(3), -np.ones(3), tmat[0,:2], tmat[1,:2])
ax.set_xlim([-5, 5])
ax.set_ylim([-5, 5])


sBpos = rot_mat.dot(saveBpos)

tr.center(saveBpos)

# Plotting the Apos and Bpos overlayed
fig = plt.figure()
ax = fig.add_subplot(111)
ax.scatter(saveBpos.T[:,0],saveBpos.T[:,1])
ax.scatter(sBpos.T[:,0],sBpos.T[:,1])
maxXAxis = np.max([Apos.max(), Bpos.max()]) + 1
ax.set_xlim([-maxXAxis, maxXAxis])
ax.set_ylim([-maxXAxis, maxXAxis])
ax.set_aspect('equal')

tBcell = np.identity(3)

tBcell[:2,:2] = Bcell[:2,:2]

Bcell = tBcell

# New stretching

dist, stretch, FofV = tr.test_ds(saveBpos,saveBpos,frac,atoms,Bcell,vec)

iBcell = np.linalg.inv(Bcell)

dist, stretch, FofTV = tr.test_ds(saveBpos,saveBpos,frac,atoms,tmat.dot(Bcell),vec)

stretchedBpos = FofV.dot(iBcell).T.dot(FofTV).dot(iBcell).dot(saveBpos)

stretchedBcell = FofV.dot(iBcell).T.dot(FofTV).dot(iBcell).dot(Bcell)

print("dist in line sure:", np.sum(np.sqrt(np.sum((saveBpos - stretchedBpos)**2,axis=0))))


# Plotting the Apos and Bpos overlayed
fig = plt.figure()
ax = fig.add_subplot(111)
ax.scatter(saveBpos.T[:,0],saveBpos.T[:,1])
ax.scatter(stretchedBpos.T[:,0],stretchedBpos.T[:,1])
maxXAxis = np.max([Apos.max(), Bpos.max()]) + 1
ax.set_xlim([-maxXAxis, maxXAxis])
ax.set_ylim([-maxXAxis, maxXAxis])
ax.set_aspect('equal')

# Plotting the cell vectors of A and B
fig = plt.figure()
ax = fig.add_subplot(111)
ax.quiver(np.ones(3), np.ones(3), Bcell[0,:2], Bcell[1,:2])
ax.quiver(-np.ones(3), -np.ones(3), stretchedBcell[0,:2], stretchedBcell[1,:2])
ax.set_xlim([-5, 5])
ax.set_ylim([-5, 5])


plt.show()
