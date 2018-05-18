from p2ptrans import transform as tr
import numpy as np
import numpy.linalg as la
# from mpl_toolkits.mplot3d import Axes3D
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from p2ptrans import tiling as t
import pickle
import time
from pylada.crystal import Structure, primitive, gruber

tol = 1e-5

def find_cell(class_list, positions, tol = 1e-5, frac_tol = 0.5):
    cell = np.zeros((3,3))
    for i in range(len(np.unique(class_list))):
        newcell = np.identity(3)
        pos = positions[:, class_list == i]
        pos = pos[:,1:] - pos[:,0:1].dot(np.ones((1,np.shape(pos)[1]-1))) # centered
        norms = la.norm(pos, axis = 0)
        idx = np.argsort(norms)
        j = 0;
        for k in idx:
            multiple = [0]
            if j == 1:
                if la.norm(np.cross(pos[:,k], newcell[:,0])) < tol:
                    continue
            elif j == 2:
                if abs(la.det(np.concatenate((newcell[:,:2],pos[:,k:k+1]), axis=1))) < tol:
                    continue
            for p in pos.T:
                if la.norm(np.cross(p,pos[:,k])) < tol:
                    multiple.append(p.dot(pos[:,k])/la.norm(pos[:,k])**2)
            if multiple != []:
                multiple.sort()
                shell = np.sum(norms < abs(multiple[0])*la.norm(pos[:,k]))/float(len(norms))
                shell = min(shell, np.sum(norms < multiple[-1]*la.norm(pos[:,k]))/float(len(norms)))
                if len(multiple) == len(np.arange(round(multiple[0]), round(multiple[-1])\
    +1)):
                    if np.allclose(multiple, np.arange(round(multiple[0]), round(multiple[-1])+1), tol) and shell > frac_tol:
                        newcell[:,j] = pos[:,k]
                        j += 1
                        if j == 2: break #In 2D only TMP TODO
        else:
            raise RuntimeError("Could not find periodic cell for displacement %d"%i)
        if i==0:
            cell = newcell
        elif la.det(cell) < la.det(newcell):
            if np.allclose(la.inv(cell).dot(newcell), np.round(la.inv(cell).dot(newcell)), tol):
                cell = newcell
            else:
                raise RuntimeError("The periodicity of the different classes of displacement is different")
        else:
            if not np.allclose(la.inv(newcell).dot(cell), np.round(la.inv(newcell).dot(cell)),\
     tol):
                raise RuntimeError("The periodicity of the different classes of displacement is different")
    return cell

def classify(disps, tol = 1.e-1):
    vec_classes = [disps[:,0:1]]
    class_list = np.zeros(np.shape(disps)[1], np.int)
    for i in range(np.shape(disps)[1]):
        classified = False
        for j, vec_class in enumerate(vec_classes):
            vec_mean = np.mean(vec_class, axis=1)
            if (abs(la.norm(vec_mean) - vec_mean.T.dot(disps[:,i])/la.norm(vec_mean)) < tol and 
                la.norm(np.cross(vec_mean, disps[:,i]))/la.norm(vec_mean) < tol):
                vec_classes[j] = np.concatenate((vec_class,disps[:,i:i+1]),axis=1)
                if j == 3:
                    print(vec_mean, disps[:,i])
                class_list[i] = j
                classified = True
                break
        if not classified:
            vec_classes.append(disps[:,i:i+1])
            class_list[i] = len(vec_classes) - 1
            
    for i, elem in enumerate(vec_classes):
        vec_classes[i] = np.mean(elem, axis=1)
        
    return class_list, vec_classes

# \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
# Begin Program

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
fig.savefig('CellVectors.png')

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
fig.savefig('Agrid.png')

# Plot gamma points of each B cell
fig = plt.figure()
ax = fig.add_subplot(111)
ax.scatter(BSC[0,:], BSC[1,:])
maxXAxis = np.max([BSC[0,:].max(),BSC[1,:].max()])
minXAxis = np.min([BSC[0,:].min(),BSC[1,:].min()])
ax.set_xlim([minXAxis-1, maxXAxis+1])
ax.set_ylim([minXAxis-1, maxXAxis+1])
ax.set_aspect('equal')
fig.savefig('Bgrid.png')

# For testing purposes
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
    # Adds atoms to A and B (for cell with different types of atoms)
    atom_Apos = np.array([[0,0]])
    atom_Bpos = np.array([[0,0]])
    atoms = np.array([1]) # One atom

    Apos = np.array([[],[]])
    Bpos = np.array([[],[]])
    for i in range(np.shape(atom_Apos)[0]):
        Apos = np.concatenate((Apos, ASC + Acell.dot(atom_Apos[i:i+1,:].T).dot(np.ones((1,np.shape(ASC)[1])))), axis=1)
    for i in range(np.shape(atom_Bpos)[0]):
        Bpos = np.concatenate((Bpos, BSC + Bcell.dot(atom_Bpos[i:i+1,:].T).dot(np.ones((1,np.shape(BSC)[1])))), axis=1)

        
    Apos = np.concatenate([Apos,np.zeros((1,np.shape(Apos)[1]))]) # TMP Only in 2D Adds the 3rd dim. 
    Bpos = np.concatenate([Bpos,np.zeros((1,np.shape(Bpos)[1]))])


frac = 0.5
nb = int(np.shape(Bpos)[1]*frac)
Acell_tmp = np.identity(3)
Acell_tmp[:2,:2] = Acell

Apos = np.asfortranarray(Apos)
Bpos = np.asfortranarray(Bpos) 
t_time = time.time()
# Apos_map, Bpos, Bposst, n_map, tmat, dmin = tr.fastoptimization(Apos, Bpos, frac, Acell_tmp, la.inv(Acell_tmp), atoms, 250, 50, 3, 5, 5e-6, 5e-6) # For dist 6
t_time = time.time() - t_time
Bpos = np.asanyarray(Bpos)
Apos = np.asanyarray(Apos)

print(dmin)
print("Mapping time:", t_time)

# pickle.dump((Apos_map, Bpos, Bposst, n_map, tmat, dmin), open("fastoptimization.dat","wb"))

# TMP for testing only -->
tr.center(Apos)
tr.center(Bpos)
Apos_map, Bpos, Bposst, n_map, tmat, dmin = pickle.load(open("fastoptimization.dat","rb"))
# <--  

Apos = Apos[:2,:]
Bpos = Bpos[:,:n_map]
Bposst = Bposst[:,:n_map]
Apos_map = Apos_map[:,:n_map]

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

# Displacements without stretching (for plotting)
disps = Apos_map - Bpos

#fig = plt.figure()
#ax = fig.add_subplot(111)
ax.quiver(Bpos.T[:,0], Bpos.T[:,1], disps.T[:,0], disps.T[:,1], scale_units='xy', scale=1)
maxXAxis = np.max([Apos.max(), Bpos.max()]) + 1
ax.set_xlim([-maxXAxis, maxXAxis])
ax.set_ylim([-maxXAxis, maxXAxis])
ax.set_aspect('equal')
fig.savefig('DispLattice.png')

# Plotting the Apos and Bposst overlayed
fig = plt.figure()
ax = fig.add_subplot(111)
ax.scatter(Apos[:,:int(frac*125)].T[:,0],Apos[:,:int(frac*125)].T[:,1], c='C0')
ax.scatter(Apos[:,int(frac*125):].T[:,0],Apos[:,int(frac*125):].T[:,1], alpha=0.2, c="C0")
ax.scatter(Bposst.T[:,0],Bposst.T[:,1], alpha=0.5, c="C1")
maxXAxis = np.max([Apos.max(), Bposst.max()]) + 1
ax.set_xlim([-maxXAxis, maxXAxis])
ax.set_ylim([-maxXAxis, maxXAxis])
ax.set_aspect('equal')

# Displacement with stretching
disps = Apos_map - Bposst

#fig = plt.figure()
#ax = fig.add_subplot(111)
ax.quiver(Bposst.T[:,0], Bposst.T[:,1],disps.T[:,0], disps.T[:,1], scale_units='xy', scale=1)
maxXAxis = np.max([Apos.max(), Bposst.max()]) + 1
ax.set_xlim([-maxXAxis, maxXAxis])
ax.set_ylim([-maxXAxis, maxXAxis])
ax.set_aspect('equal')
fig.savefig('DispLattice_stretched.png')

# Stretching Matrix
Acell3d = np.identity(3) # TMP only for 2D
Acell3d[:2,:2] = Acell
Acell = Acell3d

Bcell3d = np.identity(3) # TMP only for 2D
Bcell3d[:2,:2] = Bcell
Bcell = Bcell3d

# Stretching Matrix
stMat = la.inv(tr.canonicalize(Bcell)).dot(tr.canonicalize(tmat.dot(Bcell)))

# Rotation Matrix
rtMat = tmat.dot(Bcell).dot(la.inv(stMat)).dot(la.inv(Bcell))

print("Stretching Matrix:")
print(stMat)

print("Rotation Matrix:")
print(rtMat)


class_list, vec_classes = classify(disps)

# Display the diffrent classes of displacement 
for i in range(len(vec_classes)):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    disps_class = disps[:,class_list==i]
    Bposst_class = Bposst[:,class_list==i]
    ax.quiver(Bposst_class.T[:,0], Bposst_class.T[:,1], disps_class.T[:,0], disps_class.T[:,1], scale_units='xy', scale=1)
    maxXAxis = np.max([Apos.max(), Bposst.max()]) + 1
    ax.set_xlim([-maxXAxis, maxXAxis])
    ax.set_ylim([-maxXAxis, maxXAxis])
    ax.set_aspect('equal')
    fig.savefig('DispLattice_stretched_%d.png'%i)

# Centers the position on the first atom
pos_in_struc = Bposst- Bposst[:,0:1].dot(np.ones((1,np.shape(Bposst)[1])))

cell = find_cell(class_list, Bposst)

# Finds a squarer cell
cell = gruber(cell)
cell = cell[:,np.argsort(cell[2,:])] # Only in 2D so that the z component is last 

# Make a pylada structure
cell_coord = la.inv(cell).dot(pos_in_struc)
idx_struc = np.where(np.sum((cell_coord < 1-tol) & (cell_coord > - tol), axis = 0 ) == 3)[0]
Struc = Structure(cell)
for i, disp_type in enumerate(class_list[idx_struc]):
    Struc.add_atom(*(tuple(pos_in_struc[:,idx_struc[i]])+(str(disp_type),)))
    
# Makes sure it is the primitive cell 
Struc = primitive(Struc, tolerance = tol)

# Total displacement per unit volume a as metric
Total_disp = 0 
for disp in Struc:
    Total_disp += la.norm(vec_classes[int(disp.type)])

Total_disp = Total_disp / la.det(Struc.cell)

cell = Struc.cell

print("Displacement Lattice")
print(cell)

print("Volume stretching factor:", la.det(stMat))
print("Total displacement stretched cell:", Total_disp)

# Displays displacement with the disp cell overlayed
fig = plt.figure()
ax = fig.add_subplot(111)
ax.quiver(pos_in_struc.T[:,0], pos_in_struc.T[:,1],disps.T[:,0], disps.T[:,1], scale_units='xy', scale=1)
ax.quiver(np.zeros(2), np.zeros(2), cell[0,:2], cell[1,:2], scale_units='xy', scale=1)
maxXAxis = pos_in_struc.max() + 1
ax.set_xlim([-maxXAxis, maxXAxis])
ax.set_ylim([-maxXAxis, maxXAxis])
ax.set_aspect('equal')
fig.savefig('DispLattice_stretched_cell_primittive.png')

# Displays only the cell and the displacements in it
fig = plt.figure()
ax = fig.add_subplot(111)
for disp in Struc:
    ax.quiver(disp.pos[0],disp.pos[1], vec_classes[int(disp.type)][0],vec_classes[int(disp.type)][1], scale_units='xy', scale=1)
ax.quiver(np.zeros(2), np.zeros(2), cell[0,:2], cell[1,:2], scale_units='xy', scale=1, color = "blue", alpha = 0.3)
maxXAxis = cell.max() + 1
ax.set_xlim([-maxXAxis, maxXAxis])
ax.set_ylim([-maxXAxis, maxXAxis])
ax.set_aspect('equal')
fig.savefig('Displacement_structure.png')

plt.close('All')
    


