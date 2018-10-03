from p2ptrans import transform as tr
import numpy as np
import numpy.linalg as la
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
from p2ptrans import tiling as t
import pickle
import time
from pylada.crystal import Structure, primitive, gruber, supercell
from copy import deepcopy
import warnings

tol = 1e-2

def uniqueclose(closest, tol):
    unique = []
    idx = []
    for i,line in enumerate(closest.T):
        there = False
        for check in unique:
            if np.allclose(check, line, atol=tol):
                there = True
                break
        if not there:
            unique.append(line)
            idx.append(i)
    return (np.array(idx), np.array(unique))

def find_cell(class_list, positions, tol = 1e-5, frac_tol = 0.5):
    cell_list = []
    origin_list = []
    for i in np.unique(class_list):
        newcell = np.identity(3)
        pos = positions[:, class_list == i]
        center = np.argmin(la.norm(pos, axis = 0))
        list_in = list(range(np.shape(pos)[1]))
        list_in.remove(center)
        origin = pos[:,center:center+1]
        pos = pos[:,list_in] - origin.dot(np.ones((1,np.shape(pos)[1]-1))) # centered
        norms = la.norm(pos, axis = 0)
        idx = np.argsort(norms)
        j = 0;
        for k in idx:
            multiple = [0]
            #If there is already one cell vector (j=1) skips the candidate if it's parallel
            if j == 1: 
                if la.norm(np.cross(pos[:,k], newcell[:,0])) < tol:
                    continue
            # If there is already two cell vectors (j=2) skips the candidate if it's parallel
            # to one of the vectors
            elif j == 2:
                if abs(la.det(np.concatenate((newcell[:,:2],pos[:,k:k+1]), axis=1))) < tol:
                    continue
            # Goes through the displacements and finds the one that are parallel
            for p in pos.T:
                if la.norm(np.cross(p,pos[:,k]))/(la.norm(p) * la.norm(pos[:,k])) < tol:
                    multiple.append(p.dot(pos[:,k])/la.norm(pos[:,k])**2)

            # Find the norms of all vectors with respect to the center of the interval
            # finds all the displacements inside the 'shell' of the interval
            if multiple != []:
                multiple.sort()
                norms = la.norm(pos - 1 / 2.0 * (multiple[0]+multiple[-1]) *
                                pos[:,k:k+1].dot(np.ones((1,np.shape(pos)[1]))), axis = 0)
                shell = np.sum(norms < 1 / 2.0 * (multiple[-1]-multiple[0])*la.norm(pos[:,k]) + tol) / float(len(norms))
                # If it is the right size (to check next condition)
                if len(multiple) == len(np.arange(round(multiple[0]), round(multiple[-1])\
    +1)):
                    # If all the multiples are present and the interval cover more than a certain 
                    # fraction of displacements the displacement is added as a cell vector
                    if np.allclose(multiple, np.arange(round(multiple[0]), round(multiple[-1])+1), tol) and shell > frac_tol**3:
                        newcell[:,j] = pos[:,k]
                        j += 1
                        if j == 2:
                            for cell in cell_list:
                                if abs(la.det(cell)) < abs(la.det(newcell)):
                                    if not np.allclose(la.inv(cell).dot(newcell), np.round(la.inv(cell).dot(newcell)), tol):
                                        raise RuntimeError("The periodicity of the different classes of displacement is different")
                                else:
                                    if not np.allclose(la.inv(newcell).dot(cell), np.round(la.inv(newcell).dot(cell)),\
     tol):
                                        raise RuntimeError("The periodicity of the different classes of displacement is different")

                            # if la.det(newcell) < 0:
                            #     newcell[:,2] = -newcell[:,2]

                            cell_list.append(newcell)
                            origin_list.append(origin)
                            break
        else:
            warnings.warn("Could not find periodic cell for displacement %d. Increase sample size or use results with care."%i, RuntimeWarning)
            
    if len(cell_list) == 0:
        raise RuntimeError("Could not find periodic cell for any displacement. Increase sample size.")

    cell = cell_list[np.argmax([la.det(cell) for cell in cell_list])]
    origin = origin_list[np.argmax([la.det(cell) for cell in cell_list])]

    return cell, origin

def classify(disps, tol = 1.e-1):
    vec_classes = [disps[:,0:1]]
    class_list = np.zeros(np.shape(disps)[1], np.int)
    for i in range(np.shape(disps)[1]):
        classified = False
        for j, vec_class in enumerate(vec_classes):
            vec_mean = np.mean(vec_class, axis=1)
            # if (abs(la.norm(vec_mean) - vec_mean.T.dot(disps[:,i])/la.norm(vec_mean)) < tol and 
            #     la.norm(np.cross(vec_mean, disps[:,i]))/la.norm(vec_mean) < tol):
            if (la.norm(vec_mean - disps[:,i]) < tol):
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

def lcm(x, y):
   """This function takes two
   integers and returns the L.C.M."""

   lcm = (x*y)//gcd(x,y)
   return lcm

def gcd(x, y):
   """This function implements the Euclidian algorithm
   to find G.C.D. of two numbers"""

   while(y):
       x, y = y, x % y

   return x

# \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
# Begin Program

random = False

# Eventually this can be from pcread

# Multiple atoms example
# A = Structure(np.identity(3))
# A.add_atom(0,0,0,'Si')
# A.add_atom(0.1,0.1,0,'Be')

# B = Structure(np.array([[1.1,0,0],[1.1,1.5,0],[0,0,1]]).T)
# B.add_atom(0,0,0,'Si')
# B.add_atom(0.6,0.3,0,'Be')

# # The calssic
# A = Structure(np.identity(3))
# A.add_atom(0,0,0,'Si')

# B = Structure(np.array([[-1.0,1.0,0],[0.5,0.5,0],[0,0,1]]).T)
# B.add_atom(0,0,0,'Si')

# Classification classic

A = Structure(np.array([[19.9121208191,0,0],[9.9560632706,17.2444074280,0],[0,0,1]]).T * 0.25)
A.add_atom(0,0,0,'1')
Alabel = "Ni"

B = Structure(np.array([[7.27612877841,0,0],[3.6380643892,6.30131236331,0],[0,0,1]]).T*0.5)
B.add_atom(0,0,0,'1')
Blabel = "Zr"


mul = lcm(len(A),len(B))
mulA = mul//len(A)
mulB = mul//len(B)

if mulA*la.det(A.cell) < mulB*la.det(B.cell):
    tmp = deepcopy(B)
    tmpmul = mulB
    tmplabel = Blabel

    B = deepcopy(A)
    mulB = mulA
    Blabel = Alabel

    A = tmp
    mulA = tmpmul
    Alabel = tmplabel

ncell = 400

# Setting the unit cells of A and B
Acell = A.cell[:2,:2]
Bcell = B.cell[:2,:2] 

# Plotting the cell vectors of A and B
fig = plt.figure()
ax = fig.add_subplot(111)
ax.quiver(np.ones(3), np.ones(3), Bcell[0,:], Bcell[1,:])
ax.quiver(-np.ones(3), -np.ones(3), Acell[0,:], Acell[1,:])
ax.set_xlim([-5, 5])
ax.set_ylim([-5, 5])
fig.savefig('CellVectors.svg')


ASC = t.circle(Acell, mulA * ncell)
BSC = t.circle(Bcell, mulB * ncell)

centerA = np.mean(ASC, axis=1)
centerB = np.mean(BSC, axis=1)

print("Center A", centerA)
print("Center B", centerB)

# Plot gamma points of each A cell
fig = plt.figure()
ax = fig.add_subplot(111)
ax.scatter(ASC[0,:], ASC[1,:], c='k')
maxXAxis = np.max([ASC[0,:].max(),ASC[1,:].max()])
minXAxis = np.min([ASC[0,:].min(),ASC[1,:].min()])
ax.set_xlim([minXAxis-1, maxXAxis+1])
ax.set_ylim([minXAxis-1, maxXAxis+1])
ax.set_aspect('equal')
fig.savefig('Agrid.svg')

# Plot gamma points of each B cell
fig = plt.figure()
ax = fig.add_subplot(111)
ax.scatter(BSC[0,:], BSC[1,:], facecolor='w', edgecolor='k')
maxXAxis = np.max([BSC[0,:].max(),BSC[1,:].max()])
minXAxis = np.min([BSC[0,:].min(),BSC[1,:].min()])
ax.set_xlim([minXAxis-1, maxXAxis+1])
ax.set_ylim([minXAxis-1, maxXAxis+1])
ax.set_aspect('equal')
fig.savefig('Bgrid.svg')

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
    Apos = []
    atom_types = np.array([], np.str)
    atomsA = np.array([], np.int)
    for a in A:
        if any(atom_types == a.type):
            idx = np.where(atom_types == a.type)[0][0]
            Apos[idx] = np.concatenate((Apos[idx], ASC + Acell.dot(np.reshape(a.pos[:2],(2,1))).dot(np.ones((1,np.shape(ASC)[1])))), axis = 1) 
            atomsA[idx] += 1
        else:
            Apos.append(ASC + Acell.dot(np.reshape(a.pos[:2],(2,1))).dot(np.ones((1,np.shape(ASC)[1]))))
            atom_types = np.append(atom_types, a.type)
            atomsA = np.append(atomsA,1)

    Apos = np.concatenate(Apos, axis=1)

    Bpos = [None]*len(atom_types)
    atomsB = np.zeros(len(atom_types), np.int)
    for a in B:
        print(type(atom_types), atom_types)
        idx = np.where(atom_types == a.type)[0][0]
        if atomsB[idx] == 0:
            Bpos[idx] = BSC + Bcell.dot(np.reshape(a.pos[:2],(2,1))).dot(np.ones((1,np.shape(BSC)[1])))
        else:
            Bpos[idx] = np.concatenate((Bpos[idx], BSC + Bcell.dot(np.reshape(a.pos[:2],(2,1))).dot(np.ones((1,np.shape(BSC)[1])))), axis = 1) 
        atomsB[idx] += 1

    Bpos = np.concatenate(Bpos, axis=1)
        
    Apos = np.concatenate([Apos,np.zeros((1,np.shape(Apos)[1]))]) # TMP Only in 2D Adds the 3rd dim. 
    Bpos = np.concatenate([Bpos,np.zeros((1,np.shape(Bpos)[1]))])

    assert all(mulA*atomsA == mulB*atomsB)
    atoms = mulA*atomsA

fracB = 0.25
fracA = 0 # Set to 0 to allow skipping
Acell_tmp = np.identity(3)
Acell_tmp[:2,:2] = Acell

Apos = np.asfortranarray(Apos)
Bpos = np.asfortranarray(Bpos) 
t_time = time.time()
Apos_map, Bpos, Bposst, n_map, ttrans, rtrans, dmin = tr.fastoptimization(Apos, Bpos, fracA, fracB, Acell_tmp, la.inv(Acell_tmp), atoms, 1000, 100, 4, 4, 1e-5, 1e-5) #TMP
t_time = time.time() - t_time
Bpos = np.asanyarray(Bpos)
Apos = np.asanyarray(Apos)

print(dmin)
print("Mapping time:", t_time)

pickle.dump((Apos_map, Bpos, Bposst, n_map, ttrans, rtrans, dmin), open("fastoptimization.dat","wb"))

# # TMP for testing only -->
# tr.center(Apos)
# tr.center(Bpos)
# Apos_map, Bpos, Bposst, n_map, ttrans, rtrans, dmin = pickle.load(open("fastoptimization.dat","rb"))
# # <--  

tmat = ttrans[:,:3]

Bpos = Bpos[:,:n_map]
Bposst = Bposst[:,:n_map]
Apos_map = Apos_map[:,:n_map]

natB = np.shape(Bposst)[1] // np.sum(atoms)
nat = np.shape(Apos)[1] // np.sum(atoms)
natA = int(fracA*np.shape(Apos)[1]/np.sum(atoms))


# Plotting the Apos and Bpos overlayed
fig = plt.figure()
ax = fig.add_subplot(111)
#ax.scatter(Apos.T[:,0],Apos.T[:,1])
for i,num in enumerate(atoms):
    for j in range(num):
        ax.scatter(Apos.T[(np.sum(atoms[:i-1])+j)*nat:(np.sum(atoms[:i-1])+j)*nat + natA,0],Apos.T[(np.sum(atoms[:i-1])+j)*nat:(np.sum(atoms[:i-1])+j)*nat + natA,1], c="C%d"%(2*i))
        ax.scatter(Apos.T[(np.sum(atoms[:i-1])+j)*nat + natA:(np.sum(atoms[:i-1])+j+1)*nat,0],Apos.T[(np.sum(atoms[:i-1])+j)*nat + natA:(np.sum(atoms[:i-1])+j+1)*nat,1], c="C%d"%(2*i), alpha = 0.5)
    ax.scatter(Bpos.T[natB*num*i:natB*num*(i+1),0],Bpos.T[natB*num*i:natB*num*(i+1),1], alpha=0.5, c="C%d"%(2*i+1))
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
fig.savefig('DispLattice.svg')

# Plotting the Apos and Bposst overlayed
fig = plt.figure()
ax = fig.add_subplot(111)
#ax.scatter(Apos.T[:,0],Apos.T[:,1])
for i,num in enumerate(atoms):
    for j in range(num):
        ax.scatter(Apos.T[(np.sum(atoms[:i-1])+j)*nat:(np.sum(atoms[:i-1])+j)*nat + natA,0],Apos.T[(np.sum(atoms[:i-1])+j)*nat:(np.sum(atoms[:i-1])+j)*nat + natA,1], c="C%d"%(2*i), label=Alabel)
        ax.scatter(Apos.T[(np.sum(atoms[:i-1])+j)*nat + natA:(np.sum(atoms[:i-1])+j+1)*nat,0],Apos.T[(np.sum(atoms[:i-1])+j)*nat + natA:(np.sum(atoms[:i-1])+j+1)*nat,1], c="C%d"%(2*i), alpha=0.5)
    ax.scatter(Bposst.T[natB*num*i:natB*num*(i+1),0],Bposst.T[natB*num*i:natB*num*(i+1),1], alpha=0.5, c="C%d"%(2*i+1), label=Blabel)
ax.legend()
maxXAxis = np.max([Apos.max(), Bposst.max()]) + 1
ax.set_xlim([-maxXAxis, maxXAxis])
ax.set_ylim([-maxXAxis, maxXAxis])
ax.set_aspect('equal')

# Displacement with stretching
disps = Apos_map - Bposst


# fig = plt.figure()
# ax = fig.add_subplot(111)
ax.quiver(Bposst.T[:,0], Bposst.T[:,1],disps.T[:,0], disps.T[:,1],scale_units='xy', scale=1)
maxXAxis = np.max([Apos.max(), Bposst.max()]) + 1
ax.set_xlim([-maxXAxis, maxXAxis])
ax.set_ylim([-maxXAxis, maxXAxis])
ax.set_aspect('equal')
fig.savefig('DispLattice_stretched.svg')


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

print("vec_classes", vec_classes)


# Only the displacements
fig = plt.figure()
ax = fig.add_subplot(111)
maxXAxis = np.max(disps) + 1
ax.set_xlim([-maxXAxis, maxXAxis])
ax.set_ylim([-maxXAxis, maxXAxis])
ax.set_aspect('equal')
for i in range(len(vec_classes)):
    disps_class = disps[:,class_list==i]
    ndisps = np.shape(disps_class)[1]
    ax.quiver(np.zeros((1,ndisps)), np.zeros((1,ndisps)), disps_class.T[:,0], disps_class.T[:,1], color="C%d"%(i%10), scale_units='inches', scale=1)
fig.savefig('DispOverlayed.svg')

# Centers the position on the first atom
pos_in_struc = Bposst-Bposst[:,0:1].dot(np.ones((1,np.shape(Bposst)[1])))

print("Volume stretching factor:", la.det(tmat))
print("Cell volume ratio (should be exactly the same):", mulA * la.det(Acell)/(mulB * la.det(Bcell)))

# plt.show()

cell, origin = find_cell(class_list, Bposst)

pos_in_struc = Bposst - origin.dot(np.ones((1,np.shape(Bposst)[1])))
posA_in_struc = Apos - origin.dot(np.ones((1,np.shape(Apos)[1])))



# Make a pylada structure
DispStruc = Structure(cell)
FinalStruc = Structure(cell)

cell_coord = np.mod(la.inv(cell).dot(pos_in_struc)+tol,1)-tol

for i, disp in zip(*uniqueclose(cell_coord, tol)):
    DispStruc.add_atom(*(tuple(cell.dot(disp))+(str(class_list[i]),)))
    FinalStruc.add_atom(*(tuple(cell.dot(disp))+(Blabel,)))

cell_coord = np.mod(la.inv(cell).dot(posA_in_struc)+tol,1)-tol

for i, disp in zip(*uniqueclose(cell_coord, tol)):
    FinalStruc.add_atom(*(tuple(cell.dot(disp))+(Alabel,)))
    
if la.det(cell) < 0:
    cell[:,2] = -cell[:,2] 

# Finds a squarer cell
print("cell before gruber", cell)
cell = gruber(cell)
cell = cell[:,np.argsort(abs(cell[2,:]))] # Only in 2D so that the z component is last
cell[1,:] = cell[2,2]*cell[1,:]
cell[2,2] = abs(cell[2,2])

DispStruc = supercell(DispStruc, cell)

# Makes sure it is the primitive cell 
DispStruc = primitive(DispStruc, tolerance = tol)

FinalStruc = supercell(FinalStruc, DispStruc.cell)

# Total displacement per unit volume a as metric
Total_disp = 0 
for disp in DispStruc:
    Total_disp += la.norm(vec_classes[int(disp.type)])

Total_disp = Total_disp / la.det(DispStruc.cell)

cell = DispStruc.cell

pickle.dump((Alabel, Blabel, rtrans, ttrans, cell, centerA, centerB), open("transformation.dat","wb")) 

print("Displacement Lattice")
print(cell)

print("Volume stretching factor:", la.det(stMat))
print("Total displacement stretched cell:", Total_disp)

# Displays displacement with the disp cell overlayed
fig = plt.figure()
ax = fig.add_subplot(111)
ax.quiver(pos_in_struc.T[:,0], pos_in_struc.T[:,1],disps.T[:,0], disps.T[:,1], scale_units='xy', scale=1)
ax.quiver(np.zeros(2), np.zeros(2), cell[0,:2], cell[1,:2], scale_units='xy', scale=1, color="red")
maxXAxis = pos_in_struc.max() + 1
ax.set_xlim([-maxXAxis, maxXAxis])
ax.set_ylim([-maxXAxis, maxXAxis])
ax.set_aspect('equal')
fig.savefig('DispLattice_stretched_cell_primittive.svg')

print("DISPS", DispStruc)
print("A", FinalStruc)

# Displays only the cell and the displacements in it
fig = plt.figure()
ax = fig.add_subplot(111)
for i,disp in enumerate(DispStruc):
    ax.quiver(disp.pos[0],disp.pos[1], vec_classes[int(disp.type)][0],vec_classes[int(disp.type)][1], scale_units='xy', scale=1)

for i,a in enumerate(FinalStruc):
    if a.type == Alabel:
        ax.scatter(a.pos[0], a.pos[1], color="C0")
    else:
        ax.scatter(a.pos[0], a.pos[1], color="C1")
    
ax.quiver(np.zeros(2), np.zeros(2), cell[0,:2], cell[1,:2], scale_units='xy', scale=1, color = "blue", alpha = 0.3)
maxXAxis = abs(cell).max() + 1
ax.set_xlim([-maxXAxis, maxXAxis])
ax.set_ylim([-maxXAxis, maxXAxis])
ax.set_aspect('equal')
fig.savefig('Displacement_structure.svg')

plt.show()

plt.close('All')
    


