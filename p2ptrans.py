from p2ptrans import transform as tr
import numpy as np
import numpy.linalg as la
from mpl_toolkits.mplot3d import Axes3D
# import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import animation
from p2ptrans import tiling as t
import pickle
import time
from pylada.crystal import Structure, primitive, gruber, read
from copy import deepcopy

tol = 1e-5

def find_cell(class_list, positions, tol = 1e-5, frac_tol = 0.5):
    cell = np.zeros((3,3))
    for i in range(len(np.unique(class_list))):
        newcell = np.identity(3)
        pos = positions[:, class_list == i]
        center = np.argmin(la.norm(pos, axis = 0))
        list_in = list(range(np.shape(pos)[1]))
        list_in.remove(center)
        pos = pos[:,list_in] - pos[:,center:center+1].dot(np.ones((1,np.shape(pos)[1]-1))) # centered
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
                norms = la.norm(pos - 1 / 2.0 * (multiple[0]+multiple[-1]) *
                                pos[:,k:k+1].dot(np.ones((1,np.shape(pos)[1]))), axis = 0)
                shell = np.sum(norms < 1 / 2.0 * (multiple[-1]-multiple[0])*la.norm(pos[:,k])) / float(len(norms))
                if len(multiple) == len(np.arange(round(multiple[0]), round(multiple[-1])\
    +1)):
                    if np.allclose(multiple, np.arange(round(multiple[0]), round(multiple[-1])+1), tol) and shell > frac_tol:
                        newcell[:,j] = pos[:,k]
                        j += 1
                        if j == 3: break 
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

def classify(disps, tol = 1.e-3):
    vec_classes = [disps[:,0:1]]
    class_list = np.zeros(np.shape(disps)[1], np.int)
    for i in range(np.shape(disps)[1]):
        classified = False
        for j, vec_class in enumerate(vec_classes):
            vec_mean = np.mean(vec_class, axis=1)
            if (abs(la.norm(vec_mean) - vec_mean.T.dot(disps[:,i])/la.norm(vec_mean)) < tol and 
                la.norm(np.cross(vec_mean, disps[:,i]))/la.norm(vec_mean) < tol):
                vec_classes[j] = np.concatenate((vec_class,disps[:,i:i+1]),axis=1)
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

def animate(i):
    if i<180:
        ax.view_init(30,i)
    elif i<240:
        ax.view_init(30,360-i)
    elif i<300:
        ax.view_init(i-210,120)
    else:
        ax.view_init(390-i,120)
    return fig,

# \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
# Begin Program

random = False

# Eventually this can be from pcread
A = Structure(np.identity(3))
A.add_atom(0,0,0,'Si')

B = Structure(np.array([[1,1,0],[-0.5,0.5,0],[0,0,1]]).T)
B.add_atom(0,0,0,'Si')


mul = lcm(len(A),len(B))
mulA = mul//len(A)
mulB = mul//len(B)

if mulA*la.det(A.cell) < mulB*la.det(B.cell):
    tmp = deepcopy(B)
    tmpmul = mulB
    B = deepcopy(A)
    mulB = mulA
    A = tmp
    mulA = tmpmul

ncell = 300

Acell = A.cell
Bcell = B.cell

# Plotting the cell vectors of A and B
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.quiver(np.ones(3), np.ones(3), np.ones(3), Bcell[0,:], Bcell[1,:], Bcell[2,:])
ax.quiver(-np.ones(3), -np.ones(3), -np.ones(3), Acell[0,:], Acell[1,:], Acell[2,:])
ax.set_xlim([-5, 5])
ax.set_ylim([-5, 5])
ax.set_zlim([-5, 5])
fig.savefig('CellVectors.png')


ASC = t.sphere(Acell, mulA * ncell)
BSC = t.sphere(Bcell, mulB * ncell)

# Plot gamma points of each A cell
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(ASC[0,:], ASC[1,:], ASC[2,:])
maxXAxis = ASC.max() + 1
minXAxis = BSC.min() - 1
ax.set_xlim([minXAxis-1, maxXAxis+1])
ax.set_ylim([minXAxis-1, maxXAxis+1])
ax.set_zlim([minXAxis-1, maxXAxis+1])
ax.set_aspect('equal')
fig.savefig('Agrid.png')

# Plot gamma points of each B cell
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(BSC[0,:], BSC[1,:], BSC[2,:])
maxXAxis = BSC.max() + 1
minXAxis = BSC.min() - 1
ax.set_xlim([minXAxis-1, maxXAxis+1])
ax.set_ylim([minXAxis-1, maxXAxis+1])
ax.set_zlim([minXAxis-1, maxXAxis+1])
ax.set_aspect('equal')
fig.savefig('Bgrid.png')

# For testing purposes
if random:
    # Create a random Apos and B with random small displacement
    atoms = np.array([1]) # One atom
    n = 10
    Apos = np.random.random((3,n))*3
    
    # Transform Apos to get Bpos
    angles = 2*np.pi*np.random.random(3)
    vec = np.random.random(3)
    
    print("ANGLES:", angles)
    print("VEC:", vec)

    ttmat = np.array([[1,0,0],[0,2,0],[0,0,1]])
    Bpos = np.asfortranarray((np.array(Apos).T.dot(ttmat)).T)

    tr.trans(Bpos,angles,vec)

    randDisp = np.random.random((3,n))
    
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
            Apos[idx] = np.concatenate((Apos[idx], ASC + Acell.dot(np.reshape(a.pos,(3,1))).dot(np.ones((1,np.shape(ASC)[1])))), axis = 1) 
            atomsA[idx] += 1
        else:
            Apos.append(ASC + Acell.dot(np.reshape(a.pos,(3,1))).dot(np.ones((1,np.shape(ASC)[1]))))
            atom_types = np.append(atom_types, a.type)
            atomsA = np.append(atomsA,1)

    Apos = np.concatenate(Apos, axis=1)

    Bpos = [None]*len(atom_types)
    atomsB = np.zeros(len(atom_types), np.int)
    for a in B:
        idx = np.where(atom_types == a.type)[0][0]
        if atomsB[idx] == 0:
            Bpos[idx] = BSC + Bcell.dot(np.reshape(a.pos,(3,1))).dot(np.ones((1,np.shape(BSC)[1])))
        else:
            Bpos[idx] = np.concatenate((Bpos[idx], BSC + Bcell.dot(np.reshape(a.pos,(3,1))).dot(np.ones((1,np.shape(BSC)[1])))), axis = 1) 
        atomsB[idx] += 1

    Bpos = np.concatenate(Bpos, axis=1)
        
    assert all(mulA*atomsA == mulB*atomsB)
    atoms = mulA*atomsA

fracB = 0.4
fracA = 0.15 # fracA < fracB

# TMP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
# Slanting
tt = 30/360.*2*np.pi
ph = 90/360.*2*np.pi
sg = 90/360.*2*np.pi

v2 = np.array([np.cos(ph + np.pi/2)*np.cos(tt), np.cos(ph + np.pi/2)*np.sin(tt), np.sin(ph + np.pi/2)])
v3 = np.array([np.cos(tt + np.pi/2),np.sin(tt + np.pi/2),0])
vM = np.array([[np.cos(ph)*np.cos(tt), np.cos(ph)*np.sin(tt), np.sin(ph)], 
               np.cos(sg)*v2 + np.sin(sg)*v3,
               np.cos(sg + np.pi/2)*v2 + np.sin(sg + np.pi/2)*v3])
k = 1.5
tM = np.array([[1, k, 0], [0,1,0], [0,0,1]])

Aslant = la.inv(vM).dot(tM).dot(vM).dot(Apos)

fig = plt.figure(31)
ax = fig.add_subplot(111, projection='3d')
ax.scatter(Aslant[0,:], Aslant[1,:], Aslant[2,:])
ax.quiver([0],[0],[0],vM[0], vM[1], vM[2])
maxXAxis = np.max(Aslant.max()) + 1
ax.set_xlim([-maxXAxis, maxXAxis])
ax.set_ylim([-maxXAxis, maxXAxis])
ax.set_zlim([-maxXAxis, maxXAxis])
ax.set_aspect('equal')
ax.view_init(-90,0)

fig = plt.figure(32)
ax = fig.add_subplot(111, projection='3d')
ax.scatter(Apos[0,:], Apos[1,:], Apos[2,:])
ax.view_init(-90,0)

plt.show()
raise

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



Apos = np.asfortranarray(Apos)
Bpos = np.asfortranarray(Bpos) 
t_time = time.time()
# Apos_map, Bpos, Bposst, n_map, tmat, dmin = tr.fastoptimization(Apos, Bpos, fracA, fracB, Acell, la.inv(Acell), atoms, 1, 300, 5, 5, 1e-6, 1e-6)
t_time = time.time() - t_time
Bpos = np.asanyarray(Bpos)
Apos = np.asanyarray(Apos)

print("Mapping time:", t_time)

pickle.dump((Apos_map, Bpos, Bposst, n_map, tmat, dmin), open("fastoptimization.dat","wb"))

# # TMP for testing only -->
# tr.center(Apos)
# tr.center(Bpos)
# Apos_map, Bpos, Bposst, n_map, tmat, dmin = pickle.load(open("fastoptimization.dat","rb"))
# # <--  

print("Total distance between structures:", dmin)

Bpos = Bpos[:,:n_map]
Bposst = Bposst[:,:n_map]
Apos_map = Apos_map[:,:n_map]

natB = np.shape(Bposst)[1] // np.sum(atoms)
nat = np.shape(Apos)[1] // np.sum(atoms)
natA = int(fracA*np.shape(Apos)[1]/np.sum(atoms))


# Plotting the Apos and Bpos overlayed
fig = plt.figure(22)
ax = fig.add_subplot(111, projection='3d')
#ax.scatter(Apos.T[:,0],Apos.T[:,1])
for i,num in enumerate(atoms):
    for j in range(num):
        ax.scatter(Apos.T[(np.sum(atoms[:i-1])+j)*nat:(np.sum(atoms[:i-1])+j)*nat + natA,0],Apos.T[(np.sum(atoms[:i-1])+j)*nat:(np.sum(atoms[:i-1])+j)*nat + natA,1],Apos.T[(np.sum(atoms[:i-1])+j)*nat:(np.sum(atoms[:i-1])+j)*nat + natA,2], c="C%d"%(2*i))
        ax.scatter(Apos.T[(np.sum(atoms[:i-1])+j)*nat + natA:(np.sum(atoms[:i-1])+j+1)*nat,0],Apos.T[(np.sum(atoms[:i-1])+j)*nat + natA:(np.sum(atoms[:i-1])+j+1)*nat,1], Apos.T[(np.sum(atoms[:i-1])+j)*nat + natA:(np.sum(atoms[:i-1])+j+1)*nat,2], c="C%d"%(2*i), alpha = 0.5)
    ax.scatter(Bpos.T[natB*num*i:natB*num*(i+1),0],Bpos.T[natB*num*i:natB*num*(i+1),1], Bpos.T[natB*num*i:natB*num*(i+1),2], alpha=0.5, c="C%d"%(2*i+1))
maxXAxis = np.max([Apos.max(), Bpos.max()]) + 1
ax.set_xlim([-maxXAxis, maxXAxis])
ax.set_ylim([-maxXAxis, maxXAxis])
ax.set_zlim([-maxXAxis, maxXAxis])
ax.set_aspect('equal')


# Displacements without stretching (for plotting)
disps = Apos_map - Bpos

#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
ax.quiver(Bpos.T[:,0], Bpos.T[:,1], Bpos.T[:,2], disps.T[:,0], disps.T[:,1], disps.T[:,2])
maxXAxis = np.max([Apos.max(), Bpos.max()]) + 1
ax.set_xlim([-maxXAxis, maxXAxis])
ax.set_ylim([-maxXAxis, maxXAxis])
ax.set_zlim([-maxXAxis, maxXAxis])
ax.set_aspect('equal')
fig.savefig('DispLattice.png')

# Displacement with stretching
disps = Apos_map - Bposst

fig = plt.figure()
ax = Axes3D(fig)
maxXAxis = np.max([Apos.max(), Bposst.max()]) + 1
ax.set_xlim([-maxXAxis, maxXAxis])
ax.set_ylim([-maxXAxis, maxXAxis])
ax.set_zlim([-maxXAxis, maxXAxis])
ax.set_aspect('equal')

# Plotting the Apos and Bposst overlayed
def init_disps():
    #ax.scatter(Apos.T[:,0],Apos.T[:,1])
    for i,num in enumerate(atoms):
        for j in range(num):
            ax.scatter(Apos.T[(np.sum(atoms[:i-1])+j)*nat:(np.sum(atoms[:i-1])+j)*nat + natA,0],Apos.T[(np.sum(atoms[:i-1])+j)*nat:(np.sum(atoms[:i-1])+j)*nat + natA,1],Apos.T[(np.sum(atoms[:i-1])+j)*nat:(np.sum(atoms[:i-1])+j)*nat + natA,2], c="C%d"%(2*i))
            ax.scatter(Apos.T[(np.sum(atoms[:i-1])+j)*nat + natA:(np.sum(atoms[:i-1])+j+1)*nat,0],Apos.T[(np.sum(atoms[:i-1])+j)*nat + natA:(np.sum(atoms[:i-1])+j+1)*nat,1],Apos.T[(np.sum(atoms[:i-1])+j)*nat + natA:(np.sum(atoms[:i-1])+j+1)*nat,2], c="C%d"%(2*i), alpha=0.5)
        ax.scatter(Bposst.T[natB*num*i:natB*num*(i+1),0],Bposst.T[natB*num*i:natB*num*(i+1),1], Bposst.T[natB*num*i:natB*num*(i+1),2], alpha=0.5, c="C%d"%(2*i+1))
    
    ax.quiver(Bposst.T[:,0], Bposst.T[:,1], Bposst.T[:,2], disps.T[:,0], disps.T[:,1], disps.T[:,2])
    fig.savefig('DispLattice_stretched.png')
    return fig,

init_disps()

# anim = animation.FuncAnimation(fig, animate, init_func=init_disps,
#                                frames=490, interval=30)

# anim.save('Crystal+Disps.gif', fps=30, codec='gif')


# Stretching Matrix
stMat = la.inv(tr.canonicalize(Bcell)).dot(tr.canonicalize(tmat.dot(Bcell)))

# Rotation Matrix
rtMat = tmat.dot(Bcell).dot(la.inv(stMat)).dot(la.inv(Bcell))

print("Stretching Matrix:")
print(stMat)

print("Rotation Matrix:")
print(rtMat)

class_list, vec_classes = classify(disps)

# # Display the diffrent classes of displacement 
# for i in range(len(vec_classes)):
#     fig = plt.figure()
#     ax = fig.add_subplot(111, projection='3d')
#     disps_class = disps[:,class_list==i]
#     Bposst_class = Bposst[:,class_list==i]
#     ax.quiver(Bposst_class.T[:,0], Bposst_class.T[:,1], Bposst_class.T[:,2], disps_class.T[:,0], disps_class.T[:,1], disps_class.T[:,2], color="C0")
#     ax.scatter(Bposst_class.T[:,0], Bposst_class.T[:,1], Bposst_class.T[:,2], alpha = 0.5, s=10, color="C0")
#     maxXAxis = Bposst.max() + 1
#     ax.set_xlim([-maxXAxis, maxXAxis])
#     ax.set_ylim([-maxXAxis, maxXAxis])
#     ax.set_zlim([-maxXAxis, maxXAxis])
#     ax.set_aspect('equal')
#     fig.savefig('DispLattice_stretched_%d.png'%i)

# Centers the position on the first atom
pos_in_struc = Bposst- Bposst[:,0:1].dot(np.ones((1,np.shape(Bposst)[1])))

print("Showing")

plt.show()

cell = find_cell(class_list, Bposst)

# Finds a squarer cell
cell = gruber(cell)

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
ax = fig.add_subplot(111, projection='3d')
ax.quiver(pos_in_struc.T[:,0], pos_in_struc.T[:,1], pos_in_struc.T[:,2], disps.T[:,0], disps.T[:,1], disps.T[:,2], color = "C0")
ax.scatter(pos_in_struc.T[:,0], pos_in_struc.T[:,1], pos_in_struc.T[:,2], s=10, color = "C0")
ax.quiver(np.zeros(3), np.zeros(3), np.zeros(3), cell[0,:], cell[1,:], cell[2,:], color = "red")
maxXAxis = pos_in_struc.max() + 1
ax.set_xlim([-maxXAxis, maxXAxis])
ax.set_ylim([-maxXAxis, maxXAxis])
ax.set_zlim([-maxXAxis, maxXAxis])
ax.set_aspect('equal')
fig.savefig('DispLattice_stretched_cell_primittive.png')

# Displays only the cell and the displacements in it
fig = plt.figure()
ax = Axes3D(fig)
#ax = fig.add_subplot(111, projection='3d')

def init_struc():
    for i,disp in enumerate(Struc):
        ax.quiver(disp.pos[0], disp.pos[1], disp.pos[2], vec_classes[int(disp.type)][0],vec_classes[int(disp.type)][1], vec_classes[int(disp.type)][2], color="C%d"%i)
        ax.scatter(disp.pos[0], disp.pos[1], disp.pos[2], alpha = 0.5, s=10, color="C%d"%i)
    ax.quiver(np.zeros(3), np.zeros(3), np.zeros(3), cell[0,:], cell[1,:], cell[2,:], color = "red", alpha = 0.3)
    maxXAxis = cell.max() + 1
    ax.set_xlim([-maxXAxis, maxXAxis])
    ax.set_ylim([-maxXAxis, maxXAxis])
    ax.set_zlim([-maxXAxis, maxXAxis])
    ax.set_aspect('equal')
    fig.savefig('Displacement_structure.png')
    return fig,

anim = animation.FuncAnimation(fig, animate, init_func=init_struc,
                               frames=490, interval=30)

anim.save('DispStruc.gif', fps=30, codec='gif')
    
plt.show()

plt.close('All')
