from p2ptrans import transform as tr
import numpy as np
import numpy.linalg as la
from mpl_toolkits.mplot3d import Axes3D
import matplotlib
from matplotlib import animation
from p2ptrans import tiling as t
import pickle
import time
from pylada.crystal import Structure, primitive, gruber, read
from copy import deepcopy
import argparse
import os

# Tolerence for structure identification
tol = 1e-5

def readOptions():

    parser = argparse.ArgumentParser()
    parser.add_argument("-I","--initial",dest="A",type=str, default='./POSCAR_A', help="Initial Structure")
    parser.add_argument("-F","--final",dest="B",type=str, default='./POSCAR_B', help="Final Structure")
    parser.add_argument("-n","--ncell",dest="ncell",type=int, default=300, help="Number of cells to tile")
    parser.add_argument("-a","--fracA",dest="fracA",type=float, default=0.15, help="Fraction of the biggest structure to force to use in mapping (fracA < fracB)")
    parser.add_argument("-b","--fracB",dest="fracB",type=float, default=0.4, help="Fraction of the smallest structure to force to use in mapping")
    parser.add_argument("-r","--niter",dest="niter",type=int, default=10000, help="Number of (r)andom starts")
    parser.add_argument("-g","--nana",dest="nana",type=int, default=300, help="Number of iteration in (g)radient descent")
    parser.add_argument("-m","--remap",dest="remap",type=int, default=5, help="Number of re-(m)apping")
    parser.add_argument("-s","--adjust",dest="adjust",type=int, default=5, help="Number of (s)cale adjusting steps")
    parser.add_argument("-d","--display",dest="disp",action="store_true", default=False, help="Unable interactive display")
    parser.add_argument("-o","--outdir",dest="outdir",type=str, default='.', help="Output directory")

    options = parser.parse_args()
    
    fileA = options.A
    fileB = options.B
    ncell = options.ncell
    fracA = options.fracA
    fracB = options.fracB
    niter = options.niter
    nana = options.nana
    remap = options.remap
    adjust = options.adjust
    disp = options.disp
    outdir = options.outdir
    
    return fileA, fileB, ncell, fracA, fracB, niter, nana, remap, adjust, disp, outdir
    

def find_cell(class_list, positions, tol = 1e-5, frac_tol = 0.5):
    cell = np.zeros((3,3))
    for i in np.unique(class_list):
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
                        if j == 3: break
        else:
            raise RuntimeError("Could not find periodic cell for displacement %d"%i)
        if i==0:
            cell = newcell
        elif abs(la.det(cell)) < abs(la.det(newcell)):
            if np.allclose(la.inv(cell).dot(newcell), np.round(la.inv(cell).dot(newcell)), tol):
                cell = newcell
            else:
                raise RuntimeError("The periodicity of the different classes of displacement is different")
        else:
            if not np.allclose(la.inv(newcell).dot(cell), np.round(la.inv(newcell).dot(cell)),\
     tol):
                raise RuntimeError("The periodicity of the different classes of displacement is different")

    if la.det(cell) < 0:
        cell[:,2] = -cell[:,2]

    return cell

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

def p2ptrans(fileA, fileB, ncell, fracA, fracB, niter, nana, remap, adjust, disp, outdir):

    if not disp:
        matplotlib.use('Agg')

    import matplotlib.pyplot as plt

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    random = False
    
    A = read.poscar(fileA)
    B = read.poscar(fileB)
    
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
    fig.savefig(outdir+'/CellVectors.svg')
    
    
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
    fig.savefig(outdir+'/Agrid.svg')
    
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
    fig.savefig(outdir+'/Bgrid.svg')
    
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
        
    Apos = np.asfortranarray(Apos)
    Bpos = np.asfortranarray(Bpos) 
    t_time = time.time()
    Apos_map, Bpos, Bposst, n_map, class_list, tmat, dmin = tr.fastoptimization(Apos, Bpos, fracA, fracB, Acell, la.inv(Acell), atoms, niter, nana, remap, adjust, 1e-6, 1e-6)
    t_time = time.time() - t_time
    Bpos = np.asanyarray(Bpos)
    Apos = np.asanyarray(Apos)
    
    print("Mapping time:", t_time)
    
    pickle.dump((Apos_map, Bpos, Bposst, n_map, class_list, tmat, dmin), open(outdir+"/fastoptimization.dat","wb"))
    
    # # TMP for testing only -->
    # tr.center(Apos)
    # tr.center(Bpos)
    # Apos_map, Bpos, Bposst, n_map, class_list, tmat, dmin = pickle.load(open("fastoptimization.dat","rb"))
    # # <--
    
    
    print("Total distance between structures:", dmin)
    
    class_list = class_list[:n_map]-1
    
    Bpos = Bpos[:,:n_map]
    Bposst = Bposst[:,:n_map]
    Apos_map = Apos_map[:,:n_map]
    
    natB = np.shape(Bposst)[1] // np.sum(atoms)
    nat = np.shape(Apos)[1] // np.sum(atoms)
    natA = int(fracA*np.shape(Apos)[1]/np.sum(atoms))
    
    
    # Plotting the Apos and Bpos overlayed
    fig = plt.figure(22)
    ax = fig.add_subplot(111, projection='3d')
    ax.view_init(90,0) # TMP
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
    fig.savefig(outdir+'/DispLattice.svg')
    
    # Displacement with stretching
    disps = Apos_map - Bposst
    vec_classes = np.array([np.mean(disps[:,class_list==d_type], axis=1) for d_type in np.unique(class_list)])
    
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.view_init(90,0) # TMP
    maxXAxis = np.max([Apos.max(), Bposst.max()]) + 1
    ax.set_xlim([-maxXAxis, maxXAxis])
    ax.set_ylim([-maxXAxis, maxXAxis])
    ax.set_zlim([-maxXAxis, maxXAxis])
    ax.set_aspect('equal')

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
    
    # Plotting the Apos and Bposst overlayed
    def init_disps():
        #ax.scatter(Apos.T[:,0],Apos.T[:,1])
        for i,num in enumerate(atoms):
            for j in range(num):
                ax.scatter(Apos.T[(np.sum(atoms[:i-1])+j)*nat:(np.sum(atoms[:i-1])+j)*nat + natA,0],Apos.T[(np.sum(atoms[:i-1])+j)*nat:(np.sum(atoms[:i-1])+j)*nat + natA,1],Apos.T[(np.sum(atoms[:i-1])+j)*nat:(np.sum(atoms[:i-1])+j)*nat + natA,2], c="C%d"%(2*i))
                ax.scatter(Apos.T[(np.sum(atoms[:i-1])+j)*nat + natA:(np.sum(atoms[:i-1])+j+1)*nat,0],Apos.T[(np.sum(atoms[:i-1])+j)*nat + natA:(np.sum(atoms[:i-1])+j+1)*nat,1],Apos.T[(np.sum(atoms[:i-1])+j)*nat + natA:(np.sum(atoms[:i-1])+j+1)*nat,2], c="C%d"%(2*i), alpha=0.5)
            ax.scatter(Bposst.T[natB*num*i:natB*num*(i+1),0],Bposst.T[natB*num*i:natB*num*(i+1),1], Bposst.T[natB*num*i:natB*num*(i+1),2], alpha=0.5, c="C%d"%(2*i+1))
        
        ax.quiver(Bposst.T[:,0], Bposst.T[:,1], Bposst.T[:,2], disps.T[:,0], disps.T[:,1], disps.T[:,2])
        fig.savefig(outdir+'/DispLattice_stretched.svg')
        return fig,
    
    init_disps()
    
    anim = animation.FuncAnimation(fig, animate, init_func=init_disps,
                                   frames=490, interval=30)
    
    anim.save(outdir+'/Crystal+Disps.gif', fps=30, codec='gif')
    
    
    # Stretching Matrix
    stMat = la.inv(tr.canonicalize(Bcell)).dot(tr.canonicalize(tmat.dot(Bcell)))
    
    # Rotation Matrix
    rtMat = tmat.dot(Bcell).dot(la.inv(stMat)).dot(la.inv(Bcell))
    
    print("Stretching Matrix:")
    print(stMat)
    
    print("Rotation Matrix:")
    print(rtMat)
    
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
    #     fig.savefig(outdir+'/DispLattice_stretched_%d.svg'%i)
    
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
    fig.savefig(outdir+'/DispLattice_stretched_cell_primittive.svg')
    
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
        fig.savefig(outdir+'/Displacement_structure.svg')
        return fig,
    
    anim = animation.FuncAnimation(fig, animate, init_func=init_struc,
                                   frames=490, interval=30)
    
    anim.save(outdir+'/DispStruc.gif', fps=30, codec='gif')
        
    plt.show()
    
    plt.close('All')

if __name__=='__main__':
    p2ptrans(*readOptions())
    


