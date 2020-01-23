from fmodules import transform as tr
import numpy as np
import numpy.linalg as la
from mpl_toolkits.mplot3d import Axes3D
import matplotlib
from matplotlib import animation
from fmodules import tiling as t
import pickle
import time
from pylada.crystal import Structure, primitive, gruber, read, write, supercell, space_group
from copy import deepcopy
import argparse
import os
import warnings
from format_spglib import from_spglib, to_spglib
from spglib import get_spacegroup

# Default color styles
colorlist=['#929591', 'r', 'k','b','#06470c','#ceb301', '#9e0168', '#26f7fd', '#f97306', '#c20078']
reccolor=['blue','green','red']

# Default extra parameters +++++++
tol = 1e-5
tol_vol = 2*1e-3
tol_uvw = 1e-6
pca = False
nrep = 1
gif = False

# Crystal
ftf = True
ccell1 = np.eye(3)
ccell2 = np.eye(3)
planehkl = [1,1,0]
diruvw = [1,-1,-1]

# Steps
PoscarDirName = "/TransPOSCARS_K-S"
n_steps = 60
viewDirs = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
size = 3
habit = False
n_frames = 5
# ++++++++++++++++++++++++++++++++


try:
    from config import *
except ModuleNotFoundError:
    pass

def readOptions():

    parser = argparse.ArgumentParser()
    parser.add_argument("-I","--initial",dest="A",type=str, default='./POSCAR_A', help="Initial Structure")
    parser.add_argument("-F","--final",dest="B",type=str, default='./POSCAR_B', help="Final Structure")
    parser.add_argument("-n","--ncell",dest="ncell",type=int, default=300, help="Number of cells to tile")
    parser.add_argument("-i","--interactive",dest="interactive",action="store_true", default=False, help="Enables interactive display")
    parser.add_argument("-d","--disp",dest="savedisplay",action="store_true", default=False, help="Saves figures")
    parser.add_argument("-p","--param", dest="filename", type=str, default='./p2p.in', help="Parameter file")
    parser.add_argument("-o","--outdir",dest="outdir", type=str, default='.', help="Output directory")
    parser.add_argument("-u","--use",dest="use", type=str, default=None, help="Use previously calculated data")
    parser.add_argument("-m","--minimize",dest="minimize",action="store_true", default=False, help="Rerun minimization using the parameters stored in the folder provided with the -u option")
    parser.add_argument("-s","--switch",dest="switch",action="store_true", default=False, help="Map the larger cell on the smaller cell")
    parser.add_argument("-r","--noprim",dest="prim",action="store_false", default=True, help="Finds the primitive cell at the beginning") #TMP
    parser.add_argument("-a","--anim",dest="anim",action="store_true", default=False, help="Produce the animation") #TMP
    parser.add_argument("-v","--vol",dest="vol",action="store_true", default=False, help="Make the two (stochiometric) cells equal in volume")
    

    options = parser.parse_args()
    
    fileA = options.A
    fileB = options.B
    ncell = options.ncell
    filename = options.filename
    savedisplay = options.savedisplay
    interactive = options.interactive
    if options.use == None:
        use = False
        minimize = True
        outdir = options.outdir
    else:
        use = True
        minimize = options.minimize
        outdir = options.use
    switch = options.switch
    prim = options.prim
    anim = options.anim
    vol = options.vol
    
    return fileA, fileB, ncell, filename, interactive, savedisplay, outdir, use, switch, prim, anim, vol, minimize

def normal(A):
    return A/np.ones((3,1)).dot(la.norm(A, axis=0).reshape((1,np.shape(A)[1])))

def find_uvw(stretch_dir, basis = np.eye(3)):
    min_dist = np.ones(3)*1000
    min_uvw = np.zeros((3,np.shape(stretch_dir)[1]))
    stretch_dir = normal(stretch_dir)
    for l,vec in enumerate(stretch_dir.T):
        for i in np.arange(-10,10):
            for j in np.arange(-10,10):
                for k in np.arange(-10,10):
                    if [i,j,k] != [0,0,0]:
                        unit_vec = basis.dot(np.array([i,j,k]))
                        unit_vec = unit_vec/la.norm(unit_vec)
                        cur_dist = la.norm(unit_vec-vec)
                        if cur_dist < min_dist[l]:
                            min_dist[l] = cur_dist
                            min_uvw[:,l] = np.array([i,j,k], np.int).T

        gcd3 = gcd(min_uvw[0,l],gcd(min_uvw[1,l], min_uvw[2,l]))
        min_uvw[:,l] = min_uvw[:,l]/gcd3

    return min_uvw

def find_supercell(cell, newcell, tol):
    
    for i in range(1,10):
        for j in range(1,10):
            for k in range(1,10):
                if abs(la.det(cell)) < abs(la.det(newcell)):
                    if np.allclose(la.inv(cell).dot(newcell).dot(np.diag([i,j,k])), np.round(la.inv(cell).dot(newcell).dot(np.diag([i,j,k]))), tol):
                        newcell = newcell.dot(np.diag([i,j,k]))
                        break
                else:
                    if np.allclose(la.inv(newcell).dot(cell).dot(np.diag([i,j,k])), np.round(la.inv(newcell).dot(cell).dot(np.diag([i,j,k]))), tol):
                        cell = cell.dot(np.diag([i,j,k]))
                        break
            else:
                continue
            break
        else:
            continue
        break
    else:
        return cell, None

    return cell, newcell
                
def find_multiples(vec, pos):
    """Goes through the displacements and finds the one that are parallel"""
    multiple = [0]
    for p in pos.T:
        if la.norm(np.cross(p,vec))/(la.norm(p) * la.norm(vec)) < tol:
            multiple.append(p.dot(vec)/la.norm(vec)**2)
    return multiple

def find_cell(class_list, positions, tol = 1e-5, frac_shell = 0.5, frac_correct = 0.95, max_count=1000):
    
    cm = np.mean(positions, axis=1).reshape((3,1))

    for loop in [0,1]:
        for i in np.unique(class_list):
            pos = positions[:, class_list == i]
            center = np.argmin(la.norm(pos, axis = 0))
            list_in = list(range(np.shape(pos)[1]))
            list_in.remove(center)
            origin = pos[:,center:center+1]
            pos = pos[:,list_in] - origin.dot(np.ones((1,np.shape(pos)[1]-1))) # centered
            if not loop:
                norms = la.norm(pos, axis = 0)
                idx = np.argsort(norms)
                minj = 0
                maxj = len(idx)
            else:
                idx = np.arange(np.shape(pos)[1])
                np.random.shuffle(idx)
                minj = 3
                maxj = len(idx)
            count = 0

            for j in range(minj, maxj):
                if not loop:
                    mink = j+1
                    maxk = len(idx)
                else:
                    mink = 0
                    maxk = j-1
                for k in range(mink, maxk):
                    if not loop:
                        minl = k+1
                        maxl = len(idx)
                    else:
                        minl = 0
                        maxl = k-1
                    for l in range(minl, maxl):
                        # creates all possible cells 
                        newcell=np.concatenate([pos[:,idx[j]:idx[j]+1], 
                                                    pos[:,idx[k]:idx[k]+1], 
                                                    pos[:,idx[l]:idx[l]+1]],axis=1)

                        if abs(la.det(newcell)) > tol: # Cell as non-zero volume
                            count += 1
                            if count > max_count:
                                break

                            if la.det(newcell) < 0:
                                newcell=np.concatenate([pos[:,idx[j]:idx[j]+1],
                                                        pos[:,idx[l]:idx[l]+1],
                                                        pos[:,idx[k]:idx[k]+1]],axis=1)

                            norms = la.norm(positions - cm.dot(np.ones((1,np.shape(positions)[1]))), axis=0)
                            apos = la.inv(newcell).dot(positions - origin.dot(np.ones((1,np.shape(positions)[1]))))

                            inPos = apos[:,np.sum((apos < 1 - tol) & (apos > - tol),0)==3]
                            inType = class_list[np.sum((apos < 1 - tol) & (apos > - tol),0)==3]
                            
                            n_map = 0
                            for m, a in enumerate(apos.T):
                                for n, b in enumerate(inPos.T):
                                    # Check that the cell is repeating
                                    if (all(abs(np.mod(a+tol,1)-tol-b) < tol) and 
                                        inType[n] == class_list[m]):
                                        break
                                else:
                                    continue
                                n_map += 1
        
                            genPos = []
                            if float(n_map)/float(len(class_list)) > frac_correct:
                                xMax = int(np.max(apos[0,:]))+1
                                xMin = int(np.min(apos[0,:]))-1
                                yMax = int(np.max(apos[1,:]))+1
                                yMin = int(np.min(apos[1,:]))-1
                                zMax = int(np.max(apos[2,:]))+1
                                zMin = int(np.min(apos[2,:]))-1
                                for x in range(xMin,xMax):
                                    for y in range(yMin,yMax):
                                        for z in range(zMin,zMax):
                                            genPos.append(inPos + np.array([[x,y,z]]).T.dot(np.ones((1,np.shape(inPos)[1]))))
                            
                                genPos = newcell.dot(np.concatenate(genPos,axis=1))
                                genPos = genPos + origin.dot(np.ones((1,np.shape(genPos)[1])))
        
                                if np.sum(la.norm(genPos - cm.dot(np.ones((1,np.shape(genPos)[1]))), axis = 0) < frac_shell * np.max(norms) - tol) == np.sum(norms < frac_shell * np.max(norms) - tol):
                                    return newcell, origin
                    else:
                        continue
                    break
                else:
                    continue
                break
            print("WARNING: Could not find periodic cell using displacement %d. Increase sample size or use results with care."%i)
        print("WARNING: Could not find cell using shortest distances, trying random order") 
    print("WARNING: Could not find periodic cell for any displacement. Increase sample size.")                
    return None, None
    
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

def uniqueclose(closest, tol):
    unique = []
    idx = []
    for i,line in enumerate(closest.T):
        there = False
        for j,check in enumerate(unique):
            if np.allclose(check, line, atol=tol):
                there = True
                idx[j].append(i) 
        if not there:
            unique.append(line)
            idx.append([i])
    return (np.array(idx), np.array(unique))

def dir2angles(plane):
    plane = plane/la.norm(plane)
    angles=np.zeros(2)
    a0 = np.arccos(plane[2])
    angles[0] = np.pi/2 - a0
    angles[1] = np.arctan2(plane[1], plane[0])
    return angles*180/np.pi

def rotate(icell,fcell):
    U,S,V = la.svd(icell.dot(fcell.T))
    return V.conj().T.dot(U.conj().T).real

def rot_mat(u, theta):
    u = u/la.norm(u)

    P = u.reshape((3,1)).dot(u.reshape((1,3)))
    Q = np.array([[0,-u[2],u[1]], [u[2], 0, -u[0]], [-u[1], u[0], 0]])

    return  P + (np.eye(3) - P)*np.cos(theta) + Q*np.sin(theta)

def set_view(p,angle=0):
    p=p / la.norm(p)
    v1 = np.array([0,p[2],-p[1]])
    v1 = v1 / la.norm(v1)
    v2 = np.cross(p,v1)
    v2 = v2 / la.norm(v2)
    return la.inv(np.array([v1*np.cos(angle) + v2*np.sin(angle), -v1*np.sin(angle) + v2*np.cos(angle),
                     p]).T)

def makeSphere(A, ncell, *atom_types):

    """ 
    Creates a spherical set of atoms.
  
    The center of the set of point is defined by the origin of structures A.

    Parameters: 
    A (Structure): Crystal Structure from which to create the set of points  

    ncell (int): Number of A unit cells in the set of point

    atom_types (np.array): (optional) if provided indicates the atom type (specie) for each section of the output

    Returns: 
    np.array: 3*ncell of real number representing atomic positions. Divided into sections according to atom_types. In each section atoms are ordered from the closest to the furthest from the center of the sphere. 
    np.array: 1*natoms array indicates the number of ncell slots for each type of atoms in each section of the first output.
    np.array: (optional) if atom_types was no provided as an input it is provided as an output

    """

    Acell = A.cell*float(A.scale)
    
    if atom_types == ():

        # Adds atoms to A and B (for cell with different types of atoms)
        Apos = []
        atom_types = np.array([], np.str)
        atomsA = np.array([], np.int)
        for a in A:
            if any(atom_types == a.type):
                idx = np.where(atom_types == a.type)[0][0]
                Apos[idx] = np.concatenate((Apos[idx], t.sphere(Acell, ncell, a.pos*float(A.scale))), axis = 1) 

                # Order the atoms with respect to distance
                Apos[idx] = Apos[idx][:,np.argsort(la.norm(Apos[idx],axis=0))] 
                atomsA[idx] += 1
            else:
                Apos.append(t.sphere(Acell, ncell, a.pos*float(A.scale)))
                atom_types = np.append(atom_types, a.type)
                atomsA = np.append(atomsA,1)
        
        return (np.concatenate(Apos, axis=1), atomsA, atom_types)

    elif len(atom_types) == 1 and type(atom_types[0]) == np.ndarray:

        atom_types = atom_types[0]
        Apos = [None]*len(atom_types)
        atomsA = np.zeros(len(atom_types), np.int)
        for a in A:
            idx = np.where(atom_types == a.type)[0][0]
            if atomsA[idx] == 0:
                Apos[idx] = t.sphere(Acell, ncell, a.pos*float(A.scale))
            else:
                Apos[idx] = np.concatenate((Apos[idx], t.sphere(Acell, ncell, a.pos*float(A.scale))), axis = 1) 
                # Order the atoms with respect to distance
                Apos[idx] = Apos[idx][:,np.argsort(la.norm(Apos[idx],axis=0))]
            atomsA[idx] += 1
        
        return (np.concatenate(Apos, axis=1), atomsA)

    else:
        raise ValueError("The only optional argument must be atom_types")
    
def atCenter(pos):
    """ Align the center of mass of any 3*n array to the origin """
    n = np.shape(pos)[1]
    return  pos - (np.sum(pos,axis=1).reshape((3,1))/n).dot(np.ones((1,n)))

def makeStructures(cell, atoms, atom_types, natB, pos_in_struc, class_list, vec_classes):
    """Make the displacement structure from the repeating unit cell"""
    
    # Small function to determine type
    def whattype(pos, nat):

        pos = pos//nat + 1

        atom_tot = np.sum(np.triu(atoms.reshape((len(atoms),1)).dot(np.ones((1,len(atoms))))), axis=0)
        
        return atom_types[np.nonzero(atom_tot >= pos)[0][0]]
    
    # Make a pylada structure
    cell_coord = np.mod(la.inv(cell).dot(pos_in_struc)+tol,1)-tol
    dispStruc = Structure(cell)
    stinitStruc = Structure(cell)
    incell = []

    for idx, disp in zip(*uniqueclose(cell_coord, tol)):
        for i in idx:
            if np.allclose(pos_in_struc[:,i], cell.dot(disp), atol=tol):
                incell.append((i,pos_in_struc[:,i]))
                break
        else:
            i = np.argmin(la.norm(np.array([pos_in_struc[:,j] for j in idx]),axis=1))
            incell.append((i,cell.dot(disp)))
        
    for i, disp in incell:
        dispStruc.add_atom(*(tuple(disp)+(str(class_list[i]),)))
        stinitStruc.add_atom(*(tuple(disp)+(whattype(i, natB),)))
        
    if la.det(cell) < 0:
       cell[:,2] = -cell[:,2] 

    # Finds a squarer cell
    cell = gruber(cell)

    dispStruc = supercell(dispStruc, cell)

    # Makes sure it is the primitive cell 
    dispStruc = primitive(dispStruc, tolerance = tol)

    tmpStruc = Structure(dispStruc.cell)
    to_add = [np.mod(la.inv(dispStruc.cell).dot(a.pos)+tol,1)-tol for a in stinitStruc]
    for idx, pos in zip(*uniqueclose(np.array(to_add).T, tol)):
        tmpStruc.add_atom(*dispStruc.cell.dot(pos),stinitStruc[idx[0]].type)
    
    stinitStruc = tmpStruc

    write.poscar(stinitStruc, vasp5=True, file="POSCAR_init")

    finalStruc = Structure(dispStruc.cell)
    for i,a in enumerate(dispStruc):
        finalStruc.add_atom(*(a.pos+vec_classes[int(a.type)]),stinitStruc[i].type)

    return dispStruc, stinitStruc, finalStruc
    
def displayOptimalResult(Apos, Bpos, Bposst, disps_total, disps, class_list, vec_classes,
                         nat, natA, natB, atoms, outdir, savedisplay, interactive):

    """Displays or saves pictures of the optimal result"""
    
    fig = plt.figure("Optimal Result", figsize=(15,5))

    # Left-most pannel: All atoms and connections only rigid rotation
    ax = fig.add_subplot(131, projection='3d')
    ax.set_title('Optimal result (rigid rotation)')
    maxXAxis = np.max([Apos.max(), Bpos.max()]) + 1
    ax.set_xlim([-maxXAxis, maxXAxis])
    ax.set_ylim([-maxXAxis, maxXAxis])
    ax.set_zlim([-maxXAxis, maxXAxis])
    num_tot = 0

    for i,num in enumerate(atoms):
        ax.scatter(Apos.T[num_tot*nat:num_tot*nat+natA*num+1,0],Apos.T[num_tot*nat:num_tot*nat+natA*num+1,1],Apos.T[num_tot*nat:num_tot*nat+natA*num+1,2], c=colorlist[2*i])
        ax.scatter(Apos.T[num_tot*nat+natA*num:(num_tot + num)*nat+1,0],Apos.T[num_tot*nat+natA*num:(num_tot + num)*nat+1,1],Apos.T[num_tot*nat+natA*num:(num_tot + num)*nat+1,2], c=colorlist[2*i], alpha=0.1)
        ax.scatter(Bpos.T[natB*num_tot:natB*(num_tot+num),0],Bpos.T[natB*num_tot:natB*(num_tot+num),1], Bpos.T[natB*num_tot:natB*(num_tot+num),2], c=colorlist[2*i+1])
        num_tot = num_tot + num

    # # This could be useful in certain cases to know how much th structures were displaced
    # centerofmassA = np.mean(Apos,axis=1)
    # centerofmassB = np.mean(Bpos,axis=1)

    # ax.scatter(centerofmassA[0], centerofmassA[1], centerofmassA[2], s=60, c='red')
    # ax.scatter(centerofmassB[0], centerofmassB[1], centerofmassB[2], s=60, c='green')
    
    for i in range(len(vec_classes)):
        disps_class = disps_total[:,class_list==i]
        Bpos_class = Bpos[:,class_list==i]
        ndisps = np.shape(disps_class)[1]
        ax.quiver(Bpos_class.T[:,0], Bpos_class.T[:,1], Bpos_class.T[:,2], disps_class.T[:,0], disps_class.T[:,1], disps_class.T[:,2], color=colorlist[i%10])

    # Middle-panel
    ax = fig.add_subplot(132, projection='3d')
    ax.set_title('Optimal result (stretched)')
    maxXAxis = np.max([Apos.max(), Bposst.max()]) + 1
    ax.set_xlim([-maxXAxis, maxXAxis])
    ax.set_ylim([-maxXAxis, maxXAxis])
    ax.set_zlim([-maxXAxis, maxXAxis])
    num_tot = 0
    for i,num in enumerate(atoms):
        ax.scatter(Apos.T[num_tot*nat:num_tot*nat+natA*num+1,0],Apos.T[num_tot*nat:num_tot*nat+natA*num+1,1],Apos.T[num_tot*nat:num_tot*nat+natA*num+1,2], c=colorlist[2*i])
        ax.scatter(Apos.T[num_tot*nat+natA*num:(num_tot + num)*nat+1,0],Apos.T[num_tot*nat+natA*num:(num_tot + num)*nat+1,1],Apos.T[num_tot*nat+natA*num:(num_tot + num)*nat+1,2], c=colorlist[2*i], alpha=0.1)
        ax.scatter(Bposst.T[natB*num_tot:natB*(num_tot+num),0],Bposst.T[natB*num_tot:natB*(num_tot+num),1], Bposst.T[natB*num_tot:natB*(num_tot+num),2], c=colorlist[2*i+1])
        num_tot = num_tot + num

    for i in range(len(vec_classes)):
        disps_class = disps[:,class_list==i]
        Bposst_class = Bposst[:,class_list==i]
        ndisps = np.shape(disps_class)[1]
        ax.quiver(Bposst_class.T[:,0], Bposst_class.T[:,1], Bposst_class.T[:,2], disps_class.T[:,0], disps_class.T[:,1], disps_class.T[:,2], color=colorlist[i%10])

    # Right-panel
    ax = fig.add_subplot(133, projection='3d')
    ax.set_title('Displacement Classes')
    maxXAxis = disps.max()
    ax.set_xlim([-maxXAxis, maxXAxis])
    ax.set_ylim([-maxXAxis, maxXAxis])
    ax.set_zlim([-maxXAxis, maxXAxis])
    
    for i in range(len(vec_classes)):
        disps_class = disps[:,class_list==i]
        ndisps = np.shape(disps_class)[1]
        ax.quiver(np.zeros((1,ndisps)), np.zeros((1,ndisps)), np.zeros((1,ndisps)), disps_class.T[:,0], disps_class.T[:,1], disps_class.T[:,2], color=colorlist[i%10])

    if savedisplay:
        fig.savefig(outdir+'/optimal_result.svg')
        print("Saved display in %s"%(outdir+'/optimalRes.svg'))
        
    if interactive:
        plt.show()

def makeGif(Apos, Bposst, disps, vec_classes, nat, atoms):
    fig = plt.figure("gif")
    ax = fig.add_subplot(111, projection='3d')
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
        num_tot = 0
        for i,num in enumerate(atoms):
            ax.scatter(Apos.T[num_tot*nat:num_tot*nat+natA*num+1,0],Apos.T[num_tot*nat:num_tot*nat+natA*num+1,1],Apos.T[num_tot*nat:num_tot*nat+natA*num+1,2], c=colorlist[2*i])
            ax.scatter(Apos.T[num_tot*nat+natA*num:(num_tot + num)*nat+1,0],Apos.T[num_tot*nat+natA*num:(num_tot + num)*nat+1,1],Apos.T[num_tot*nat+natA*num:(num_tot + num)*nat+1,2], c=colorlist[2*i], alpha=0.1)
            ax.scatter(Bposst.T[natB*num_tot:natB*(num_tot+num),0],Bposst.T[natB*num_tot:natB*(num_tot+num),1], Bposst.T[natB*num_tot:natB*(num_tot+num),2], c=colorlist[2*i+1])
            num_tot = num_tot + num
    
        for i in range(len(vec_classes)):
            disps_class = disps[:,class_list==i]
            Bposst_class = Bposst[:,class_list==i]
            ndisps = np.shape(disps_class)[1]
            ax.quiver(Bposst_class.T[:,0], Bposst_class.T[:,1], Bposst_class.T[:,2], disps_class.T[:,0], disps_class.T[:,1], disps_class.T[:,2], color=colorlist[i%10])
        return fig,
    
    anim = animation.FuncAnimation(fig, animate, init_func=init_disps,
                                   frames=490, interval=30)
    anim.save(outdir+'/Crystal+Disps.gif', fps=30, codec='gif')

def displayTransCell(disps, dispStruc, finalStruc, foundcell,
                     pos_in_struc, vec_classes, interactive, savedisplay):   

    cell = dispStruc.cell
    
    # Displays only the cell and the displacements in it
    fig = plt.figure("Transformation Cell", figsize = [10,5])
    
    ax = fig.add_subplot(121, projection='3d')
    ax.set_title('Transformation cell with displacements\n and atomic positions')
    for i,disp in enumerate(dispStruc):
        ax.quiver(disp.pos[0], disp.pos[1], disp.pos[2], vec_classes[int(disp.type)][0],vec_classes[int(disp.type)][1], vec_classes[int(disp.type)][2], color=colorlist[i%10])
        ax.scatter(disp.pos[0], disp.pos[1], disp.pos[2], alpha = 0.5, s=10, color=colorlist[i%10])
        ax.scatter(finalStruc[i].pos[0], finalStruc[i].pos[1], finalStruc[i].pos[2], alpha = 1, s=10, color=colorlist[i%10])
    ax.quiver(np.zeros(3), np.zeros(3), np.zeros(3), cell[0,:], cell[1,:], cell[2,:], color = "red", alpha = 0.3)
    # ax.quiver(np.zeros(3), np.zeros(3), np.zeros(3), Acell[0,:], Acell[1,:], Acell[2,:], color = "blue", alpha = 0.3)
    maxXAxis = abs(cell).max() + 1
    ax.set_xlim([-maxXAxis, maxXAxis])
    ax.set_ylim([-maxXAxis, maxXAxis])
    ax.set_zlim([-maxXAxis, maxXAxis])

    # Displays displacement with the disp cell overlayed
    ax = fig.add_subplot(122, projection='3d')
    ax.set_title('Transformation cell in displacement crystal \n as found (red), primitive (green)')
    ax.quiver(pos_in_struc.T[:,0], pos_in_struc.T[:,1], pos_in_struc.T[:,2], disps.T[:,0], disps.T[:,1], disps.T[:,2], color = "C0")
    ax.scatter(pos_in_struc.T[:,0], pos_in_struc.T[:,1], pos_in_struc.T[:,2], s=10, color = "C0")
    ax.quiver(np.zeros(3), np.zeros(3), np.zeros(3), foundcell[0,:], foundcell[1,:], foundcell[2,:], color = "red")
    ax.quiver(np.zeros(3), np.zeros(3), np.zeros(3), cell[0,:], cell[1,:], cell[2,:], color = "green")
    maxXAxis = pos_in_struc.max() + 1
    ax.set_xlim([-maxXAxis, maxXAxis])
    ax.set_ylim([-maxXAxis, maxXAxis])
    ax.set_zlim([-maxXAxis, maxXAxis])
    
    if savedisplay:
        fig.savefig(outdir+'/transformation_cell.svg')

    if interactive:
        print("(Close the display window to continue)")
        plt.show()
    
def printMatAndDir(A,ccell):
    print("--------Matrix--------|-----Closest uvw------")
    print("    v1    v2    v3    |    d1    d2    d3    ")
    B = find_uvw(A, basis=ccell)
    A = la.inv(ccell).dot(A)
    for i, rowA in enumerate(A):
          rowB = B[i]
          print(' '.join(["% 5.3f"%(val) for val in rowA]), " |", ' '.join(["%5d"%(val) for val in rowB]))

def strainDirs(tmat):

    if ftf:
        eigval, P = la.eig(tmat.T.dot(tmat))
        eigval = np.sqrt(eigval)

        idx = np.argsort(eigval)
        eigval = eigval[idx]
        P = P[:,idx]
        
        U = P.dot(np.diag(eigval)).dot(P.T)

        invEigval, Q = la.eig(la.inv(tmat).T.dot(la.inv(tmat)))
        invEigval = np.sqrt(invEigval)
        idx = np.argsort(1/invEigval)
        Q = Q[:,idx]
        
    else:
        U = rotate(tmat, np.eye(3)).dot(tmat)

        eigval, P = la.eig(U)

        idx = np.argsort(eigval)
        eigval = eigval[idx]
        P = P[:,idx]
        
        iU = rotate(la.inv(tmat), np.eye(3)).dot(la.inv(tmat))
        
        invEigval, Q = la.eig(iU)
        invEigval = np.sqrt(invEigval)
        idx = np.argsort(1/invEigval)
        Q = Q[:,idx]
        
    P = normal(P)
    
    Q = normal(Q)

    return eigval, U, P, Q

def findHabit(U, P, eigval):
    
    # Using uniformly strained plane
    ratio = np.sqrt(abs((eigval[2]**2 - eigval[1]**2)/(eigval[1]**2 - eigval[0]**2)))
    
    # uvw direction of normal in cartesian coord
    planeHab = np.zeros((3,2))
    planeHab[:,0] = P[:,0] + ratio*P[:,2]
    planeHab[:,1] = P[:,0] - ratio*P[:,2]
    return planeHab, ratio

def findR(U, P=None, planeHab=None, ratio=None):

    if planeHab is None or ratio is None or P is None:
        eigval, P = la.eig(mat)
        planeHab, ratio = findHabit(U, P, eigval)
    
    V = np.zeros((2,3,3))
    M = np.zeros((2,3,3))
    R = np.zeros((2,3,3))
    
    for i in range(2):
        V[i,:,0] = ratio*P[:,0] + (1-2*i)*P[:,2]
        V[i,:,1] = P[:,1]
        V[i,:,2] = planeHab[:,i]/la.norm(planeHab[:,i]) 
    
        M[i,:,:2] = U.dot(V[i,:,:2])
        M[i,:,2] = np.cross(M[i,:,0], M[i,:,1])
        M[i,:,2] = M[i,:,2] / la.norm(M[i,:,2])

        R[i,:,:] = M[i,:,:].dot(M[i,:,:].T)
        
    return R

def crystallography(tmat, A, B, fileA, fileB, ccellA, ccellB):
    print("----------CRYSTALLOGRAPHY----------")
    print()
    eigval, U, P, Q = strainDirs(tmat)

    print("Strain Directions in %s (%s) coordinates:"%(B.name, fileB))
    print("    d1    d2    d3    ")
    printMatAndDir(P, ccellA)
    print()
    print("Strain Directions in %s (%s) coordinates:"%(A.name, fileA))
    print("    d1    d2    d3    ")
    printMatAndDir(Q, ccellB)
    print()
    print("Strains + 1 (eigenvalues)")
    print("    e1    e2    e3    ")
    print(' '.join(["% 5.3f"%(val) for val in eigval])) 
    print()

    rcellA = la.inv(ccellA) #Reciprocal cell
    rcellB = la.inv(ccellB) #Reciprocal cell

    planeHab, ratio = findHabit(U, P, eigval)
        
    print("Uniformly strained planes:")
    print()
    print("Exact plane in %s (%s) coordinates:"%(B.name, fileB))
    print("(+): (% 6.4f, % 6.4f, % 6.4f)"%(*ccellB.dot(planeHab[:,0]),))
    print("(-): (% 6.4f, % 6.4f, % 6.4f)"%(*ccellB.dot(planeHab[:,1]),))
    print()
    print("Closest uvw:")
    print("(+): (% d, % d, % d)"%(*find_uvw(planeHab[:,1:2], rcellB),))
    print("(-): (% d, % d, % d)"%(*find_uvw(planeHab[:,0:1], rcellB),))
    print()

    R = findR(U, P=P, planeHab=planeHab, ratio=ratio)
    
    print("Orientation Relationship with habit plane:")
    print()
    orMat = np.zeros((2,3,3))
    resPlanehkl = np.zeros((2,3))
    resDiruvw = np.zeros((2,3))
    for i in range(2):
        orMat[i,:,:] = Q.dot(P.T).dot(R[i,:,:].T)
        resPlanehkl[i,:] = ccellB.dot(orMat[i,:,:].dot(rcellA.dot(planehkl)))
        resDiruvw[i,:] = la.inv(ccellB).dot(orMat[i,:,:].dot(ccellA.dot(diruvw)))
    print("%s (%s) // %s (%s)"%(B.name, fileB, A.name, fileA))
    print("(+): (% 2d, % 2d, % 2d) [% 2d, % 2d, % 2d] // (% 6.4f, % 6.4f, % 6.4f) [% 6.4f, % 6.4f, % 6.4f]"%(*planehkl, *diruvw, *resPlanehkl[0,:], *resDiruvw[0,:]))
    print("(-): (% 2d, % 2d, % 2d) [% 2d, % 2d, % 2d] // (% 6.4f, % 6.4f, % 6.4f) [% 6.4f, % 6.4f, % 6.4f]"%(*planehkl, *diruvw, *resPlanehkl[1,:], *resDiruvw[1,:]))
    print()
    print("Approximate low index OR")
    resPlanehklClose = np.zeros((2,3))
    resDiruvwClose = np.zeros((2,3))
    for i in range(2):
        resPlanehklClose[i,:] = find_uvw(rcellB.dot(resPlanehkl[i,:].reshape((3,1))), rcellB).T[0]
        resDiruvwClose[i,:] = find_uvw(ccellB.dot(resDiruvw[i,:].reshape((3,1))), ccellB).T[0]
    print("(+): (% 2d, % 2d, % 2d) [% 2d, % 2d, % 2d] //  (% d, % d, % d) [% d, % d, % d]"%(*planehkl, *diruvw, *resPlanehklClose[0,:], *resDiruvwClose[0,:]))
    print("(-): (% 2d, % 2d, % 2d) [% 2d, % 2d, % 2d] //  (% d, % d, % d) [% d, % d, % d]"%(*planehkl, *diruvw, *resPlanehklClose[1,:], *resDiruvwClose[1,:]))
    print()
    
    print("Orientation Relationship in thin films:")
    print()
    orMat = Q.dot(P.T)
    resPlanehkl = ccellB.dot(orMat.dot(rcellA.dot(planehkl)))
    resDiruvw = la.inv(ccellB).dot(orMat.dot(ccellA.dot(diruvw)))
    print("%s (%s) // %s (%s)"%(B.name, fileB, A.name, fileA))
    print("(% 2d, % 2d, % 2d) [% 2d, % 2d, % 2d] // (% 6.4f, % 6.4f, % 6.4f) [% 6.4f, % 6.4f, % 6.4f]"%(*planehkl, *diruvw, *resPlanehkl, *resDiruvw))
    print()
    print("Approximate low index OR")
    resPlanehklClose = find_uvw(rcellB.dot(resPlanehkl.reshape((3,1))), rcellB).T[0]
    resDiruvwClose = find_uvw(ccellB.dot(resDiruvw.reshape((3,1))), ccellB).T[0]
    print("(% 2d, % 2d, % 2d) [% 2d, % 2d, % 2d] //  (% d, % d, % d) [% d, % d, % d]"%(*planehkl, *diruvw, *resPlanehklClose, *resDiruvwClose))
    print()

    return eigval, U, P, Q, planeHab 

def produceTransition(tmat, dispStruc, finalStruc, vec_classes, outdir, atom_types,
                      anim, savedisplay, interactive):

    os.makedirs(outdir+PoscarDirName, exist_ok=True)

    itmat = rotate(la.inv(tmat).dot(finalStruc.cell), finalStruc.cell).dot(la.inv(tmat))

    spgList = []
    transStruc = []
    color_array = []
    Tpos = [] 
    for i in range(n_steps+1):
        if habit:
            curMat = find_R_RU(curMat).dot(curMat)
        else:
            curMat = (itmat-np.eye(3))*i/n_steps + np.eye(3)
        curStruc = Structure(curMat.dot(finalStruc.cell))
        for j,a in enumerate(dispStruc):
            curDisp = vec_classes[int(a.type)]*i/n_steps
            curPos = curMat.dot((finalStruc[j].pos - curDisp).reshape((3,1)))
            curStruc.add_atom(*(curPos.T.tolist()[0]),finalStruc[j].type)
            # curStruc = supercell(curStruc, curStruc.cell)
        write.poscar(curStruc, vasp5=True, file=outdir+PoscarDirName+"/POSCAR_%03d"%i) # Write resulting structure in POSCAR
        spgList.append(get_spacegroup(to_spglib(curStruc), symprec=0.3, angle_tolerance=3.0))
        transStruc.append(curStruc)

        if savedisplay or anim or interactive:
            color_array.append([])
            Tpos.append([])
            for l,pl in enumerate(viewDirs):
            
                if isinstance(pl[0],(list, np.ndarray)): 
                    plane = normal(curStruc.cell).dot(pl[0] + i/n_steps*(pl[1] - pl[0]))
                else:
                    plane = normal(curStruc.cell).dot(pl)
                tickness = 5 * la.norm(plane)
                plane = plane/tickness
                
                lattices = []
            
                Tpos_tmp = []
                types = []
            
                sizes = np.round(max(la.norm(curStruc.cell,axis=0))*size/la.norm(curStruc.cell,axis=0))
            
                # Cut along the current plane
                for l in np.arange(-sizes[0],sizes[0]):
                    for j in np.arange(-sizes[1],sizes[1]):
                        for k in np.arange(-sizes[2],sizes[2]):
                            for a in curStruc:
                                pos = curStruc.cell.dot(np.array([l,j,k])) + a.pos
                                if plane.dot(pos) < 0 and plane.dot(pos) > - tickness:
                                    Tpos_tmp.append(pos)
                                    types.append(np.where(atom_types == a.type)[0][0])
                                    
                Tpos[-1].append(np.array(Tpos_tmp).T)
                color_array[-1].append(types)
                
    return transStruc, spgList, Tpos, color_array

def add_panel(fig,g,plane, state, p, anchor, Tpos, color_array, transStruc):
    Rectangle = matplotlib.patches.Rectangle
    
    ax = fig.add_subplot(g, projection='3d', proj_type = 'ortho')
    ax.set_anchor(anchor)
    # ax.view_init(*angles)
    ax.view_init(azim=-90, elev=90) # x-y plane view
    maxXAxis = np.abs([c for a in Tpos for b in a for c in b.flatten()]).max() + 1
    ax.set_xlim([-maxXAxis, maxXAxis])
    ax.set_ylim([-maxXAxis, maxXAxis])
    ax.set_zlim([-maxXAxis, maxXAxis])
    
    toplot = set_view(-plane).dot(Tpos[state][p])
    color_to_plot = np.array(color_array[state][p])
    idxx = np.argsort(toplot[2,:])
    toplot = toplot[:,idxx]
    color_to_plot = color_to_plot[idxx]
    ax.set_axis_off()
    ax.dist = 2
    axlims = [a for b in ax.get_position().get_points() for a in b]
    rec = Rectangle((axlims[0],axlims[1]),(axlims[2]-axlims[0]),(axlims[3]-axlims[1]), transform = fig.transFigure, fill=False,lw=3, color=reccolor[p])
    fig.patches.append(rec)
    for i,point in enumerate(toplot.T):
        ax.scatter(*point, c=colorlist[color_to_plot[i]], s=40, depthshade=False, alpha = float(i)/len(toplot.T))

    ax2d = fig.add_subplot(g)
    ax2d.set_axis_off()
    ax2d.patch.set_alpha(0)
    ax2d.set_anchor(anchor)
    ax2d.set_xlim([-maxXAxis, maxXAxis])
    ax2d.set_ylim([-maxXAxis, maxXAxis])
    ax2d.set_aspect('equal')
    
    for a,ar in enumerate(viewDirs):
    # tmpdir = np.zeros((3,3))
    # tmpdir[0,0] = 0
    # tmpdir[0,2] = ratio
    # tmpdir[0,1] = -1
    # tmpdir[1,0] = 0
    # tmpdir[1,2] = 1
    # tmpdir[1,1] = ratio
    # tmpdir[2,0] = 1
    # tmpdir[2,2] = 0
    # tmpdir[2,1] = 0
    # for a,ar in enumerate(tmpdir):
        if a!=p:
            if isinstance(ar[0],(list, np.ndarray)): 
                arrow = transStruc[state].cell.dot(ar[0] + state/n_steps*(ar[1] - ar[0]))
                # arrow = normal(transStruc[state].cell).dot(ar[0] + state/n_steps*(ar[1] - ar[0]))
            else:
                arrow = transStruc[state].cell.dot(ar)
                # arrow = normal(transStruc[state].cell).dot(ar)
            arrow = set_view(-plane).dot(arrow)
            ax2d.quiver(*np.zeros(2), *arrow[:2], scale_units='inches', scale=10, color=reccolor[a])
            #ax.quiver(*np.zeros(3), *arrow, color=reccolor[a])
    #ax.quiver(0, 0, 0, 0, viewDirs[p][2]*100, -viewDirs[p][1]*100, color="g")
    

def all_panels(fig, gs, state, label, Tpos, color_array, transStruc, atom_types, spgList):
    Rectangle = matplotlib.patches.Rectangle
    
    for p,pl in enumerate(viewDirs):
        if isinstance(pl[0],(list, np.ndarray)): 
            plane = normal(transStruc[state].cell).dot(pl[0] + state/n_steps*(pl[1] - pl[0]))
        else:
            plane = normal(transStruc[state].cell).dot(pl)
        plane = plane/la.norm(plane)
        #angles = dir2angles(plane)
        if p ==0:
            add_panel(fig,gs[0,0], plane, state, p, 'E', Tpos, color_array, transStruc)
        elif p == 1:
            add_panel(fig,gs[0,1], plane, state, p, 'W', Tpos, color_array, transStruc)
        else:
            add_panel(fig,gs[1,0], plane, state, p, 'E', Tpos, color_array, transStruc)
        
    ax = fig.add_subplot(gs[1,1], projection='3d', proj_type = 'ortho')
    ax.set_anchor('W')
    ax.set_axis_off()
    # ax.view_init(*(np.array(dir2angles(transStruc[state].cell[:,1]))+np.array([30,30])))
    ax.view_init(*(np.array(dir2angles(transStruc[state].cell[:,0]-transStruc[state].cell[:,1])) + np.array([10,10])))
    ax.dist = 5
    maxXAxis = abs(np.array([s.cell for s in transStruc])).max() + 1
    ax.set_xlim([-maxXAxis, maxXAxis])
    ax.set_ylim([-maxXAxis, maxXAxis])
    ax.set_zlim([-maxXAxis, maxXAxis])
    
    axlims = [a for b in ax.get_position().get_points() for a in b]
    rec = Rectangle((axlims[0],axlims[1]),(axlims[2]-axlims[0]),(axlims[3]-axlims[1]), transform = fig.transFigure, fill=False,lw=3, color="k")
    fig.patches.append(rec)

    origin = np.sum(transStruc[state].cell, axis=1)/2

    for i in range(3):
        base = np.array([np.zeros(3), transStruc[state].cell[:,(i+1)%3],
                         transStruc[state].cell[:,(i+2)%3], 
                         transStruc[state].cell[:,(i+1)%3] + transStruc[state].cell[:,(i+2)%3]])
        vec = transStruc[state].cell[:,i:i+1].dot(np.ones((1,4)))
        ax.quiver(base[:,0]-origin[0], base[:,1]-origin[1], base[:,2]-origin[2], vec[0,:], vec[1,:], vec[2,:], arrow_length_ratio=0, color="k", alpha=0.5)
    # origin2 = np.sum(transStruc[state].cell.dot(np.diag([5,5,5])),axis=1)/2
    # for a in supercell(transStruc[state], transStruc[state].cell.dot(np.diag([5,5,5]))):
    #     ax.scatter(a.pos[0]-origin2[0], a.pos[1]-origin2[1], a.pos[2]-origin2[2], alpha = 0.05, s=400, color=colorlist[(np.where(atom_types == a.type)[0][0])%10])
    a_list = []
    first = True
    for a in transStruc[state]:
        if first:
            first = False
            apos1 = a.pos
        if a.type in a_list or not label:
            ax.scatter(*(a.pos-origin-apos1), alpha = 1, s=200, color=colorlist[(np.where(atom_types == a.type)[0][0])%10])
        else:
            a_list.append(a.type)
            ax.scatter(*(a.pos-origin-apos1), alpha = 1, s=200, color=colorlist[(np.where(atom_types == a.type)[0][0])%10], label=a.type)
        
    for p,pl in enumerate(viewDirs):
        if isinstance(pl[0],(list, np.ndarray)): 
            plane = normal(transStruc[state].cell).dot(pl[0] + state/n_steps*(pl[1] - pl[0]))
        else:
            # plane = normal(transStruc[state].cell).dot(pl)
            plane = transStruc[state].cell.dot(pl)
        # plane = 5*plane/la.norm(plane)
        ax.quiver(*(-origin), *plane, color=reccolor[p])
        # ax.quiver(*np.zeros(3), *plane, pivot='tip', color=reccolor[p])
        
    fig.suptitle("Space group: " + spgList[state], fontsize=16)
    for x in ax.get_children():
        if isinstance(x, matplotlib.legend.Legend):
            break
    else:
        fig.legend()

def make_fig(state, Tpos, color_array, transStruc, atom_types, spgList, outdir, savedisplay):
    fig = plt.figure(figsize=[7.2,7.2])
    gs = matplotlib.gridspec.GridSpec(2, 2)
    gs.update(wspace=0.03, hspace=0.03)
    all_panels(fig,gs, state, True, Tpos, color_array, transStruc, atom_types, spgList)
    if savedisplay:
        fig.savefig(outdir+"/Trans_%d.svg"%state)

def make_anim(n_states, Tpos, color_array, transStruc, atom_types, spgList, outdir):
    fig = plt.figure(figsize=[12.8,7.2])
    gs = matplotlib.gridspec.GridSpec(2, 2)
    gs.update(wspace=0.03, hspace=0.03)
    def animate_trans(state):
        all_panels(fig,gs, state, state==0, Tpos, color_array, transStruc, atom_types, spgList)

    # animation.verbose.set_level('debug')

    plt.rcParams['animation.ffmpeg_path'] = '/home/felixt/projs/bin/ffmpeg'
    # Writer = animation.writers['ffmpeg']
    
    writer = animation.FFMpegWriter(fps=int(n_states/6.0),codec='prores', extra_args=['-loglevel', 'verbose','-f','mov'])

    anim = animation.FuncAnimation(fig, animate_trans,
                               frames=n_states+1, interval=1)
    anim.save(outdir + '/Trans3.mov', writer=writer)

def PCA(disps):
    # This is just kind of cool, but useless for now
    n = np.shape(disps)[1]
    M = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            # M[i,j] = disps[:,i].dot(disps[:,j])
            M[i,j] = la.norm(disps[:,i]-disps[:,j])

    M = np.exp(-M/(1 - M/M.max()))
    # M = np.exp(M.max() - M) - 1

    eigval,eigvec = la.eig(M)
    
    idx = np.argsort(-eigval)

    eigval = eigval[idx]
    eigvec = eigvec[:,idx]

    logdiffs = np.log(eigval[:-1]) - np.log(eigval[1:])

    n_class = np.argmax(logdiffs)+1

    plt.figure()
    plt.plot(eigval,".")
    
    plt.figure()
    plt.semilogy(eigval,".")

    plt.figure()
    plt.plot(logdiffs,".-")
    
    return n_class

def optimizationLoop(A, Acell, mulA, B, Bcell, mulB, ncell, filename): 
    """ This loop will repeat the entire minimization if no periodic cell can be found
        the final tmat is used to retile the structures (off by default) """
    tmat = np.eye(3)
    foundcell = None
    rep = 0
    while (rep < nrep and foundcell is None):
        if rep != 0:
            print("_____RESTARTING_____")
        
        rep += 1
    
        Apos, atomsA, atom_types = makeSphere(A, mulA*ncell) # Create Apos
        
        # Temporarly stretching Bcell, for tiling
        Btmp = deepcopy(B)
        Btmp.cell = tmat.dot(B.cell)
        for b in Btmp:
            b.pos = tmat.dot(b.pos)
    
        Bpos, atomsB = makeSphere(Btmp, mulB*ncell, atom_types) # Create Bpos
    
        Bpos = la.inv(tmat).dot(Bpos)
    
        assert all(mulA*atomsA == mulB*atomsB)
        atoms = mulA*atomsA
        
        Apos = np.asfortranarray(Apos)
        Bpos = np.asfortranarray(Bpos)
        t_time = time.time()
        Apos_map, Bpos, Bposst, n_map, natA, class_list, tmat, dmin, vec = tr.fastoptimization(Apos, Bpos, Acell, la.inv(Acell), mulA * la.det(Acell)/(mulB * la.det(Bcell)), atoms, filename)
        t_time = time.time() - t_time
        Bpos = np.asanyarray(Bpos)
        Apos = np.asanyarray(Apos)
        
        print("Mapping time:", t_time)
        
        class_list = class_list[:n_map]-class_list[:n_map].min()
    
        Bpos = Bpos[:,:n_map]
        Bposst = Bposst[:,:n_map]
        Apos_map = Apos_map[:,:n_map] 
    
        if abs(abs(la.det(tmat)) - abs(mulA * la.det(Acell)/(mulB * la.det(Bcell)))) > tol_vol and rep < nrep:
            print("Warning: The volume factor is wrong.")
            
        print("Looking for periodic cell...")        
        foundcell, origin = find_cell(class_list, Bposst)

        print(class_list)
        
        if foundcell is not None:
            print("Found cell!")
        else:
            print("Could not find periodic cell")

    return Apos, Apos_map, Bpos, Bposst, n_map, natA, class_list, tmat, dmin, atoms, atom_types, foundcell, origin
            
            
def p2ptrans(fileA, fileB, ncell, filename, interactive, savedisplay,
             outdir, use, switch, prim, anim, vol, minimize):
    
    on_top = None

    os.makedirs(outdir, exist_ok=True)

    if not interactive:
        matplotlib.use('Agg')

    global plt
        
    import matplotlib.pyplot as plt

    plt.rcParams["figure.figsize"] = [5, 5]

    # START
    print(" ________    _______  ________   ")   
    print("|\\   __  \\  /  ___  \\|\\   __  \\  ")
    print("\\ \\  \\|\\  \\/__/|_/  /\\ \\  \\|\\  \\ ")
    print(" \\ \\   ____\\__|//  / /\\ \\   ____\\")
    print("  \\ \\  \\___|   /  /_/__\\ \\  \\___|")
    print("   \\ \\__\\     |\\________\\ \\__\\   ")
    print("    \\|__|      \\|_______|\\|__|   ")
    print()
    print("__________TRANSFORMATIONS__________")
    print()
    
    # If reusing a result load the info
    if use:

        fileA = "File 1"
        fileB = "File 2"
        
        A, B, ncell, filecontent, switch, prim, vol = pickle.load(open(outdir+"/param.dat","rb"))

        print("==>Using information from %s<=="%(outdir))
        print()
        print("The inputs for that run were:")
        print("-----------------------------------")
        print("File 1:", A)
        print("File 2:", B)
        print("ncell:", ncell)
        print("Param file:")
        for l in filecontent:
            print(l.rstrip())
        print("switch:", switch)
        print("prim:", prim)
        print("vol:", vol)
        print("-----------------------------------")
        print()

    else:

        # Set up the output directory
        if not os.path.exists(outdir):
            os.makedirs(outdir)                
                
        # Read the structure files
        A = read.poscar(fileA)
        B = read.poscar(fileB)

        with open(filename, "r") as f:
            filecontent = f.readlines()

        # Save the parameters for the use function
        pickle.dump((A, B, ncell, filecontent, switch, prim, vol), open(outdir+"/param.dat","wb"))
        
    # Make the structure primitive (default)
    if prim:
        lenA = len(A)
        lenB = len(B)
        print("Making the cells primitive.")
        A = primitive(A, tol)
        B = primitive(B, tol)
        if len(A) == lenA:
            print("%s (%s) did not change."%(A.name, fileA))
        else:
            print("The size of %s (%s) changed from %d to %d."%(A.name, fileA, lenA, len(A)))
        if len(B) == lenB:
            print("%s (%s) did not change."%(B.name, fileB))
        else:
            print("The size of %s (%s) changed from %d to %d."%(B.name, fileB, lenB, len(B)))
        print()
            
    # Make sure they have the same number of atoms
    mul = lcm(len(A),len(B))
    mulA = mul//len(A)
    mulB = mul//len(B)

    Acell = A.cell*float(A.scale)
    Bcell = B.cell*float(B.scale)

    ccellA = ccell1
    ccellB = ccell2
    
    # Transformations is always from A to B. Unless switch is True, the more dense structure
    # is set as B
    if (abs(mulA*la.det(Acell)) < abs(mulB*la.det(Bcell))) != switch: # (is switched?) != switch
        A, mulA, Acell, fileA, ccellA, B, mulB, Bcell, fileB, ccellB = B, mulB, Bcell, fileB, ccellB, A, mulA, Acell, fileA, ccellA
    print("Transition from %s (%s) to %s (%s)"%(A.name, fileA, B.name, fileB))
    
    
    # If vol is true the volumes are made equal
    if vol:
        normalize = (abs(mulA*la.det(Acell)) / abs(mulB*la.det(Bcell)))**(1./3.)
        Bcell = normalize*Bcell
        B.cell = Bcell
        for b in B:
            b.pos = normalize*b.pos*float(B.scale)
        B.scale = 1
        
    initSpg = get_spacegroup(to_spglib(A), symprec=0.3, angle_tolerance=3.0)
    finalSpg = get_spacegroup(to_spglib(B), symprec=0.3, angle_tolerance=3.0)
            
    print("Initial SpaceGroup:", initSpg)
    print("Final SpaceGroup:", finalSpg)

    print("Number of %s (%s) cells in sphere"%(A.name, fileA), mulA*ncell)
    print("Number of %s (%s) cells in sphere:"%(B.name, fileB), mulB*ncell)
    print("Total number of atoms in each sphere:", mulA*ncell*len(A))
    print()
    
    if minimize:
        print("==>Ready to start optimmization<==")

        result = optimizationLoop(A, Acell, mulA, B, Bcell, mulB, ncell, filename)
        pickle.dump(result, open(outdir+"/fastoptimization.dat","wb"))
        
    else:
        print("==>Gathering optimization data from %s<=="%(outdir))
        result = pickle.load(open(outdir+"/fastoptimization.dat","rb"))

    Apos, Apos_map, Bpos, Bposst, n_map, natA, class_list, tmat, dmin, atoms, atom_types, foundcell, origin = result
        
    natB = n_map // np.sum(atoms)
    nat_map = n_map // np.sum(atoms)
    nat = np.shape(Apos)[1] // np.sum(atoms)

    print()
    print("-------OPTIMIZATION RESULTS--------")
    print()
    print("Number of classes:", len(np.unique(class_list)))
    print("Number of mapped atoms:", n_map)
    print("Total distance between structures:", dmin) 
    print("Volume stretching factor:", la.det(tmat))
    print("Cell volume ratio (should be exactly the same):", mulA * la.det(Acell)/(mulB * la.det(Bcell)))

    print()
    print("-----------PERIODIC CELL-----------")
    print()
    # Displacements without stretching (for plotting)
    disps_total = Apos_map - Bpos

    # Displacement with stretching
    disps = Apos_map - Bposst

    # Classes of vectors and their value
    vec_classes = np.array([np.mean(disps[:,class_list==d_type], axis=1) for d_type in np.unique(class_list)])

    # Run the PCA analysis if turned on
    if pca:
        print("PCA found %d classes"%PCA(disps))

    # Show and/or save interactive display of optimal result
    if savedisplay or interactive:
        print("Displaying Optimal Connections...")
        if interactive:
            print("(Close the display window to continue)")
        displayOptimalResult(Apos, Bpos, Bposst, disps_total, disps, class_list, vec_classes,
                                 nat, natA, natB, atoms, outdir, savedisplay, interactive)
        print()

    # Create gif if turned on
    if gif:
        makeGif(Apos, Bposst, disps, vec_classes, nat, atoms)

    # End program if a periodic cell couldn't be found earlier
    if foundcell is None:
        raise RuntimeError("Could not find good displacement cell. Increase system size")
    
    pos_in_struc = Bposst - origin.dot(np.ones((1,np.shape(Bposst)[1])))

    dispStruc, stinitStruc, finalStruc = makeStructures(foundcell, atoms, atom_types,
                                                        natB, pos_in_struc, class_list, vec_classes)

    print("Size of the transformation cell (TC):", len(dispStruc))

    print("Number of %s (%s) cells in TC:"%(A.name, fileA), abs(la.det(dispStruc.cell)/(mulA*la.det(Acell))))
    print("Number of %s (%s) cells in TC:"%(B.name, fileB), abs(la.det(dispStruc.cell)/(mulB*la.det(tmat.dot(Bcell)))))
    # Total displacement per unit volume a as metric
    Total_disp = 0
    for disp in dispStruc:
        Total_disp += la.norm(vec_classes[int(disp.type)])
    
    Total_disp = Total_disp / la.det(dispStruc.cell)
    
    print("Total displacement in stretched cell:", Total_disp)

    print()
    print("TC in %s (%s) coordinates:"%(B.name, fileB))
    printMatAndDir(la.inv(tmat).dot(dispStruc.cell), ccellB)
    print()
    print("TC in %s (%s) coordinates:"%(A.name, fileA))
    printMatAndDir(dispStruc.cell, ccellA)
    print()    

    if interactive or savedisplay:
        print("Displaying the Transition cell...")
        displayTransCell(disps, dispStruc, finalStruc, foundcell,
                         pos_in_struc, vec_classes, interactive, savedisplay)
        print()
    
    eigval, U, P, Q, planeHab = crystallography(tmat, A, B, fileA, fileB, ccellA, ccellB)

    print("=>Producing the steps along the transition<=")
    transStruc, spgList, Tpos, color_array = produceTransition(tmat, dispStruc, finalStruc,
                                                               vec_classes, outdir, atom_types,
                                                               anim, savedisplay, interactive)
                
    print("Spacegroups along the transition:")
    print(" -> ".join([e for i,e in enumerate(spgList) if i==len(spgList)-1 or e!=spgList[i+1]])) #This will remove repetitions that are next to each other
    print()
        
    # Showing some of the steps
    if interactive or savedisplay:
        print("Displaying Frames...")
        for i in range(n_frames):
            make_fig(int(i*n_steps/n_frames-1), Tpos, color_array,
                     transStruc, atom_types, spgList, outdir, savedisplay)
            
        if interactive:
            print("(Close the display windows to continue)")
            plt.show()
            print()

    if anim:
        print("Producing the animation...(this may take several hours)")
        make_anim(n_steps, Tpos, color_array, transStruc, atom_types, spgList, outdir)
        print()
        
    print("p2ptrans finished successfully")

if __name__=='__main__':
    p2ptrans(*readOptions())
