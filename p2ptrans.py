from fmodules import transform as tr
import numpy as np
import numpy.linalg as la
from matplotlib import animation
from fmodules import tiling as t
import pickle
import time
from pylada.crystal import Structure, primitive, gruber, read, write, supercell, space_group
from copy import deepcopy
import argparse
import os
from format_spglib import from_spglib, to_spglib
from spglib import get_spacegroup
from display import displayOptimalResult, makeGif, displayTransCell, make_anim, make_fig, printMatAndDir
from utils import lcm, find_uvw, normal
from config import *
        
def readOptions():

    parser = argparse.ArgumentParser()
    parser.add_argument("-I","--initial",dest="A",type=str, default='./POSCAR_A', help="Initial Structure")
    parser.add_argument("-F","--final",dest="B",type=str, default='./POSCAR_B', help="Final Structure")
    parser.add_argument("-n","--ncell",dest="ncell",type=int, default=300, help="Number of cells to tile")
    parser.add_argument("-i","--interactive",dest="interactive",action="store_true", default=False, help="Enables interactive display")
    parser.add_argument("-d","--disp",dest="savedisplay",action="store_true", default=False, help="Saves figures")
    parser.add_argument("-p","--param", dest="filename", type=str, default='./p2p.in', help="Parameter file")
    parser.add_argument("-c","--crystal", dest="crystfile", type=str, default='./cryst.in', help="Parameter file for crystllography analysis")
    parser.add_argument("-o","--outdir",dest="outdir", type=str, default='.', help="Output directory")
    parser.add_argument("-u","--use",dest="use", type=str, default=None, help="Use previously calculated data")
    parser.add_argument("-m","--minimize",dest="minimize",action="store_true", default=False, help="Rerun minimization using the parameters stored in the folder provided with the -u option")
    parser.add_argument("-s","--switch",dest="switch",action="store_true", default=False, help="Map the larger cell on the smaller cell")
    parser.add_argument("-r","--noprim",dest="prim",action="store_false", default=True, help="Finds the primitive cell at the beginning") #TMP
    parser.add_argument("-a","--anim",dest="anim",action="store_true", default=False, help="Produce the animation") #TMP
    parser.add_argument("-v","--vol",dest="vol",action="store_true", default=False, help="Make the two (stochiometric) cells equal in volume")
    parser.add_argument("-t","--test",dest="test",action="store_true", default=False, help="Tests the input file and prepares the run, you can continue this run with the -u [directory] -r option")

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
    test = options.test
    crystfile = options.crystfile
    
    return fileA, fileB, ncell, filename, interactive, savedisplay, outdir, use, switch, prim, anim, vol, minimize, test, crystfile

def readCrystParam(crystfile):
    """ Reads the 4 possible parameters of the cystfile"""
    
    # Default values
    ccell1 = np.eye(3)
    ccell2 = np.eye(3)
    planehkl = [1,0,0]
    diruvw = [0,1,0]
    
    try:
        with open(crystfile,"r") as f:
            content = f.readlines()
    except FileNotFoundError:
            content = []

    for l in content:
        if l[0].rstrip() == "#":
            continue
        line = l.split('=')
        if len(line) == 2:
            if line[0].rstrip()=="ccell1":
                ccell1 = eval(line[1].rstrip())
            elif line[0].rstrip()=="ccell2":
                ccell2 = eval(line[1].rstrip())
            elif line[0].rstrip()=="planehkl":
                planehkl = eval(line[1].rstrip())
            elif line[0].rstrip()=="diruvw":
                diruvw = eval(line[1].rstrip())
            else:
                print("WARNING: %s is not a supported input"%(line[0].rstrip()))
        elif len(line) > 2:
            raise SyntaxError(l)

    return ccell1, ccell2, planehkl, diruvw

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
    
def uniqueclose(closest, tol):
    """ For a 3xn matrix of coordinates returns an array of unique cooridnates and their
    index in the original matrix"""
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

def rotate(icell,fcell):
    U,S,V = la.svd(icell.dot(fcell.T))
    return V.conj().T.dot(U.conj().T).real

def rot_mat(u, theta):
    u = u/la.norm(u)

    P = u.reshape((3,1)).dot(u.reshape((1,3)))
    Q = np.array([[0,-u[2],u[1]], [u[2], 0, -u[0]], [-u[1], u[0], 0]])

    return  P + (np.eye(3) - P)*np.cos(theta) + Q*np.sin(theta)

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

def makeInitStruc(dispStruc, vec_classes):
    initStruc = Structure(dispStruc.cell)
    for a in dispStruc:
        initStruc.add_atom(*(a.pos+vec_classes[int(a.type)]),a.atom)
    return initStruc

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
    incell = []

    # Goes thtough the unique atomic positions found in 1 cell
    for idx, disp in zip(*uniqueclose(cell_coord, tol)):
        for i in idx:
            if np.allclose(pos_in_struc[:,i], cell.dot(disp), atol=tol):
                incell.append((i,pos_in_struc[:,i]))
                break
        else:
            i = np.argmin(la.norm(np.array([pos_in_struc[:,j] for j in idx]),axis=1))
            incell.append((idx[i],cell.dot(disp)))

    # Adds the atoms to the structure
    atom_list = []
    for i, disp in incell:
        dispStruc.add_atom(*(tuple(disp)+(str(class_list[i]),)))
        atom_list.append(whattype(i, natB))

    # Adds a new label to each displacement called atom
    for i,a in enumerate(dispStruc):
        a.atom = atom_list[i]
        
    if la.det(cell) < 0:
       cell[:,2] = -cell[:,2] 

    # Finds a squarer cell
    cell = gruber(cell)

    dispStruc = supercell(dispStruc, cell)

    # Makes sure it is the primitive cell 
    dispStruc = primitive(dispStruc, tolerance = tol)

    return dispStruc

def strainDirs(tmat, ftf=True):

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

def crystallography(tmat, A, B, ccellA, ccellB, planehkl, diruvw, fileA="input 1", fileB="input 2", ftf=True):
    print("----------CRYSTALLOGRAPHY----------")
    print()
    eigval, U, P, Q = strainDirs(tmat, ftf=ftf)

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

def produceTransition(n_steps, tmat, dispStruc, vec_classes, outdir,
                      display):

    initStruc = makeInitStruc(dispStruc, vec_classes)

    if display:
        atom_types = np.array(list(set([a.type for a in initStruc])),np.str)
        color_array = []
        Tpos = []
    else:
        atom_types = None
        color_array = None
        Tpos = None
    
    os.makedirs(outdir+PoscarDirName, exist_ok=True)

    itmat = rotate(la.inv(tmat).dot(initStruc.cell), initStruc.cell).dot(la.inv(tmat))
    
    spgList = []
    transStruc = [] 
    for i in range(n_steps+1):
        if habit:
            curMat = find_R_RU(curMat).dot(curMat)
        else:
            curMat = (itmat-np.eye(3))*i/n_steps + np.eye(3)
        curStruc = Structure(curMat.dot(initStruc.cell))
        for j,a in enumerate(dispStruc):
            curDisp = vec_classes[int(a.type)]*i/n_steps
            curPos = curMat.dot((initStruc[j].pos - curDisp).reshape((3,1)))
            curStruc.add_atom(*(curPos.T.tolist()[0]),initStruc[j].type)
            # curStruc = supercell(curStruc, curStruc.cell)
        write.poscar(curStruc, vasp5=True, file=outdir+PoscarDirName+"/POSCAR_%03d"%i) # Write resulting structure in POSCAR
        spgList.append(get_spacegroup(to_spglib(curStruc), symprec=0.3, angle_tolerance=3.0))
        transStruc.append(curStruc)

        if display:
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
                
    return transStruc, spgList, Tpos, color_array, atom_types
        
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

def optimizationLoop(A, Acell, mulA, B, Bcell, mulB, ncell, filename, outdir): 
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
        print("Optimizing... (this may take several hours)")
        print("Check progress in %s"%(outdir+"/progress.txt"))
        result = tr.fastoptimization(Apos, Bpos, Acell, la.inv(Acell),
                                     mulA * la.det(Acell)/(mulB * la.det(Bcell)),
                                     atoms, filename, outdir)
        Apos_map, Bpos, Bposst, n_map, natA, class_list, tmat, dmin, vec = result
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
        
        if foundcell is not None:
            print("Found cell!")
        else:
            print("Could not find periodic cell")

    return Apos, Apos_map, Bpos, Bposst, n_map, natA, class_list, tmat, dmin, atoms, atom_types, foundcell, origin
            
            
def p2ptrans(A, B, ncell,
             fileA='Input 1', fileB='Input 2',
             ccellA=np.eye(3), ccellB=np.eye(3), 
             filename='p2p.in',
             interactive=False, savedisplay=False,
             outdir='output',
             use=False, switch= False, prim=True,
             vol=False, minimize=True, test= False):

    os.makedirs(outdir, exist_ok=True)

    setplt(interactive)
    
    # Make the structure primitive (default)
    if prim:
        lenA = len(A)
        lenB = len(B)
        print("Making the cells primitive.")
        A = primitive(A, tol)
        B = primitive(B, tol)
        if len(A) == lenA:
            print("%s (%s) did not change. It has %d atoms."%(A.name, fileA, len(A)))
        else:
            print("The size of %s (%s) changed from %d to %d."%(A.name, fileA, lenA, len(A)))
        if len(B) == lenB:
            print("%s (%s) did not change. It has %d atoms."%(B.name, fileB, len(B)))
        else:
            print("The size of %s (%s) changed from %d to %d."%(B.name, fileB, lenB, len(B)))
        print()
        
    # Make sure they have the same number of atoms
    mul = lcm(len(A),len(B))
    mulA = mul//len(A)
    mulB = mul//len(B)

    Acell = A.cell*float(A.scale)
    Bcell = B.cell*float(B.scale)
    
    # Transformations is always from A to B. Unless switch is True, the more dense structure
    # is set as B
    switched = False
    if (abs(mulA*la.det(Acell)) < abs(mulB*la.det(Bcell))) != switch: # (is switched?) != switch
        A, mulA, Acell, fileA, ccellA, B, mulB, Bcell, fileB, ccellB = B, mulB, Bcell, fileB, ccellB, A, mulA, Acell, fileA, ccellA
        switched = True
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

    print("Number of %s (%s) cells in sphere:"%(A.name, fileA), mulA*ncell)
    print("Number of %s (%s) cells in sphere:"%(B.name, fileB), mulB*ncell)
    print("Total number of atoms in each sphere:", mulA*ncell*len(A))
    print()

    if test:
        return None, None, None, None
    
    if minimize:
        print("==>Ready to start optimmization<==")

        result = optimizationLoop(A, Acell, mulA, B, Bcell, mulB, ncell, filename, outdir)
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

    dispStruc = makeStructures(foundcell, atoms, atom_types,
                               natB, pos_in_struc, class_list, vec_classes)
    
    print("Size of the transformation cell (TC):", len(dispStruc))

    print("Number of %s (%s) cells in TC:"%(A.name, fileA), abs(la.det(dispStruc.cell)/(la.det(Acell))))
    print("Number of %s (%s) cells in TC:"%(B.name, fileB), abs(la.det(dispStruc.cell)/(la.det(tmat.dot(Bcell)))))
    
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
        displayTransCell(disps, dispStruc, foundcell,
                         pos_in_struc, vec_classes, outdir, interactive, savedisplay)
        print()

    return switched, tmat, dispStruc, vec_classes
        
def main():

    (fileA, fileB, ncell, filename, interactive, savedisplay, outdir,
     use, switch, prim, anim, vol, minimize, test, crystfile) = readOptions()

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
        print("Param file (at time of running):")
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

        A = read.poscar(fileA)
        B = read.poscar(fileB)
            
        try:
            with open(filename, "r") as f:
                filecontent = f.readlines()
        except FileNotFoundError:
            filecontent = ""

        # Save the parameters for the use function
        pickle.dump((A, B, ncell, filecontent, switch, prim, vol), open(outdir+"/param.dat","wb"))

    print("=>Reading crsytal analysis parameters<=") 
    ccell1, ccell2, planehkl, diruvw = readCrystParam(crystfile)
    print()
          
    switched, tmat, dispStruc, vec_classes = p2ptrans(A, B, ncell, fileA=fileA, fileB=fileB,
                                                      ccellA=ccell1, ccellB=ccell2,
                                                      filename=filename, interactive=interactive,
                                                      savedisplay=savedisplay, outdir=outdir,
                                                      use=use, switch=switch, prim=prim, vol=vol,
                                                      minimize=minimize, test=test)
    if test:
        return
    
    if switched:
        eigval, U, P, Q, planeHab = crystallography(la.inv(tmat), B, A, ccell2, ccell1, planehkl,
                                                    diruvw, fileA=fileB, fileB=fileA,)
    else:
        eigval, U, P, Q, planeHab = crystallography(tmat, A, B, ccell1, ccell2, planehkl,
                                                    diruvw, fileA=fileA, fileB=fileB,)
        
    print("=>Producing the steps along the transition<=")
    result = produceTransition(n_steps, tmat, dispStruc, vec_classes,
                               outdir, anim or savedisplay or interactive)
    transStruc, spgList, Tpos, color_array, atom_types = result
    
    print("Spacegroups along the transition:")
    print(" -> ".join([e for i,e in enumerate(spgList) if i==len(spgList)-1 or e!=spgList[i+1]]))
    #This will remove repetitions that are next to each other
    print()
        
    # Showing some of the steps
    if interactive or savedisplay:
        print("Displaying Frames...")
        for i in range(n_frames):
            make_fig(int(i*n_steps/(n_frames-1)), Tpos, color_array,
                     transStruc, atom_types, spgList, outdir, savedisplay, False)
            
        if interactive:
            print("(Close the display windows to continue)")
            plt.show()
            print()
    
    if anim:
        print("Producing the animation...(this may take several hours)")
        make_anim(n_steps, Tpos, color_array, transStruc, atom_types, spgList, outdir)
        print()
        
    print("p2ptrans finished successfully")

if __name__=="__main__":
    main()
