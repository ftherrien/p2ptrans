#!/usr/bin/env python3
import time
from copy import deepcopy
import os

from .config import *
from .fmodules import transform as tr
from .fmodules import tiling as t
from .display import displayStats, displayOptimalResult, displayTransCell, printMatAndDir
from .core import makeSphere, find_cell, makeStructures, switchDispStruc 
from .utils import scale, is_same_struc, lcm

def find_layer(strucList, rule):
    """This function replaces the names of the atoms according to the rule in a list of 
    structures.
    """
    A = []
    ruleset = set.union(*list(rule.values())) # Puts all possible values in one set
    for struc in strucList:
        struc = scale(struc)
        types = set([a.type for a in struc]) # Create a set of types
        if not ruleset.isdisjoint(types): # If ruleset is contained in types
            tmpStruc = Structure(struc.cell)
            tmpStruc.name = struc.name
            for a in struc:
                if a.type in ruleset:
                    for key, val in rule.items():
                        if a.type in val:
                            newtype = key
                            break
                    tmpStruc.add_atom(*a.pos, newtype) # Changes the types according to the rule 
            A.append(tmpStruc)
    return A
        
def find_basis(diruvw, ccell = np.eye(3), tol=tol, maxi=10):
    """ Finds a new cell for the structure where the specified plane is in the z direction """

    vectors = [diruvw*1000]
    inplane = []
    rangeijk = np.arange(-maxi,maxi+1)
    rangeijk = rangeijk[np.argsort(abs(rangeijk))]
    for i in rangeijk:
        for j in rangeijk:
            for k in rangeijk:
                if [i,j,k] != [0,0,0]: #Non-zero?
                    pos = ccell.dot(np.array([i,j,k]))
                    if abs(diruvw.dot(pos)) < tol: #In plane?
                        inplane.append(pos)
                    elif la.norm(np.cross(diruvw, pos)) < tol: #Is it parallel to uvw?
                        if la.norm(pos) < la.norm(vectors[0]) and pos.dot(vectors[0]) > 0: # Is it shorter?
                            vectors[0] = pos
                            
    # Make an array
    inplane = np.array(inplane).T

    idx = np.argsort(la.norm(inplane, axis = 0))
    for i in idx:
        if len(vectors) < 3 and la.norm(np.cross(vectors[-1], inplane[:,i])) > tol: #Form a basis?
            vectors.append(inplane[:,i])
        if len(vectors) == 3:
            break
    
    cell3D = np.array(vectors[::-1]).T
    
    # Expressed in terms of the cell
    inplane_cell = la.inv(cell3D).dot(inplane)

    # Vertices inside the cell
    inplane_cell = inplane_cell[:,np.sum((inplane_cell < 1 - tol) & (inplane_cell > - tol),0)==3]
    
    if len(inplane_cell[0,:]) > 0:

        # Make sure the cell is primitive (not sure this is necessary)
        
        idd = np.argmin(la.norm(cell3D[:,:2], axis=0))
    
        cell3D[:,int(not idd)] = cell3D.dot(inplane_cell[:,np.argmin(abs(inplane_cell[int(not idd),:]))])

    if la.det(cell3D) < 0:
        cell3D[:,:2] = cell3D[:,1::-1] 

    cell3D = gruber(cell3D)    
    
    cell2D = np.array([[la.norm(cell3D[:,0]),0,0],
                       [cell3D[:,0].dot(cell3D[:,1])/la.norm(cell3D[:,0]),la.norm(np.cross(cell3D[:,0],cell3D[:,1]))/la.norm(cell3D[:,0]),0],
                       [0,0,la.norm(cell3D[:,2])]]).T
    
    return cell2D, cell3D

def readSurface(A, planehkl, rule, ccell=None, tol=tol):
    """Creates a structure for all possible combinations and renames the atoms according to the rule"""

    A = scale(A)
    
    if ccell is None:
        ccell = A.cell
        
    diruvw = la.inv(ccell).dot(planehkl) # Transforms hkl into uvw

    diruvw = diruvw/la.norm(diruvw)

    cell2D, cell3D = find_basis(diruvw, ccell = ccell, tol=tol) #Create a new cell with diruvw as z
    
    A = supercell(A, cell3D) # Reorganize the atoms in the new cell
    
    idx = list(range(len(A)))
    
    inplane = []
    # Creates a list of list of indices where each list of indices correspond to a possible termination
    while len(idx)>0:
        inplane.append([idx[0]])
        for i in idx[1:]:
            if abs(diruvw.dot(A[idx[0]].pos) - diruvw.dot(A[i].pos)) < tol:
                inplane[-1].append(i)
                idx.remove(i)
        idx = idx[1:]
        
    strucs = []
    # Creates a structure for each possible termination
    for j,ix in enumerate(inplane):
        tmpstruc = Structure(cell2D)
        tmpstruc.add_atom(0,0,0,A[ix[0]].type)
        for i in ix[1:]:
            pos = cell2D.dot(la.inv(cell3D)).dot(A[i].pos-A[ix[0]].pos)
            tmpstruc.add_atom(*(pos),A[i].type)

        tmpstruc = supercell(tmpstruc,tmpstruc.cell) #Puts atoms inside the cell

        # Removes structures that are trivially the same
        for i, s in enumerate(strucs):
            for a in s:
                tmps = deepcopy(s)
                for b in tmps:
                    b.pos = b.pos - a.pos
                tmps = supercell(tmps,tmps.cell)
                if is_same_struc(tmpstruc, tmps, tol=tol):
                    break
            else:
                continue
            break
        else:
            tmpstruc.name = A.name + " " + "(%d %d %d)"%(tuple(planehkl)) + " " + str(len(strucs))
            strucs.append(tmpstruc)
            write.poscar(strucs[-1], vasp5=True, file="POSCAR_%03d"%j)

    for i in range(len(strucs)):
        strucs[i] = primitive(strucs[i])
        reshift(strucs[i])
            
    return find_layer(strucs, rule)    

def reshift(struc):
    """Reshift the z axis to the last position in the matrix and makes it positive"""
    
    idx = np.argmax(abs(struc.cell[2,:]))
    if idx == 2:
        struc.cell[2,2] = np.sign(struc.cell[2,2])*struc.cell[2,2]
        if (struc.cell[2,2] < 0) ^ (la.det(struc.cell) < 0): 
            struc.cell[:,:2] = struc.cell[:,1::-1]
    else:
        struc.cell[:, idx], struc.cell[:, 2] = -np.sign(la.det(struc.cell))*np.sign(struc.cell[2, idx])*struc.cell[:,2], np.sign(struc.cell[2,idx])*struc.cell[:, idx]
        
    return supercell(struc,struc.cell)
    

def optimization2D(A, mulA, B, mulB, ncell, n_iter, sym, filename, outdir): 
    """ Creates the structures and finds the mapping between them """
    
    Apos, atomsA, atom_types = makeSphere(A, mulA*ncell, twoD=True) # Create Apos
    
    Bpos, atomsB = makeSphere(B, mulB*ncell, atom_types, twoD=True) # Create Bpos
    
    assert all(mulA*atomsA == mulB*atomsB)
    atoms = mulA*atomsA
        
    Apos = np.asfortranarray(Apos)
    Bpos = np.asfortranarray(Bpos)
    t_time = time.time()
    print("Optimizing... (this may take several hours)")
    print("Check progress in %s"%(outdir+"/progress.txt"))
    result = tr.intoptimization(sym, n_iter, Apos, Bpos, A.cell, la.inv(A.cell),
                                atoms, filename, outdir)
    Apos_map, Bpos, Bposst, n_map, natA, class_list, ttrans, rtrans, dmin, stats, n_peaks, peak_thetas = result
    t_time = time.time() - t_time
    Bpos = np.asanyarray(Bpos)
    Apos = np.asanyarray(Apos)
        
    print("Mapping time:", t_time)
    
    class_list = class_list[:n_peaks,:n_map]-class_list[:n_peaks,:n_map].min()
    
    Bpos = Bpos[:n_peaks,:,:n_map]
    Bposst = Bposst[:n_peaks,:,:n_map]
    Apos_map = Apos_map[:n_peaks,:,:n_map]    

    ttrans = ttrans[:n_peaks,:,:]
    rtrans = rtrans[:n_peaks,:,:]
    dmin = dmin[:n_peaks]
    
    foundcell = [None]*n_peaks
    origin = [None]*n_peaks
    
    for i in range(n_peaks): 
            
        print("Looking for periodic cell for peak %d..."%(i))        

        foundcell[i] = np.eye(3)
        origin[i] = np.zeros((3,1))
        
        foundcell[i][:2,:2], origin[i][:2,:] = find_cell(class_list[i,:], Bposst[i,:2,:])

        
        
        if foundcell[i] is not None:
            print("Found cell!")
        else:
            print("Could not find periodic cell")

    return (Apos, Apos_map, Bpos, Bposst, n_map, natA, class_list,
            ttrans, rtrans, dmin, stats, n_peaks, peak_thetas,
            atoms, atom_types, foundcell, origin)

def findMatchingInterfaces(A, B, ncell, n_iter, sym=1, filename="p2p.in", interactive=False,
                           savedisplay=False, outdir='.',
                           minimize=True, test=False):

    # Make sure the scale is 1 for both structures
    A = scale(A)
    B = scale(B)
    
    os.makedirs(outdir, exist_ok=True)

    # Make sure structures have the same number of atoms
    mul = lcm(len(A),len(B))
    mulA = mul//len(A)
    mulB = mul//len(B)
    
    # Setting the unit cells of A and B in 2D
    Acell = A.cell[:2,:2]
    Bcell = B.cell[:2,:2]

    switched = False
    if abs(mulA*la.det(Acell)) > abs(mulB*la.det(Bcell)):
        A, mulA, Acell, B, mulB, Bcell = B, mulB, Bcell, A, mulA, Acell
        switched = True

    print("Number of %s cells in disk:"%(A.name), mulA*ncell)
    print("Number of %s cells in disk:"%(B.name), mulB*ncell)
    print("Total number of atoms in each disk:", mulA*ncell*len(A))
    print()

    if test:
        return

    if minimize:
        print("==>Ready to start optimmization<==")

        result = optimization2D(A, mulA, B, mulB, ncell, n_iter, sym, filename, outdir)
        pickle.dump(result, open(outdir+"/intoptimization.dat","wb"))
        
    else:
        print("==>Gathering optimization data from %s<=="%(outdir))
        result = pickle.load(open(outdir+"/intoptimization.dat","rb"))
        
    Apos, Apos_map, Bpos, Bposst, n_map, natA, class_list, ttrans, rtrans, dmin, stats, n_peaks, peak_thetas, atoms, atom_types, foundcell, origin = result

    natB = n_map // np.sum(atoms)
    nat_map = n_map // np.sum(atoms)
    nat = np.shape(Apos)[1] // np.sum(atoms)
    
    if interactive or savedisplay:
        print("Displaying Statistics...")
        if interactive:
            print("(Close the display window to continue)")
        displayStats(stats, n_iter, peak_thetas, ttrans, dmin, n_peaks, sym, interactive, savedisplay, outdir)
        print()

    dispStruc = [None]*n_peaks
    vec_classes = [None]*n_peaks
        
    for k in range(n_peaks):


        print("ANGLE!", peak_thetas[k])
        
        outcur = outdir+"/peak_%03d"%k
        
        os.makedirs(outcur, exist_ok=True)
        
        print()
        print("-------OPTIMIZATION RESULTS FOR PEAK %d--------"%k)
        print()
        print("Number of classes:", len(np.unique(class_list[k,:])))
        print("Number of mapped atoms:", n_map)
        print("Total distance between structures:", dmin[k]) 
        print("Volume stretching factor (det(T)):", la.det(ttrans[k,:2,:2]))
        print("Cell volume ratio (initial cell volume)/(final cell volume):", mulA * la.det(Acell)/(mulB * la.det(Bcell)))

        print()
        print("-----------PERIODIC CELL-----------")
        print()
        
        # Displacements without stretching (for plotting)
        bonds_total = Apos_map[k,:,:] - Bpos[k,:,:]

        # Displacement with stretching
        bonds = Apos_map[k,:,:] - Bposst[k,:,:]

        # Classes of vectors and their value
        vec_classes_estimate = np.array([np.mean(bonds[:,class_list[k,:]==d_type], axis=1) for d_type in np.unique(class_list[k,:])])

        # Show and/or save interactive display of optimal result
        if savedisplay or interactive:
            print("Displaying Optimal Connections...")
            if interactive:
                print("(Close the display window to continue)")
            displayOptimalResult(Apos, Bpos[k,:,:], Bposst[k,:,:], bonds_total, bonds, class_list[k,:],
                                 vec_classes_estimate, nat, natA, natB, atoms, outcur, savedisplay, interactive)
            print()
                         
        # End loop if a periodic cell couldn't be found earlier
        if foundcell[k] is None:
            print("Could not find good displacement cell. Increase system size")
            continue

        pos_in_struc = [None]*2
        pos_in_struc[0] = Apos_map[k,:,:] - origin[k].reshape((3,1)).dot(np.ones((1,np.shape(Apos_map[k,:,:])[1])))
        pos_in_struc[1] = Bposst[k,:,:] - origin[k].reshape((3,1)).dot(np.ones((1,np.shape(Bposst[k,:,:])[1])))

        dispStruc[k], vec_classes[k] = makeStructures(foundcell[k], atoms, atom_types,
                                                  natB, pos_in_struc, class_list[k,:])

        dispStruc[k] = reshift(dispStruc[k])

        if len(vec_classes[k]) > len(vec_classes_estimate):
            print("WARNING: The number of classes found during the optimization is not sufficent to describe the transformation; increase the classification tolerence (tol_class)")
            print()
        elif len(vec_classes[k]) < len(vec_classes_estimate):
            print("WARNING: The number of classes found during the optimization is unnecessarily large; decrease the classification tolerence (tol_class)")
            print()
        
        print("Number of bonds in Interface Cell (IC):", len(dispStruc[k]))
        
        print("Number of %s cells in IC:"%(A.name), abs(la.det(dispStruc[k].cell[:2,:2])/(la.det(Acell))))
        print("Number of %s cells in IC:"%(B.name), abs(la.det(dispStruc[k].cell[:2,:2])/(la.det(ttrans[k,:,:3].dot(B.cell)[:2,:2]))))
        print()
        print("IC in %s coordinates:"%(B.name))
        printMatAndDir(la.inv(ttrans[k,:,:3]).dot(dispStruc[k].cell), np.eye(3))
        print()
        print("IC in %s coordinates:"%(A.name))
        printMatAndDir(dispStruc[k].cell, np.eye(3))
        print()    

        if interactive or savedisplay:
            if switched:
                direction = "reverse"
            else:
                direction = "direct"
            print("Displaying the interface cell as optimized (%s)..."%(direction))
            displayTransCell(bonds, dispStruc[k], foundcell[k],
                             pos_in_struc[0], vec_classes[k], outdir, interactive, savedisplay)
            print()
            
        if switched:
            dispStruc[k], ttrans[k,:,:3], vec_classes[k] = switchDispStruc(dispStruc[k], ttrans[k,:,:3], vec_classes[k])        

    return ttrans, dispStruc, vec_classes
