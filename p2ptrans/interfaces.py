#!/usr/bin/env python3
import time
from copy import deepcopy
import os

from .config import *
from .fmodules import transform as tr
from .fmodules import tiling as t
from .display import displayStats, displayOptimalResult, displayTransCell, printMatAndDir
from .core import makeSphere, find_cell, makeStructures, switchDispStruc, find_periodicity 
from .utils import scale, is_same_struc, lcm, lccell, superize, reshift

def find_type(a, rule):
    newtype = []
    for key, val in rule.items():
        for i,v in enumerate(val):
            if v == a.type:
                newtype.append(key+"-"+str(i))
    return newtype
        
def find_basis(dircart, ccell, tol=tol, maxi=10, primitive = False):
    """ Finds a new cell for the structure where a and b are orthogonal to dircart"""

    vectors = [dircart*1000]
    tmpz = dircart*1000
    inplane = []
    rangeijk = np.arange(-maxi,maxi+1)
    rangeijk = rangeijk[np.argsort(abs(rangeijk))]
    for i in rangeijk:
        for j in rangeijk:
            for k in rangeijk:
                if [i,j,k] != [0,0,0]: #Non-zero?
                    pos = ccell.dot(np.array([i,j,k]))
                    if abs(dircart.dot(pos)) < tol: #In plane?
                        inplane.append(pos)
                    elif la.norm(np.cross(dircart, pos)) < tol: #Is it parallel to uvw?
                        if la.norm(pos) < la.norm(vectors[0]) and pos.dot(vectors[0]) > 0: # Is it shorter?
                            vectors[0] = pos
                    else:
                        if dircart.dot(pos) > 0:
                            if (primitive and (dircart.dot(pos) + tol < dircart.dot(tmpz) or (abs(dircart.dot(pos) - dircart.dot(tmpz)) < tol and la.norm(pos) < la.norm(pos)))) or (not primitive and la.norm(pos) < la.norm(tmpz)):
                                tmpz = pos
                            
    if np.allclose(vectors[0], dircart*1000, 1e-12):
        print("WARNING: Could not find lattice point in the specified z direction:", dircart)
        vectors[0] = tmpz

    if len(inplane) < 3:
        raise RuntimeError("The termination plane specified could not be found with maxi = %d"%maxi)
    
    # Make an array
    inplane = np.array(inplane).T
    
    if primitive:
        for i in range(len(inplane)):
            for j in range(i):
                area = la.norm(np.cross(inplane[:,i], inplane[:,j]))
                if area > tol:
                    if len(vectors) == 1:
                        vectors.append(inplane[:,i])
                        vectors.append(inplane[:,j])
                    elif area < la.norm(np.cross(vectors[1], vectors[2])):
                        vectors[1] = inplane[:,i]
                        vectors[2] = inplane[:,j]
    else:
        idx = np.argsort(la.norm(inplane, axis = 0))
        for i in idx:
            if len(vectors) < 3 and la.norm(np.cross(vectors[-1], inplane[:,i])) > tol: #Form a basis?
                vectors.append(inplane[:,i])
            if len(vectors) == 3:
                break
        else:
            raise RuntimeError("Could not form a basis with specified uvw with and a maximum of %d unit cells"%(maxi))
    
    cell3D = np.array(vectors[::-1]).T

    # I don't think this does anything
    # Expressed in terms of the cell
    inplane_cell = la.inv(cell3D).dot(inplane)

    # Vertices inside the cell
    inplane_cell = inplane_cell[:,np.sum((inplane_cell < 1 - tol) & (inplane_cell > - tol),0)==3]
    
    if len(inplane_cell[0,:]) > 0:

        # Make sure the cell is primitive (not sure this is necessary)
        
        idd = np.argmin(la.norm(cell3D[:,:2], axis=0))
    
        cell3D[:,int(not idd)] = cell3D.dot(inplane_cell[:,np.argmin(abs(inplane_cell[int(not idd),:]))])

    cell3D[:,2] = dircart
    cell3D = gruber(cell3D)
    
    for i in range(2):
        if np.allclose(cell3D[:,i], dircart, atol=1e-12):
            cell3D[:,i], cell3D[:,2] = -cell3D[:,2], vectors[0]
            break

        if np.allclose(cell3D[:,i], -dircart, atol=1e-12):
            cell3D[:,i], cell3D[:,2] = deepcopy(cell3D[:,2]), vectors[0]
            break
    else:
        cell3D[:,2] = vectors[0]
            
    if la.det(cell3D) < 0:
        cell3D[:,0] = -cell3D[:,0]

    vtmp = cell3D[:,1] - cell3D[:,0].dot(cell3D[:,1])*cell3D[:,0] / la.norm(cell3D[:,0])**2
        
    cell2D = np.array([[la.norm(cell3D[:,0]),0,0],
                       [cell3D[:,0].dot(cell3D[:,1])/la.norm(cell3D[:,0]),la.norm(np.cross(cell3D[:,0],cell3D[:,1]))/la.norm(cell3D[:,0]),0],
                       [cell3D[:,0].dot(cell3D[:,2])/la.norm(cell3D[:,0]),
                        cell3D[:,2].dot(vtmp)/la.norm(vtmp),
                        la.det(cell3D)/la.norm(np.cross(cell3D[:,0],cell3D[:,1]))]]).T
    
    return cell2D, cell3D

def readSurface(A, planehkl, rule, ccell=None, tol=tol, primtol=1e-3,
                maxi=10, surface = None, max_thickness = None):
    """Creates a structure for all possible combinations and renames the atoms according to the rule"""

    if max_thickness is None:
        max_thickness = tol
    
    for key, val in rule.items():
        rule[key] = []
        for v in val:
            for i, symbol in  enumerate(v):
                if not symbol.isdigit():
                    break
            if i:
                rule[key].extend([v[i:]]*int(v[:i]))
            else:
                rule[key].append(v)
                
    ruleset = set([a for b in list(rule.values()) for a in b]) # Set of atoms named in the rule

    A = scale(A)
    
    if ccell is None:
        ccell = A.cell
        
    dircart = la.inv(ccell.T).dot(planehkl) # Transforms hkl into cartesian coord

    dircart = dircart/la.norm(dircart)
    
    if surface is None:
        A = primitive(A, primtol) # Find the primitive cell

    cell2D, cell3D = find_basis(dircart, A.cell, tol=tol, maxi=maxi) #Create a new cell with dircart as z
    
    A = supercell(A, cell3D) # Reorganize the atoms in the new cell 

    z = [dircart.dot(a.pos) for a in A]
    
    idx = list(np.argsort(z))

    if surface == "bottom":
        idx = idx[::-1]
        
    inplane = []
    # Creates a list of lists of indices where each list of indices correspond to a possible termination considering the rule
    c = 0
    z = []
    while len(idx)>0:
        c=c+1
        if A[idx[0]].type in ruleset:
            inplane.append([idx[0]])
            z.append([dircart.dot(A[idx[0]].pos)])
            for i in idx[1:]:
                if abs(dircart.dot(A[idx[0]].pos) - dircart.dot(A[i].pos)) < max_thickness:
                    if A[i].type in ruleset:
                        inplane[-1].append(i)
                        z[-1].append(dircart.dot(A[i].pos))
                    idx.remove(i)

            if surface is not None:
                if surface == "bottom":
                    inplane[-1] = inplane[-1][::-1]
                break
            
        idx = idx[1:]

    if len(z) > 1:
        if z[0][-1] - z[-1][0] + cell2D[2,2] < max_thickness:
            shift = A[inplane[-1][0]].pos
            for a in A:
                a.pos = into_cell(a.pos - shift, A.cell)
            inplane[0] = inplane.pop() + inplane[0]
    elif len(z) == 0:
        print("WARNING: The structure does not contain the elements specified in the rule")
    
    strucs = []
    recMap = []
    reconstructure = [None]*len(inplane)
    # Creates a structure for each possible termination
    for j,ix in enumerate(inplane):

        tmpstruc = Structure(cell2D)
        reconstructure[j] = Structure(cell2D)
        for t in find_type(A[ix[0]], rule):
            tmpstruc.add_atom(0,0,0,t)
        for i in ix[1:]:
            pos = cell2D.dot(la.inv(cell3D)).dot(A[i].pos-A[ix[0]].pos)
            for t in find_type(A[i], rule):
                tmpstruc.add_atom(*(into_cell(pos,tmpstruc.cell)),t)
        for k,a in enumerate(A):
            pos = cell2D.dot(la.inv(cell3D)).dot(a.pos-A[ix[0]].pos)
            if (surface is None or
                ((surface == "top" and pos[2] > 0 - tol) or
                 (surface == "bottom" and pos[2] < tmpstruc[-1].pos[2] + tol))):
                reconstructure[j].add_atom(*(into_cell(pos,reconstructure[j].cell)), a.type)
                
                
        # Removes structures that are trivially the same
        for i, s in enumerate(strucs):
            for a in s:
                tmps = deepcopy(s)
                for b in tmps:
                    b.pos = into_cell(b.pos - a.pos, tmps.cell)
                if is_same_struc(tmpstruc, tmps, tol=tol):
                    for b in reconstructure[j]:
                        b.pos = into_cell(b.pos + a.pos,reconstructure[j].cell)
                    for rec in recMap[i]:
                        if is_same_struc(rec, reconstructure[j], tol=tol):
                            break
                    else:
                        recMap[i].append(reconstructure[j])
                    break
            else:
                continue
            break
        else:
            tmpstruc.name = A.name + " " + "(%d %d %d)"%(tuple(planehkl)) + " " + str(len(strucs))
            strucs.append(tmpstruc)

            recMap.append([reconstructure[j]])

    for i in range(len(strucs)):
        strucs[i] = primitive(strucs[i], primtol)
        for a in strucs[i]:
            a.type=a.type.split("-")[0]
        strucs[i] = supercell(strucs[i], reshift(strucs[i].cell))

    list3D = [cell3D.dot(la.inv(cell2D))]*len(strucs)
        
    return list(zip(strucs, recMap, list3D))
    

def optimization2D(A, mulA, B, mulB, ncell, n_iter, sym, switched, filename, outdir, max_cell_size): 
    """ Creates the structures and finds the mapping between them."""

    # Make sure the scale is 1 for both structures
    A = scale(A)
    B = scale(B)
    
    Apos, atomsA, atom_types = makeSphere(A, mulA*ncell, twoD=True) # Create Apos
    
    Bpos, atomsB = makeSphere(B, mulB*ncell, atom_types, twoD=True) # Create Bpos

    centroidA = np.mean(Apos, axis=1)
    centroidB = np.mean(Bpos, axis=1)
    
    assert all(mulA*atomsA == mulB*atomsB)
    atoms = mulA*atomsA
    
    Apos = np.asfortranarray(Apos)
    Bpos = np.asfortranarray(Bpos)
    Bpos_in = deepcopy(Bpos)
    t_time = time.time()
    print("Optimizing... (this may take several hours)")
    print("Check progress in %s"%(outdir+"/progress.txt"), flush=True)
    result = tr.intoptimization(sym, n_iter, Apos, Bpos, A.cell, la.inv(A.cell),
                                atoms, switched, filename, outdir)
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
    dmin = dmin[:n_peaks,:]
    
    foundcell = [None]*n_peaks
    origin = [None]*n_peaks
    
    for i in range(n_peaks): 

        ttrans[i,:,3] = ttrans[i,:,3] - ttrans[i,:,:3].dot(centroidB) + centroidA
        
        print("Looking for periodic cell for peak %d..."%(i))        
        
        print("Trying to find periodicity directly (find_cell):")
        origin[i] = np.zeros((3,1))
        foundcell[i], origin[i][:2,:] = find_cell(class_list[i,:], Bposst[i,:2,:], minvol=abs(la.det(B.cell[:2,:2]))/len(B))
        
        if foundcell[i] is None:

            print("Trying to find the closest tmat that is commensurate with both cells (find_periodicity)")
            
            ttrans[i,:2,:2], foundcell[i] = find_periodicity(ttrans[i,:2,:2], A.cell[:2,:2], B.cell[:2,:2], n=max_cell_size)

            origin[i] = np.zeros((3,1))
            
            if n_map < la.det(foundcell)/la.det(A.cell[:2,:2])*len(A):

                print("WARNING: The cell found is larger the the number of mapped atoms!. If you can affort it, try to optimize with a larger ncell (-n).")

        else:
            foundcell[i] = A.cell[:2,:2].dot(np.round(la.inv(A.cell[:2,:2]).dot(foundcell[i])))
            bfoundcell = np.round(la.inv(ttrans[i,:2,:2].dot(B.cell[:2,:2])).dot(foundcell[i]))
            ttrans[i,:2,:2] = foundcell[i].dot(la.inv(B.cell[:2,:2].dot(bfoundcell)))
        
        if foundcell[i] is not None:
            print("Found cell!")

            # Readjusting the vec after changing tmat 
            
            dmin_old = dmin[i]

            print("Readjusting for new tmat...", flush = True) 
            Apos_in = np.asfortranarray(Apos)
            result = tr.fixed_tmat_int(Apos_in, Bpos_in, ttrans[i,:,:], atoms, switched, filename, outdir)
            Apos_map_out, Bpos_out, Bposst_out, _, _, class_list_out, ttrans[i,:,3], dmin[i] = result
            
            ttrans[i,:,3] = ttrans[i,:,3] - ttrans[i,:,:3].dot(centroidB) + centroidA
            
            class_list[i,:] = class_list_out[:n_map]-class_list_out[:n_map].min()
            
            Bpos[i,:,:] = Bpos_out[:,:n_map]
            Bposst[i,:,:] = Bposst_out[:,:n_map]
            Apos_map[i,:,:] = Apos_map_out[:,:n_map]
            
            print("PEAK %d: Change in distance after tmat adjustment: %f"%(i, dmin[i][1]-dmin_old[1]))
            
            tmpcell = np.eye(3)
            tmpcell[:,2] = A.cell[:,2] + B.cell[:,2]
            tmpcell[:2,:2] = foundcell[i]
            foundcell[i] = deepcopy(tmpcell) 
        else:
            print("Could not find periodic cell")
            
    return (Apos, Apos_map, Bpos, Bposst, n_map, natA, class_list,
            ttrans, rtrans, dmin, stats, n_peaks, peak_thetas,
            atoms, atom_types, foundcell, origin)

def createPoscar(A, B, reconA, reconB, ttrans, dispStruc, outdir=".", layers=1, vacuum=10, tol=1e-3):
    """ Creates the interfacial POSCARS """
    
    reconA = deepcopy(reconA)
    reconB = deepcopy(reconB)
    
    max_thickness = np.array([a.pos[2] for a in A]).max()
    
    reconB.cell = ttrans[:,:3].dot(reconB.cell)
    shift = deepcopy(ttrans[:,3])
    shift[2] = 0 
    for b in reconB:
        b.pos = ttrans[:,:3].dot(b.pos) + shift

    if la.det(A.cell[:2,:2]) <  la.det(reconA.cell[:2,:2]):
        dispA = deepcopy(dispStruc.cell)
        dispA[:,2] = A.cell[:,2]
        SdA = la.inv(A.cell).dot(dispA)
        SrA = la.inv(A.cell).dot(reconA.cell)
    
        NA = lccell(SdA, SrA, tol)
        SA = la.inv(SdA).dot(NA) 
    else:
        SA = np.eye(3)

    if la.det(B.cell[:2,:2]) <  la.det(reconB.cell[:2,:2]):
        dispB = deepcopy(dispStruc.cell)
        dispB[:,2] = B.cell[:,2]
        SdB = la.inv(ttrans[:,:3].dot(B.cell)).dot(dispB)
        SrB = la.inv(ttrans[:,:3].dot(B.cell)).dot(reconB.cell)
    
        NB = lccell(SdB, SrB, tol)
        SB = la.inv(SdB).dot(NB)
    else:
        SB = np.eye(3)
          
    N = lccell(SA[:2,:2], SB[:2,:2], tol)

    newA = -deepcopy(A.cell)
    
    newA[:2,:2] = dispStruc.cell[:2,:2].dot(N)
    
    newB = deepcopy(B.cell)

    newB[:2,:2] = newA[:2,:2]
    
    for a in reconA:
        a.pos[2] = a.pos[2] - max_thickness
    
    reconA = supercell(reconA, superize(reconA.cell, newA))
        
    reconB = supercell(reconB, superize(reconB.cell, newB))
    
    write.poscar(reconA, vasp5=True, file= outdir + "/POSCAR_bottom")
    
    write.poscar(reconB, vasp5=True, file= outdir + "/POSCAR_top")

    reconA = supercell(reconA, reconA.cell.dot(np.array([[1,0,0],[0,1,0],[0,0,layers]])))
        
    reconB = supercell(reconB, reconB.cell.dot(np.array([[1,0,0],[0,1,0],[0,0,layers]])) )

    interface = Structure(newA)
    interface.cell[2,2] = reconB.cell[2,2] - reconA.cell[2,2] + ttrans[2,3] - max_thickness + vacuum
    for b in reconB:
        pos = b.pos + np.array([0,0,1])*(ttrans[2,3] - max_thickness + vacuum/2 - reconA.cell[2,2])
        interface.add_atom(*pos, b.type)

    for a in reconA:
        pos = a.pos + np.array([0,0,1])*(vacuum/2 - reconA.cell[2,2])
        interface.add_atom(*pos, a.type)

    write.poscar(supercell(interface, interface.cell), vasp5=True, file=outdir + "/POSCAR_interface")


def findMatchingInterfaces(A, B, ncell, n_iter, sym=1, filename="p2p.in", interactive=False,
                           savedisplay=False, outdir='.',
                           minimize=False, test=False, A3D=np.eye(3), B3D=np.eye(3), max_cell_size=1000):

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
        A, mulA, Acell, A3D, B, mulB, Bcell, B3D = B, mulB, Bcell, B3D, A, mulA, Acell, A3D
        switched = True

    print("Number of %s cells in disk:"%(A.name), mulA*ncell)
    print("Number of %s cells in disk:"%(B.name), mulB*ncell)
    print("Total number of atoms in each disk:", mulA*ncell*len(A))
    print()

    if test:
        return None, None, None, None

    if minimize:
        print("==>Ready to start optimization<==")

        result = optimization2D(A, mulA, B, mulB, ncell, n_iter, sym, switched, filename, outdir, max_cell_size)
        pickle.dump(result, open(outdir+"/intoptimization.dat","wb"))
        
    else:
        try:
            result = pickle.load(open(outdir+"/intoptimization.dat","rb"))
            print("==>Gathered optimization data from %s<=="%(outdir))
        except FileNotFoundError:
            result = optimization2D(A, mulA, B, mulB, ncell, n_iter, sym, switched, filename, outdir, max_cell_size)
            pickle.dump(result, open(outdir+"/intoptimization.dat","wb"))
        
    Apos, Apos_map, Bpos, Bposst, n_map, natA, class_list, ttrans, rtrans, dmin, stats, n_peaks, peak_thetas, atoms, atom_types, foundcell, origin = result
    
    natB = n_map // np.sum(atoms)
    nat_map = n_map // np.sum(atoms)
    nat = np.shape(Apos)[1] // np.sum(atoms)
    
    if interactive or savedisplay:
        print("Displaying statistics...")
        if interactive:
            print("(Close the display window to continue)")
        displayStats(stats, n_iter, peak_thetas, ttrans, dmin[:,1], n_peaks, sym, interactive, savedisplay, outdir)
        print()

    dispStruc = [None]*n_peaks
    vec_classes = [None]*n_peaks
        
    for k in range(n_peaks):

        outcur = outdir+"/peak_%03d"%k

        if savedisplay:
            os.makedirs(outcur, exist_ok=True)

        area = n_map/len(B)*la.det(B.cell[:2,:2])
            
        print()
        print("-------OPTIMIZATION RESULTS FOR PEAK %d--------"%k)
        print()
        print("Number of classes:", len(np.unique(class_list[k,:])))
        print("Number of mapped atoms:", n_map)
        print("Bonding energy between structures:")
        print("----Total bonding energy (not-strained):", dmin[k,0])
        print("--------Total bonding energy (strained):", dmin[k,1])
        print("----Average energy per bond (n-s and s):", dmin[k,0]/n_map, dmin[k,1]/n_map)
        print("---Average energy per Ang^2 (n-s and s):", dmin[k,0]/area, dmin[k,1]/area)
        print("Optimal angle between structures:", np.mod(peak_thetas[k]*180/np.pi,360/sym))
        print("Volume stretching factor (det(T)):", la.det(ttrans[k,:2,:2]))
        print("Cell volume ratio (initial cell volume)/(final cell volume):", mulA * la.det(Acell)/(mulB * la.det(Bcell)))
        eig, _ = la.eig(ttrans[k,:2,:2].T.dot(ttrans[k,:2,:2]))        
        print("In-plane stretching: %f %f"%(np.sqrt(eig[0]), np.sqrt(eig[1])))
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
            print("Displaying optimal connections...")
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

        dispStruc[k] = supercell(dispStruc[k], reshift(dispStruc[k].cell))

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
        printMatAndDir(B3D.dot(la.inv(ttrans[k,:,:3]).dot(dispStruc[k].cell)), np.eye(3))
        print()
        print("IC in %s coordinates:"%(A.name))
        printMatAndDir(A3D.dot(dispStruc[k].cell), np.eye(3))
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
            ttrans[k,:,3] = -ttrans[k,:,:3].dot(ttrans[k,:,3])
                
    return ttrans, dispStruc, vec_classes, dmin[:,0].min()
