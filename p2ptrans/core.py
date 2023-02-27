import time
from copy import deepcopy
import os
from spglib import get_spacegroup
from scipy.optimize import linear_sum_assignment

from .config import *
from .fmodules import transform as tr
from .fmodules import hungarian as lap
from .format_spglib import from_spglib, to_spglib
from .fmodules import tiling as t
from .display import displayOptimalResult, makeGif, displayTransCell, make_anim, make_fig, printMatAndDir
from .utils import lcm, find_uvw, normal, rotate, PCA, makeInitStruc, reshift
from .analysis import findR, strainDirs


def find_cell(class_list, positions, tol = 1e-5, frac_shell = 0.5, frac_correct = 0.95, max_count=1000, minvol = 1e-5):
    
    cm = np.mean(positions, axis=1).reshape((np.shape(positions)[0],1)) # Centroid

    for loop in [0,1]:
        missed = 0
        for i in np.unique(class_list):
            pos = positions[:, class_list == i] # All positions of type i
            center = np.argmin(la.norm(pos, axis = 0)) # Index of most centered atom
            
            list_in = list(range(np.shape(pos)[1])) # List of indices of positions
            list_in.remove(center) # Remove the center atoms
            origin = pos[:,center:center+1] # Position of the most centered atom (the origin)
            pos = pos[:,list_in] - origin.dot(np.ones((1,np.shape(pos)[1]-1))) # Position relative to the origin
            if not loop: # First loop
                # Loop goes through in order of distances
                # Goes through all z>y, for all y>x, for all x
                norms = la.norm(pos, axis = 0)
                idx = np.argsort(norms)
                minj = 0
                maxj = len(idx)
            else:
                # Loop goes through randomly
                # Goes through all x, y smaller than z for all z 
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
                    if np.shape(pos)[0] == 3:
                        if not loop:
                            minl = k+1
                            maxl = len(idx)
                        else:
                            minl = 0
                            maxl = k-1
                    elif np.shape(pos)[0] == 2:
                        minl = 0
                        maxl = 1
                    for l in range(minl, maxl):
                        # creates all possible cells for each combinations
                        if np.shape(pos)[0] == 3:
                            newcell=np.concatenate([pos[:,idx[j]:idx[j]+1], 
                                                    pos[:,idx[k]:idx[k]+1], 
                                                    pos[:,idx[l]:idx[l]+1]],axis=1)
                        elif np.shape(pos)[0] == 2:
                            newcell=np.concatenate([pos[:,idx[j]:idx[j]+1], 
                                                    pos[:,idx[k]:idx[k]+1]],axis=1)

                        if abs(la.det(newcell)) > minvol: # Cell as non-zero volume
                            count += 1
                            if count > max_count: # Stops if reached max tries
                                break

                            if la.det(newcell) < 0:
                                newcell[:,:2] = newcell[:,1::-1] # Make sure volume is positive

                            # Norms of all positions of all types of atoms with respect to centroid
                            norms = la.norm(positions - cm.dot(np.ones((1,np.shape(positions)[1]))), axis=0)
                            # Expressed in coordinates of the new cell
                            apos = la.inv(newcell).dot(positions - origin.dot(np.ones((1,np.shape(positions)[1]))))
                            
                            # Positions and types inside the cell
                            inPos = apos[:,np.sum((apos < 1 - tol) & (apos > - tol),0)==np.shape(apos)[0]]
                            inType = class_list[np.sum((apos < 1 - tol) & (apos > - tol),0)==np.shape(apos)[0]]
                            
                            n_map = 0
                            for m, a in enumerate(apos.T):
                                # Can this pos and type be found in inPos?
                                for n, b in enumerate(inPos.T):
                                    if (all(abs(np.mod(a+tol,1)-tol-b) < tol) and 
                                        inType[n] == class_list[m]):
                                        break
                                else:
                                    continue
                                n_map += 1
        
                            genPos = []
                            # If the fractions on generated atoms is greater than frac_correct
                            if float(n_map)/float(len(class_list)) >= frac_correct:
                                # Generate a parallelepieped that contains the apos sphere
                                xMax = int(np.max(apos[0,:]))+1
                                xMin = int(np.min(apos[0,:]))-1
                                yMax = int(np.max(apos[1,:]))+1
                                yMin = int(np.min(apos[1,:]))-1
                                if np.shape(apos)[0] == 3:
                                    zMax = int(np.max(apos[2,:]))+1
                                    zMin = int(np.min(apos[2,:]))-1
                                elif np.shape(apos)[0] == 2:
                                    zMax = 1
                                    zMin = 0
                                for x in range(xMin,xMax):
                                    for y in range(yMin,yMax):
                                        for z in range(zMin,zMax):
                                            if np.shape(apos)[0] == 3:
                                                genPos.append(inPos + np.array([[x,y,z]]).T.dot(np.ones((1,np.shape(inPos)[1]))))
                                            if np.shape(apos)[0] == 2:
                                                genPos.append(inPos + np.array([[x,y]]).T.dot(np.ones((1,np.shape(inPos)[1]))))

                                genPos = newcell.dot(np.concatenate(genPos,axis=1)) #in coordinates of newcell
                                genPos = genPos + origin.dot(np.ones((1,np.shape(genPos)[1])))

                                # Cut a smaller sphere in genPos
                                genSphere = np.sum(la.norm(genPos - cm.dot(np.ones((1,np.shape(genPos)[1]))), axis = 0) < frac_shell * np.max(norms) - tol)
                                # Cut the same sphere in Apos
                                aposSphere = np.sum(norms < frac_shell * np.max(norms) - tol)
                                
                                if  genSphere == aposSphere:
                                    if missed > 0:
                                        print("Could not find periodic cell using %d of the displacements."%(missed))
                                    return newcell, origin
                                # else:
                                    # print("There should be %d atoms within %f but the cell produces %d (frac_shell)"%(aposSphere, frac_shell*np.max(norms), genSphere)) 
                            # else:
                                # print(count, "Not enough atoms are mapped to this cell %f (frac_correct = %f)"%(float(n_map)/float(len(class_list)), frac_correct))
                    else:
                        continue
                    break
                else:
                    continue
                break
            missed+=1
        if not loop:    
            print("Could not find cell using order of shortest distances, trying random order.") 
    print("Could not find periodic cell for any displacement using find_cell.")                
    return None, None

def find_periodicity(tmat, Acell, Bcell, ratio=None, n=1000):

    if ratio is not None:
        ratio = int(ratio)
    
    if len(Acell) == 3:
        A = t.sphere(Acell, n, [0,0,0])
        B = t.sphere(tmat.dot(Bcell), n, [0,0,0])

        Bin = np.round(la.inv(tmat.dot(Bcell)).dot(B)).astype(int)
    else:
        Acell3 = np.eye(3)
        Bcell3 = np.eye(3)

        Acell3[:2,:2] = Acell
        Bcell3[:2,:2] = tmat.dot(Bcell)
        
        A = t.circle(Acell3, n, [0,0,0])
        B = t.circle(Bcell3, n, [0,0,0])
        
        Bin = np.round(la.inv(Bcell3).dot(B)).astype(int)

    cost = tr.closest(np.asfortranarray(A), np.asfortranarray(B))

    cost_orig = cost
    
    # cost = cost / np.maximum(la.norm(A, axis=0).reshape((-1,1)).dot(np.ones((1,A.shape[1]))).T, la.norm(B, axis=0).reshape((-1,1)).dot(np.ones((1,B.shape[1]))))
    
    idx = np.argsort(cost, axis = None)
    
    volA = la.det(Acell)
    for ii, i in enumerate(idx):
        ia = i//n
        ib = i%n
        for jj, j in enumerate(idx[:ii]):
            ja = j//n
            jb = j%n
            if len(Acell) == 3:
                if la.norm(np.cross(A[:,ia], A[:,ja])) > tol and la.norm(np.cross(B[:,ia], B[:,ja])) > tol:
                    for kk, k in enumerate(idx[:jj]):
                        ka = k//n
                        kb = k%n
                        newcellA = np.concatenate([A[:,ia:ia+1], 
                                                   A[:,ja:ja+1],
                                                   A[:,ka:ka+1]],axis=1)
                        newcellB = np.concatenate([Bin[:,ib:ib+1], 
                                                   Bin[:,jb:jb+1], 
                                                   Bin[:,kb:kb+1]],axis=1)
                        detA = la.det(newcellA)/volA
                        detB = la.det(newcellB)
                        if abs(detA) > 1 - tol and abs(detB) > 1 - tol and (ratio is None or int(round(detB/detA)) == ratio):
                            break
                    else:
                        continue
                    break
            else:
                newcellA = np.concatenate([A[:2,ia:ia+1], 
                                           A[:2,ja:ja+1]],axis=1)
                newcellB = np.concatenate([Bin[:2,ib:ib+1], 
                                           Bin[:2,jb:jb+1]],axis=1)
            
                detA = la.det(newcellA)/volA
                detB = la.det(newcellB)
                if abs(detA) > 1 - tol and abs(detB) > 1 - tol and (ratio is None or int(round(detB/detA)) == ratio):
                    break
        else:
            continue
        break
    else:
        print("WARNING: Could not find any suitable cell.")
        return tmat, None
    
    # print("Final distances:", la.norm(newcellA - tmat.dot(Bcell).dot(newcellB), axis=0))

    tmat = newcellA.dot(la.inv(Bcell.dot(newcellB)))

    if len(Acell) == 2:
        tmpcell = np.eye(3)
        tmpcell[:2,:2] = newcellA
        tmpcell = reshift(gruber(tmpcell))
        newcellA = tmpcell[:2,:2]
    else:
        newcellA = gruber(newcellA)

    if la.det(newcellA) < 0:
       newcellA[:,0] = -newcellA[:,0]
    
    return tmat, newcellA
    
def makeSphere(A, ncell, *atom_types, twoD=False):

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
        atom_types = np.array([], str)
        atomsA = np.array([], int)
        for a in A:
            if any(atom_types == a.type):
                idx = np.where(atom_types == a.type)[0][0]
                if twoD:
                    Apos[idx] = np.concatenate((Apos[idx], t.circle(Acell, ncell, a.pos*float(A.scale))), axis = 1)
                else:
                    Apos[idx] = np.concatenate((Apos[idx], t.sphere(Acell, ncell, a.pos*float(A.scale))), axis = 1) 

                # Order the atoms with respect to distance
                Apos[idx] = Apos[idx][:,np.argsort(la.norm(Apos[idx],axis=0))] 
                atomsA[idx] += 1
            else:
                if twoD:
                    Apos.append(t.circle(Acell, ncell, a.pos*float(A.scale)))
                else:
                    Apos.append(t.sphere(Acell, ncell, a.pos*float(A.scale)))
                atom_types = np.append(atom_types, a.type)
                atomsA = np.append(atomsA,1)
        
        return (np.concatenate(Apos, axis=1), atomsA, atom_types)

    elif len(atom_types) == 1 and type(atom_types[0]) == np.ndarray:

        atom_types = atom_types[0]
        Apos = [None]*len(atom_types)
        atomsA = np.zeros(len(atom_types), int)
        for a in A:
            idx = np.where(atom_types == a.type)[0][0]
            if atomsA[idx] == 0:
                if twoD:
                    Apos[idx] = t.circle(Acell, ncell, a.pos*float(A.scale))
                else:
                    Apos[idx] = t.sphere(Acell, ncell, a.pos*float(A.scale))
            else:
                if twoD:
                    Apos[idx] = np.concatenate((Apos[idx], t.circle(Acell, ncell, a.pos*float(A.scale))), axis = 1)
                else:
                    Apos[idx] = np.concatenate((Apos[idx], t.sphere(Acell, ncell, a.pos*float(A.scale))), axis = 1)
                    # Order the atoms with respect to distance
                Apos[idx] = Apos[idx][:,np.argsort(la.norm(Apos[idx],axis=0))]
            atomsA[idx] += 1
        
        return (np.concatenate(Apos, axis=1), atomsA)

    else:
        raise ValueError("The only optional argument must be atom_types")

def uniqueclose(closest, tol):
    """ For a 3xn matrix of coordinates returns an array of unique coordinates and their
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
    return (idx, unique)

def makeStructures(cell, atoms, atom_types, natB, pos_in_struc, class_list):
    """Make the displacement structure from the repeating unit cell"""
    
    # Small function to determine type
    def whattype(pos, nat):

        pos = pos//nat + 1

        atom_tot = np.sum(np.triu(atoms.reshape((len(atoms),1)).dot(np.ones((1,len(atoms))))), axis=0)
        
        return atom_types[np.nonzero(atom_tot >= pos)[0][0]]
    
    # Make a pylada structure
    cell_coord = np.mod(la.inv(cell).dot(pos_in_struc[1])+tol,1)-tol
    cell_coord_dest = np.mod(la.inv(cell).dot(pos_in_struc[0])+tol,1)-tol
    dispStruc = Structure(cell)
    incell = []
    dest = []

    dest_id,_ = uniqueclose(cell_coord_dest, tol)
    start_id,_ = uniqueclose(cell_coord, tol)
    
    destinations = cell.dot(cell_coord_dest[:,[idx[0] for idx in dest_id]])

    start_locs = np.array([idx[0] for idx in start_id], dtype=int)

    size = max(len(dest_id), len(start_id))
    
    cost = np.zeros((size, size))

    cost_index = np.zeros((size, size), dtype=int) - 1

    # Finds the best mapping among the possible combinations of mappings
    mappings_start = [[] for i in range(len(dest_id))]
    for i, group_dest in enumerate(dest_id):
        for k, group_start in enumerate(start_id):
            inter = set(group_dest).intersection(set(group_start))
            if len(inter) > 0:                    
                mappings_start[i].append(k)
                min_dist = None
                for j in inter:
                    dist = la.norm(pos_in_struc[1][:,j] - pos_in_struc[0][:,j])
                    if dist is not None or dist < min_dist:
                        min_dist = dist
                        index = j
                cost[i,k] = min_dist
                cost_index[i,k] = index
            else:
                cost[i,k] = 1000
  
    dist, mapping = lap.munkres(cost)
    
    mapping = mapping - 1
    
    incell = np.array([cost_index[i,j] for j,i in enumerate(mapping)], dtype = int)
    
    incell = incell[incell!=-1]
        
    # Adds the atoms to the structure and creates a new vec_classes based on the actual
    # atomic positions in the final structure
    atom_list = []
    vec_classes = [None]*len(set(class_list))
    for i in incell:
        disp = pos_in_struc[1][:,i]
        vec = pos_in_struc[0][:,i] - pos_in_struc[1][:,i]
        ivec = class_list[i]
        if vec_classes[ivec] is None:
            vec_classes[ivec] = vec 
        elif not np.allclose(vec, vec_classes[ivec], atol=tol):
            vec_classes.append(vec)
            ivec = len(vec_classes)-1
                
        dispStruc.add_atom(*(tuple(disp)+(str(ivec),)))
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
    
    return dispStruc, vec_classes

def switchDispStruc(dispStruc, tmat, vec_classes):
    """ Switches the direction of the output so that the transformation is in the order
    provided by the user"""
    
    initStruc = makeInitStruc(dispStruc, vec_classes)
    newStruc = Structure(la.inv(tmat).dot(initStruc.cell))
    for j,a in enumerate(dispStruc):
        vec_classes[int(a.type)] = la.inv(tmat).dot(a.pos) - la.inv(tmat).dot(initStruc[j].pos)
        newStruc.add_atom(*(la.inv(tmat).dot(initStruc[j].pos)),a.type)
        newStruc[j].atom = initStruc[j].type
            
    return newStruc, la.inv(tmat), vec_classes
        
    
def produceTransition(n_steps, tmat, dispStruc, vec_classes, outdir,
                      display, habit=habit, a_name='Initial', b_name='Final'):

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

    if habit is None:
        itmat = rotate(la.inv(tmat).dot(initStruc.cell), initStruc.cell).dot(la.inv(tmat))
    else:
        _, itmat, _, _ = strainDirs(la.inv(tmat))
        
        
    spgList = []
    transStruc = [] 
    for i in range(n_steps+1):
        curMat = (itmat-np.eye(3))*i/n_steps + np.eye(3)
        if habit is not None:
            curMat = findR(curMat)[habit].dot(curMat)
        curStruc = Structure(curMat.dot(initStruc.cell))
        for j,a in enumerate(dispStruc):
            curDisp = vec_classes[int(a.type)]*i/n_steps
            curPos = curMat.dot((initStruc[j].pos - curDisp).reshape((3,1)))
            curStruc.add_atom(*(curPos.T.tolist()[0]),initStruc[j].type)
            curStruc.name = '{:.0%} {}, {:.0%} {}'.format(1 - i/n_steps, a_name, i/n_steps, b_name)
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

def optimization(A, Acell, mulA, B, Bcell, mulB, ncell, filename, outdir, max_cell_size): 
    """ Creates the spheres, runs the optimization in finds the cell """

    foundcell = None
    
    Apos, atomsA, atom_types = makeSphere(A, mulA*ncell) # Create Apos
    
    Bpos, atomsB = makeSphere(B, mulB*ncell, atom_types) # Create Bpos
    
    assert all(mulA*atomsA == mulB*atomsB)
    atoms = mulA*atomsA
    
    Apos = np.asfortranarray(Apos)
    Bpos = np.asfortranarray(Bpos)
    Bpos_in = deepcopy(Bpos)
    t_time = time.time()
    print("Optimizing... (this may take several hours)")
    print("Check progress in %s"%(outdir+"/progress.txt"), flush=True)
    result = tr.fastoptimization(Apos, Bpos, Acell, la.inv(Acell),
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

    tmat_old = tmat
        
    print("Trying to find periodicity directly (find_cell):")        
    foundcell, _ = find_cell(class_list, Bposst, minvol=abs(la.det(Bcell))/len(B), frac_correct = 0.95)

    if foundcell is None:

        print("Trying to find the closest tmat that is commensurate with both cells (find_periodicity)")
        
        tmat, foundcell = find_periodicity(tmat, Acell, Bcell, n = max_cell_size) 
        
        if n_map < la.det(foundcell)/la.det(Acell)*len(A):
        
            print("WARNING: The cell found is larger the the number of mapped atoms! Initial and final structures will have missing atoms. If you can't affort to optimize with a larger ncell, you can find the mapping with the current tmat on a larger system by setting map_ncell (-cn)")

    else:
        foundcell = Acell.dot(np.round(la.inv(Acell).dot(foundcell)))
        bfoundcell = np.round(la.inv(tmat.dot(Bcell)).dot(foundcell))
        tmat = foundcell.dot(la.inv(Bcell.dot(bfoundcell)))
        
    if foundcell is not None:

        old_distp =  np.sum(la.norm(Bposst - Apos_map, axis=0))

        Bposst = tmat.dot(la.inv(tmat_old)).dot(Bposst - vec.reshape((3,1)).dot(np.ones((1,Bposst.shape[1]))))
        
        vec = tr.optimize_vec(Apos_map, Bposst, vec, filename) 

        Bposst = Bposst + vec.reshape((3,1)).dot(np.ones((1,Bposst.shape[1])))
        
        print("Changed stretched dist from %f to %f (%f %%)"%(old_distp, np.sum(la.norm(Bposst - Apos_map, axis=0)), (np.sum(la.norm(Bposst - Apos_map, axis=0)) - old_distp)/old_distp*100))
        
        print("Found cell!")
        if np.any(abs(tmat-tmat_old) > tol):
            print("WARNING: tmat changed by more then the set tolerence (diff: %e, tol: %e). Visually inspect the matrices below, if the differences are small and the change in stretched distance (above) is small you can often ignore this warning."%(np.max(abs(tmat-tmat_old)),tol))
            print("Optimized tmat:")
            print(tmat_old)
            print("Exact tmat:")
            print(tmat)
        if abs(abs(la.det(tmat)) - abs(mulA * la.det(Acell)/(mulB * la.det(Bcell)))) > 1e-12:
            print("WARNING: the optimal mapping is not periodic. This might be physical, check the actual mapping using the interactive mode. A transformation *involving vacancies* will be produced.") 
        
    else:
        print("WARNING: Could not find periodic cell")

    return Apos, Apos_map, Bpos, Bposst, n_map, natA, class_list, tmat, dmin, atoms, atom_types, foundcell, vec
                       
def findMatching(A, B, ncell,
                 fileA='Input 1', fileB='Input 2',
                 ccellA=np.eye(3), ccellB=np.eye(3), 
                 filename='p2p.in',
                 interactive=False, savedisplay=False,
                 outdir='output',
                 switch= False, prim=True,
                 vol=False, minimize=False, test= False, primtol=primtol, map_ncell = None, max_cell_size = 1000):
    """ 
    This function finds the best matching between two given structures and returns the transformation matrix and the transformation cell.

    Parameters: 
    A,B (Structure): Crystal structures in the Pylada format  

    ncell (int): Minimal number of unit cells in the set of point (automatically adjusted so that there is the same number of atoms in each sphere)

    Optional Parameters:
    fileA, fileB (str): Location of file associated with the structure (display purposes only)

    ccellA, ccellB (np.array((3,3))): Conventional cell of A and B in cartesian coord

    filename (str): Location of the fortran namelist containing the parameters for the optimization

    interactive (bool): Turn on interactive mode

    savedisplay (bool): Save displayed images

    outdir (str): Location of the output directory

    switch (bool): Switch the direction of the transformation during the minimization. By default the minimization will be done from the structure with the lowest density to the structure with the highest desnity.

    prim (bool): Make the structures primitive before creating the spheres

    vol (bool): Make the volumes of the spheres equal

    minimize (bool): Run the minimization (optimization()) or read a previous result

    test (bool): Display all the parameter without starting the optimization

    Returns:
    (np.array((3,3)): Transformation matrix from Final to Initial  
    
    (Structure): Transformtation cell. a.type indicates the types of displacement, a.atom indicates the atom
    
    (list): List of 3D vectors corresponding to different types of displacements

    (np.array(3)): Final distance (dmin), Estimate of G_1, estimate of K_1
    """
    
    os.makedirs(outdir, exist_ok=True)
    
    # Make the structure primitive (default)
    if prim:
        lenA = len(A)
        lenB = len(B)
        print("Making the cells primitive.")
        A = primitive(A, primtol)
        B = primitive(B, primtol)
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
        print("The optimization will be performed in the reversed direction.")
    else:
        print("The optimization will be performed in the provided direction.")
    print("Transition is *optimized* from %s (%s) to %s (%s)"%(A.name, fileA, B.name, fileB))
    
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
            
    print("Initial spacegroup in optimization:", initSpg)
    print("Final spacegroup in optimization:", finalSpg)

    print("Number of %s (%s) cells in sphere:"%(A.name, fileA), mulA*ncell)
    print("Number of %s (%s) cells in sphere:"%(B.name, fileB), mulB*ncell)
    print("Total number of atoms in each sphere:", mulA*ncell*len(A))
    if (mulA*ncell*len(A) > 2000):
        print("WARNING: You are attempting to optimize a very large number of atoms (>2000). If this is not what you intended, reduce the minimal number of unit cells (default: 300) using the `-n` flag.")
    print()

    if test:
        return None, None, None, None
    
    if minimize:
        print("==>Ready to start optimization<==")

        result = optimization(A, Acell, mulA, B, Bcell, mulB, ncell, filename, outdir, max_cell_size)
        pickle.dump(result, open(outdir+"/fastoptimization.dat","wb"))
        
    else:
        try:
            result = pickle.load(open(outdir+"/fastoptimization.dat","rb"))
            print("==>Gathered optimization data from %s<=="%(outdir))
        except FileNotFoundError:
            result = optimization(A, Acell, mulA, B, Bcell, mulB, ncell, filename, outdir, max_cell_size)
            pickle.dump(result, open(outdir+"/fastoptimization.dat","wb"))

    Apos, Apos_map, Bpos, Bposst, n_map, natA, class_list, tmat, dmin, atoms, atom_types, foundcell, vec = result
        
    if map_ncell is not None:

        # Remapping with the new exact tmat
        
        Apos, atomsA, atom_types = makeSphere(A, Apos.shape[1]/len(A)) # recreate original Apos
        
        Bpos, atomsB = makeSphere(B, Apos.shape[1]/len(B), atom_types) # recreate original Bpos
        
        center_old =  tmat.dot(np.mean(Bpos, axis=1)) - np.mean(Apos, axis=1) # Diff between centers 
            
        Apos, atomsA, atom_types = makeSphere(A, mulA*map_ncell) # Create Apos
        
        Bpos, atomsB = makeSphere(B, mulB*map_ncell, atom_types) # Create Bpos

        center_new =  tmat.dot(np.mean(Bpos, axis=1)) - np.mean(Apos, axis=1)
        
        atoms = mulA*atomsA

        dmin_old = dmin
        
        vec = center_new - center_old + vec # Adjusting old vec to new centroids
        
        Apos_in = np.asfortranarray(Apos)
        Bpos_in = np.asfortranarray(Bpos)
        result = tr.fixed_tmat(Apos_in, Bpos_in, tmat, vec, atoms, filename, outdir)
        Apos_map, Bpos, Bposst, n_map, natA, class_list, tmat, dmin, vec = result
        
        class_list = class_list[:n_map]-class_list[:n_map].min()
        
        Bpos = Bpos[:,:n_map]
        Bposst = Bposst[:,:n_map]
        Apos_map = Apos_map[:,:n_map]

        new_result = Apos, Apos_map, Bpos, Bposst, n_map, natA, class_list, tmat, dmin, atoms, atom_types, foundcell, vec

        pickle.dump(new_result, open(outdir+"/fastoptimization.dat","wb"))
        
        print("Change in distance after tmat adjustment: %f"%(dmin[0]-dmin_old[0]))
    
    
    natB = n_map // np.sum(atoms)
    nat_map = n_map // np.sum(atoms)
    nat = np.shape(Apos)[1] // np.sum(atoms)

    print()
    print("-------OPTIMIZATION RESULTS--------")
    print()
    print("Number of classes:", len(np.unique(class_list)))
    print("Number of mapped atoms:", n_map)
    print("Total distance between structures:", dmin[0],
          "(%f N^(4/3) + %f N)"%(dmin[1],dmin[2])) 
    print("Volume stretching factor (det(T)):", la.det(tmat))
    print("Cell volume ratio (initial cell volume)/(final cell volume):", mulA * la.det(Acell)/(mulB * la.det(Bcell)))

    # Displacements without stretching (for plotting)
    disps_total = Apos_map - Bpos

    # Displacement with stretching
    disps = Apos_map - Bposst

    # Classes of vectors and their value
    vec_classes_estimate = np.array([np.mean(disps[:,class_list==d_type], axis=1) for d_type in np.unique(class_list)])    
    
    # Run the PCA analysis if turned on
    if pca:
        print("PCA found %d classes"%PCA(disps))

    # Show and/or save interactive display of optimal result
    if savedisplay or interactive:
        print("Displaying optimal connections...")
        if interactive:
            print("(Close the display window to continue)")
        displayOptimalResult(Apos, Bpos, Bposst, disps_total, disps, class_list, vec_classes_estimate,
                                 nat, natA, natB, atoms, outdir, savedisplay, interactive)
        print()

    # Create gif if turned on
    if gif:
        makeGif(Apos, Bposst, disps, vec_classes_estimate, nat, atoms, outdir)

    # End program if a periodic cell couldn't be found earlier
    if foundcell is None:
        if switched:
            tmat = la.inv(tmat)
        return tmat, None, None, dmin

    print()
    print("-----------PERIODIC CELL-----------")
    print()

    origin = np.array([[0],[0],[0]])
    
    pos_in_struc = [None]*2
    pos_in_struc[0] = Apos_map - origin.dot(np.ones((1,np.shape(Apos_map)[1])))
    pos_in_struc[1] = Bposst - origin.dot(np.ones((1,np.shape(Bposst)[1])))
    
    dispStruc, vec_classes = makeStructures(foundcell, atoms, atom_types,
                               natB, pos_in_struc, class_list)
    
    if len(vec_classes) > len(vec_classes_estimate):
        print("Note: the number of classes found during the optimization was smaller than what was necessary to create the interface cell.")
        print()
    elif len(vec_classes) < len(vec_classes_estimate):
        print("Note: the number of classes found during the optimization was larger than what was necessary to create the interface cell.")
        print()
        
    print("Size of the transformation cell (TC):", len(dispStruc))

    print("Number of %s (%s) cells in TC:"%(A.name, fileA), abs(la.det(dispStruc.cell)/(la.det(Acell))))
    print("Number of %s (%s) cells in TC:"%(B.name, fileB), abs(la.det(dispStruc.cell)/(la.det(tmat.dot(Bcell)))))
    print()
    
    expectedA = int(np.round(abs(la.det(dispStruc.cell)/(la.det(Acell)))))*len(A)
    expectedB = int(np.round(abs(la.det(dispStruc.cell)/(la.det(tmat.dot(Bcell))))))*len(B)
    
    vacA = expectedA - len(dispStruc)
    vacB = expectedB - len(dispStruc)

    if vacA > 0 and vacB > 0:
        print("WARNING: There are missing atoms in both structures:")
        print("%s (%s) is missing %d atoms out of %d"%(A.name, fileA, vacA, expectedA))
        print("%s (%s) is missing %d atoms out of %d"%(B.name, fileB, vacB, expectedB))
    elif vacA > 0:
        print("WARNING: %d vacancies out of %d atomic positions were created in %s (%s) in order to produce a transition"%(vacA, expectedA, A.name, fileA))
    elif vacB > 0:
        print("WARNING: %d vacancies out of %d atomic positions were created in %s (%s) in order to produce a transition"%(vacB, expectedB, B.name, fileB))
    print()
        
    print("TC in %s (%s) coordinates:"%(B.name, fileB))
    printMatAndDir(la.inv(tmat).dot(dispStruc.cell), ccellB)
    print()
    print("TC in %s (%s) coordinates:"%(A.name, fileA))
    printMatAndDir(dispStruc.cell, ccellA)
    print()
        
    if interactive or savedisplay:
        if switched:
            direction = "reverse"
        else:
            direction = "direct"
        print("Displaying the Transition cell as optimized (%s)..."%(direction))
        displayTransCell(disps, dispStruc, foundcell,
                         pos_in_struc[0], vec_classes, outdir, interactive, savedisplay)
        print()

    if switched:
        dispStruc, tmat, vec_classes = switchDispStruc(dispStruc, tmat, vec_classes)
        
    return tmat, dispStruc, vec_classes, dmin
