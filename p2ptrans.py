from p2ptrans import transform as tr
import numpy as np
import numpy.linalg as la
from mpl_toolkits.mplot3d import Axes3D
import matplotlib
from matplotlib import animation
from p2ptrans import tiling as t
import pickle
import time
from pylada.crystal import Structure, primitive, gruber, read, supercell
from copy import deepcopy
import argparse
import os
import warnings

# Tolerence for structure identification
tol = 1e-5
tol_vol = 1e-3

def readOptions():

    parser = argparse.ArgumentParser()
    parser.add_argument("-I","--initial",dest="A",type=str, default='./POSCAR_A', help="Initial Structure")
    parser.add_argument("-F","--final",dest="B",type=str, default='./POSCAR_B', help="Final Structure")
    parser.add_argument("-n","--ncell",dest="ncell",type=int, default=300, help="Number of cells to tile")
    parser.add_argument("-d","--display",dest="display",action="store_true", default=False, help="Unable interactive display")
    parser.add_argument("-p","--param", dest="filename", type=str, default='./p2p.in', help="Parameter file")
    parser.add_argument("-o","--outdir",dest="outdir",type=str, default='.', help="Output directory")
    parser.add_argument("-u","--use",dest="use",action="store_true", default=False, help="Use previously calculated data")
    parser.add_argument("-s","--switch",dest="switch",action="store_true", default=False, help="Map the larger cell on the smaller cell")


    options = parser.parse_args()
    
    fileA = options.A
    fileB = options.B
    ncell = options.ncell
    filename = options.filename
    display = options.display
    outdir = options.outdir
    use = options.use
    switch = options.switch
    
    return fileA, fileB, ncell, filename, display, outdir, use, switch
    

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
                        if j == 3:
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
        for check in unique:
            if np.allclose(check, line, atol=tol):
                there = True
                break
        if not there:
            unique.append(line)
            idx.append(i)
    return (np.array(idx), np.array(unique))

def p2ptrans(fileA, fileB, ncell, filename, display, outdir, use, switch):

    if not display:
        matplotlib.use('Agg')

    import matplotlib.pyplot as plt

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    random = False
    
    A = read.poscar(fileA)
    B = read.poscar(fileB)
    
    print("POSCARS")
    print(A)
    print(B)

    mul = lcm(len(A),len(B))
    mulA = mul//len(A)
    mulB = mul//len(B)

    Acell = A.cell*float(A.scale)
    Bcell = B.cell*float(B.scale)

    if (abs(mulA*la.det(Acell)) < abs(mulB*la.det(Bcell))) != switch: # (is switched?) != switch
        tmp = deepcopy(B)
        tmpmul = mulB
        tmpcell = Bcell
        B = deepcopy(A)
        mulB = mulA
        Bcell = Acell
        A = tmp
        mulA = tmpmul
        Acell = tmpcell

    print(mulA, Acell, la.det(Acell))
    print(mulB, Bcell, la.det(Bcell))
        
    
    # Plotting the cell vectors of A and B
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.quiver(np.ones(3), np.ones(3), np.ones(3), Bcell[0,:], Bcell[1,:], Bcell[2,:])
    ax.quiver(-np.ones(3), -np.ones(3), -np.ones(3), Acell[0,:], Acell[1,:], Acell[2,:])
    ax.set_xlim([-5, 5])
    ax.set_ylim([-5, 5])
    ax.set_zlim([-5, 5])
    fig.savefig(outdir+'/CellVectors.svg')
    
    tmat = np.eye(3)
    found = False
    rep = 0
    while (not found and rep < 2):
        rep += 1
        
        # Adds atoms to A and B (for cell with different types of atoms)
        Apos = []
        atom_types = np.array([], np.str)
        atomsA = np.array([], np.int)
        for a in A:
            if any(atom_types == a.type):
                idx = np.where(atom_types == a.type)[0][0]
                Apos[idx] = np.concatenate((Apos[idx], t.sphere(Acell, mulA*ncell, a.pos*float(A.scale))), axis = 1) 

                # Order the atoms in terms of distance
                Apos[idx] = Apos[idx][:,np.argsort(la.norm(Apos[idx],axis=0))] 
                atomsA[idx] += 1
            else:
                Apos.append(t.sphere(Acell, mulA*ncell, a.pos*float(A.scale)))
                atom_types = np.append(atom_types, a.type)
                atomsA = np.append(atomsA,1)
        
        Apos = np.concatenate(Apos, axis=1)
        
        # Temporarly stretching Bcell, for tiling
        Bcell = tmat.dot(Bcell)

        Bpos = [None]*len(atom_types)
        atomsB = np.zeros(len(atom_types), np.int)
        for a in B:
            idx = np.where(atom_types == a.type)[0][0]
            if atomsB[idx] == 0:
                Bpos[idx] = t.sphere(Bcell, mulB*ncell, tmat.dot(a.pos)*float(B.scale))
            else:
                Bpos[idx] = np.concatenate((Bpos[idx], t.sphere(Bcell, mulB*ncell, tmat.dot(a.pos)*float(B.scale))), axis = 1) 
                # Order the atoms in terms of distance
                Bpos[idx] = Bpos[idx][:,np.argsort(la.norm(Bpos[idx],axis=0))]
            atomsB[idx] += 1
        
        Bpos = np.concatenate(Bpos, axis=1)

        Bpos = np.linalg.inv(tmat).dot(Bpos)
        Bcell = np.linalg.inv(tmat).dot(Bcell)
            
        assert all(mulA*atomsA == mulB*atomsB)
        atoms = mulA*atomsA
        
        if not use:
            Apos = np.asfortranarray(Apos)
            Bpos = np.asfortranarray(Bpos) 
            t_time = time.time()
            Apos_map, Bpos, Bposst, n_map, natA, class_list, tmat, dmin = tr.fastoptimization(Apos, Bpos, Acell, la.inv(Acell), mulA * la.det(Acell)/(mulB * la.det(Bcell)), atoms, filename)
            t_time = time.time() - t_time
            Bpos = np.asanyarray(Bpos)
            Apos = np.asanyarray(Apos)
        
            print("Mapping time:", t_time)
        
            pickle.dump((Apos_map, Bpos, Bposst, n_map, natA, class_list, tmat, dmin), open(outdir+"/fastoptimization.dat","wb"))
        
        else:
            print("Using data from "+outdir+"/fastoptimization.dat")
            Apos = np.asfortranarray(Apos)
            Bpos = np.asfortranarray(Bpos) 
            tr.center(Apos)
            tr.center(Bpos)
            Apos_map, Bpos, Bposst, n_map, natA , class_list, tmat, dmin = pickle.load(open(outdir+"/fastoptimization.dat","rb"))
            Bpos = np.asanyarray(Bpos)
            Apos = np.asanyarray(Apos)
        
        print("Total distance between structures:", dmin)
        
        class_list = class_list[:n_map]-1
        
        Bpos = Bpos[:,:n_map]
        Bposst = Bposst[:,:n_map]
        Apos_map = Apos_map[:,:n_map]
        
        natB = n_map // np.sum(atoms)
        nat_map = n_map // np.sum(atoms)
        nat = np.shape(Apos)[1] // np.sum(atoms)
        print("NAT", nat)
        print("N_MAP",n_map)


        try:
            foundcell, origin = find_cell(class_list, Bposst)
            if abs(la.det(tmat) - mulA * la.det(Acell)/(mulB * la.det(Bcell))) > tol_vol:
                found = False
                print("The volume factor is wrong.")
                print("_____RESTARTING_____")
        except RuntimeError:
            print("Could not find periodic cell")
            print("_____RESTARTING_____")
            found = False

    # Plotting the Apos and Bpos overlayed
    fig = plt.figure(22)
    ax = fig.add_subplot(111, projection='3d')
    ax.view_init(90,0) # TMP
    #ax.scatter(Apos.T[:,0],Apos.T[:,1])
    num_tot = 0

    for i,num in enumerate(atoms):
        ax.scatter(Apos.T[num_tot*nat:num_tot*nat+natA*num+1,0],Apos.T[num_tot*nat:num_tot*nat+natA*num+1,1],Apos.T[num_tot*nat:num_tot*nat+natA*num+1,2], c="C%d"%(2*i))
        ax.scatter(Apos.T[num_tot*nat+natA*num:(num_tot + num)*nat+1,0],Apos.T[num_tot*nat+natA*num:(num_tot + num)*nat+1,1],Apos.T[num_tot*nat+natA*num:(num_tot + num)*nat+1,2], c="C%d"%(2*i), alpha=0.1)
        ax.scatter(Bpos.T[natB*num*i:natB*num*(i+1),0],Bpos.T[natB*num*i:natB*num*(i+1),1], Bpos.T[natB*num*i:natB*num*(i+1),2], c="C%d"%(2*i+1))
        num_tot = num_tot + num
    
    centerofmassA = np.mean(Apos,axis=1)
    centerofmassB = np.mean(Bpos,axis=1)

    print("Center of mass A", centerofmassA)
    print("Center of mass B", centerofmassB)

    ax.scatter(centerofmassA[0], centerofmassA[1], centerofmassA[2], s=60, c='red')
    ax.scatter(centerofmassB[0], centerofmassB[1], centerofmassB[2], s=60, c='green')

    maxXAxis = np.max([Apos.max(), Bpos.max()]) + 1
    ax.set_xlim([-maxXAxis, maxXAxis])
    ax.set_ylim([-maxXAxis, maxXAxis])
    ax.set_zlim([-maxXAxis, maxXAxis])
    ax.set_aspect('equal')
    
    
    # Displacements without stretching (for plotting)
    print("SHAPE", np.shape(Apos_map), np.shape(Bpos)) 
    disps = Apos_map - Bpos
    
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
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

    fig =plt.figure()
    ax = Axes3D(fig)
    centerofmassA = np.mean(Apos,axis=1)
    centerofmassB = np.mean(Bpos,axis=1)

    print("Center of mass A", centerofmassA)
    print("Center of mass B", centerofmassB)

    ax.scatter(centerofmassA[0], centerofmassA[1], centerofmassA[2], s=60, c='red')
    ax.scatter(centerofmassB[0], centerofmassB[1], centerofmassB[2], s=60, c='green')
    ax.quiver(centerofmassA[0], centerofmassA[1], centerofmassA[2], Acell[0,:], Acell[1,:], Acell[2,:])
    ax.quiver(centerofmassB[0], centerofmassB[1], centerofmassB[2], tmat.dot(Bcell[0,:]), tmat.dot(Bcell[1,:]), tmat.dot(Bcell[2,:]))

    maxXAxis = np.max([Apos.max(), Bpos.max()]) + 1
    ax.set_xlim([-maxXAxis, maxXAxis])
    ax.set_ylim([-maxXAxis, maxXAxis])
    ax.set_zlim([-maxXAxis, maxXAxis])
    ax.set_aspect('equal')

    
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
        num_tot = 0
        for i,num in enumerate(atoms):
            ax.scatter(Apos.T[num_tot*nat:num_tot*nat+natA*num+1,0],Apos.T[num_tot*nat:num_tot*nat+natA*num+1,1],Apos.T[num_tot*nat:num_tot*nat+natA*num+1,2], c="C%d"%(2*i))
            ax.scatter(Apos.T[num_tot*nat+natA*num:(num_tot + num)*nat+1,0],Apos.T[num_tot*nat+natA*num:(num_tot + num)*nat+1,1],Apos.T[num_tot*nat+natA*num:(num_tot + num)*nat+1,2], c="C%d"%(2*i), alpha=0.1)
            ax.scatter(Bposst.T[natB*num*i:natB*num*(i+1),0],Bposst.T[natB*num*i:natB*num*(i+1),1], Bposst.T[natB*num*i:natB*num*(i+1),2], c="C%d"%(2*i+1))
            num_tot = num_tot + num

        for i in range(len(vec_classes)):
            disps_class = disps[:,class_list==i]
            Bposst_class = Bposst[:,class_list==i]
            ndisps = np.shape(disps_class)[1]
            ax.quiver(Bposst_class.T[:,0], Bposst_class.T[:,1], Bposst_class.T[:,2], disps_class.T[:,0], disps_class.T[:,1], disps_class.T[:,2], color="C%d"%(i%10))
        fig.savefig(outdir+'/DispLattice_stretched.svg')
        return fig,

    if not display:
        anim = animation.FuncAnimation(fig, animate, init_func=init_disps,
                                       frames=490, interval=30)
        anim.save(outdir+'/Crystal+Disps.gif', fps=30, codec='gif')
    else:
        init_disps()
    
    
    # Stretching Matrix
    stMat = la.inv(tr.canonicalize(Bcell)).dot(tr.canonicalize(tmat.dot(Bcell)))
    
    # Rotation Matrix
    rtMat = tmat.dot(Bcell).dot(la.inv(stMat)).dot(la.inv(Bcell))
    
    print("Stretching Matrix:")
    print(stMat)
    
    print("Rotation Matrix:")
    print(rtMat)

    # Just the displacements
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    maxXAxis = disps.max()
    ax.set_xlim([-maxXAxis, maxXAxis])
    ax.set_ylim([-maxXAxis, maxXAxis])
    ax.set_zlim([-maxXAxis, maxXAxis])
    ax.set_aspect('equal')
    for i in range(len(vec_classes)):
        disps_class = disps[:,class_list==i]
        ndisps = np.shape(disps_class)[1]
        ax.quiver(np.zeros((1,ndisps)), np.zeros((1,ndisps)), np.zeros((1,ndisps)), disps_class.T[:,0], disps_class.T[:,1], disps_class.T[:,2], color="C%d"%(i%10))
    fig.savefig(outdir+'/DispOverlayed.svg')

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
    print("Volume stretching factor:", la.det(stMat), la.det(tmat))
    print("Cell volume ratio (should be exactly the same):", mulA * la.det(Acell)/(mulB * la.det(Bcell)))
        
    print("Showing")
    
    plt.show()

    if not found:
        raise RuntimeError("Could not find good displacement cell. Increase system size")

    cell = foundcell
    
    pos_in_struc = Bposst - origin.dot(np.ones((1,np.shape(Bposst)[1])))

    # Make a pylada structure
    cell_coord = np.mod(la.inv(cell).dot(pos_in_struc)+tol,1)-tol
    Struc = Structure(cell)
    print(zip(uniqueclose(cell_coord, tol)))
    print(uniqueclose(cell_coord, tol))
    for i, disp in zip(*uniqueclose(cell_coord, tol)):
        Struc.add_atom(*(tuple(cell.dot(disp))+(str(class_list[i]),)))

    if la.det(cell) < 0:
        cell[:,2] = -cell[:,2] 

    Struc = supercell(Struc, cell)

    # Finds a squarer cell
    cell = gruber(cell)

    # Makes sure it is the primitive cell 
    Struc = primitive(Struc, tolerance = tol)
    
    
    print("Is it a supercell?")
    print(la.inv(Acell)*Struc.cell)
    print(la.inv(Bcell)*Struc.cell)
    

    # Total displacement per unit volume a as metric
    Total_disp = 0 
    for disp in Struc:
        Total_disp += la.norm(vec_classes[int(disp.type)])
    
    Total_disp = Total_disp / la.det(Struc.cell)
    
    cell = Struc.cell
    
    print("Displacement Structure")
    print(Struc)
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
            ax.quiver(disp.pos[0], disp.pos[1], disp.pos[2], vec_classes[int(disp.type)][0],vec_classes[int(disp.type)][1], vec_classes[int(disp.type)][2], color="C%d"%(i%10))
            ax.scatter(disp.pos[0], disp.pos[1], disp.pos[2], alpha = 0.5, s=10, color="C%d"%(i%10))
        ax.quiver(np.zeros(3), np.zeros(3), np.zeros(3), cell[0,:], cell[1,:], cell[2,:], color = "red", alpha = 0.3)
        maxXAxis = abs(cell).max() + 1
        ax.set_xlim([-maxXAxis, maxXAxis])
        ax.set_ylim([-maxXAxis, maxXAxis])
        ax.set_zlim([-maxXAxis, maxXAxis])
        ax.set_aspect('equal')
        fig.savefig(outdir+'/Displacement_structure.svg')
        return fig,
    
    if not display:
        anim = animation.FuncAnimation(fig, animate, init_func=init_struc,
                                       frames=490, interval=30)
        anim.save(outdir+'/DispStruc.gif', fps=30, codec='gif')
    else:
        init_struc()
        
    plt.show()
    
    plt.close('All')

if __name__=='__main__':
    p2ptrans(*readOptions())
    


