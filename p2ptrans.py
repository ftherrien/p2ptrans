from p2ptrans import transform as tr
import numpy as np
import numpy.linalg as la
from mpl_toolkits.mplot3d import Axes3D
import matplotlib
from matplotlib import animation
from p2ptrans import tiling as t
import pickle
import time
from pylada.crystal import Structure, primitive, gruber, read, write, supercell, space_group
from copy import deepcopy
import argparse
import os
import warnings

colorlist=['#929591', 'r', 'k','b','#06470c','#ceb301', '#9e0168', '#26f7fd', '#f97306', '#c20078']

pca = False

# Tolerence for structure identification
tol = 1e-5
tol_vol = 2*1e-3
tol_uvw = 1e-6
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
    parser.add_argument("-r","--primitive",dest="prim",action="store_true", default=False, help="Finds the primitive cell at the beginning") #TMP


    options = parser.parse_args()
    
    fileA = options.A
    fileB = options.B
    ncell = options.ncell
    filename = options.filename
    display = options.display
    outdir = options.outdir
    use = options.use
    switch = options.switch
    prim = options.prim
    
    return fileA, fileB, ncell, filename, display, outdir, use, switch, prim
    
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

    return cell, newcell
                

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
                            repeat = True
                            while repeat:
                                for l, cell in enumerate(cell_list):
                                    prev_cell = cell
                                    prev_newcell = newcell
                                    cell_list[l], newcell = find_supercell(cell, newcell, tol)
                                else:
                                    repeat = False
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

def lcm3(x, y, z):
   """This function takes two
   integers and returns the L.C.M."""

   gcd3 = gcd(x,gcd(y,z))
   lcm = (x*y*z)/ (gcd3**2 * gcd(x/gcd3,y/gcd3) * gcd(x/gcd3,z/gcd3) * gcd(z/gcd3,y/gcd3)) 
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
                class_list[i] = j
                classified = True
                break
        if not classified:
            vec_classes.append(disps[:,i:i+1])
            class_list[i] = len(vec_classes) - 1
            
    for i, elem in enumerate(vec_classes):
        vec_classes[i] = np.mean(elem, axis=1)
        
    return class_list, vec_classes

def dir2angles(plane):
    angles=np.zeros(2)
    a0 = np.arccos(plane[2])
    angles[0] = np.pi/2 - a0
    angles[1] = np.arctan2(plane[1], plane[0])
    return angles*180/np.pi


def p2ptrans(fileA, fileB, ncell, filename, display, outdir, use, switch, prim):

    on_top = None

    if not display:
        matplotlib.use('Agg')

    import matplotlib.pyplot as plt

    import matplotlib.gridspec as gridspec
    from matplotlib.patches import Rectangle

    def add_panel(fig,g,angles, state, p, label):
        color_list=[]
        ax = fig.add_subplot(g, projection='3d', proj_type = 'ortho')
        ax.view_init(*angles)
        maxXAxis = np.max([Apos.max(), Bpos.max()]) + 1
        ax.set_xlim([-maxXAxis, maxXAxis])
        ax.set_ylim([-maxXAxis, maxXAxis])
        ax.set_zlim([-maxXAxis, maxXAxis])
        ax.set_aspect('equal')
        toplot = Apos_map[:,idx[p]] - disps_total[:,idx[p]]*1/(np.exp(10*(state-n_states/2)/n_states) + 1)
        color_to_plot = color_array[idx[p]]
        idxx = np.argsort(toplot[2,:])
        toplot = toplot[:,idxx]
        color_to_plot = color_to_plot[idxx]
        ax.set_axis_off()
        ax.dist = 2
        axlims = [a for b in ax.get_position().get_points() for a in b]
        rec = Rectangle((axlims[0],axlims[1]),(axlims[2]-axlims[0]),(axlims[3]-axlims[1]), transform = fig.transFigure, fill=False,lw=1, color="k")
        fig.patches.append(rec)
        for i,point in enumerate(toplot.T):
            if any(color_list==color_to_plot[i]) or not label or p:
                ax.scatter(*point, c=colorlist[color_to_plot[i]], s=(2-color_to_plot[i])*40, depthshade=False)
            else:
                color_list.append(color_to_plot[i])
                ax.scatter(*point, c=colorlist[color_to_plot[i]], s=(2-color_to_plot[i])*40, depthshade=False, label=atom_types[color_to_plot[i]])

    def all_panels(fig, gs, state, label):
        def func(start,end):
            return start + (end-start)*1/(np.exp(10*(state-n_states/2)/n_states) + 1)
        for p,pl in enumerate(planes):
            plane = dispStruc.cell.dot(pl)
            plane = plane/la.norm(plane)
            angles = dir2angles(plane)
            if p ==0:
                add_panel(fig,gs[0:2,0:2],(func(angles[0],angles[0]),func(angles[1],angles[1])),state, p, label)
            elif p == 1:
                add_panel(fig,gs[0,2],(func(angles[0],angles[0]),func(angles[1],angles[1])),state, p, label)
            else:
                add_panel(fig,gs[1,2],(func(angles[0],angles[0]),func(angles[1],angles[1])),state, p, label)
        fig.legend()

    def make_fig(state,n_states):
        fig = plt.figure(figsize=[12.8,7.2])
        gs = gridspec.GridSpec(2, 3)
        gs.update(wspace=0.01, hspace=0.01)
        all_panels(fig,gs, state, True)
    
    def make_anim(n_states):
        fig = plt.figure(figsize=[12.8,7.2])
        gs = gridspec.GridSpec(2, 3)
        gs.update(wspace=0.01, hspace=0.01)
        def animate_trans(state):
            all_panels(fig,gs, state, state==1)

        animation.verbose.set_level('debug')
    
        plt.rcParams['animation.ffmpeg_path'] = '/home/felixt/bin/ffmpeg'
        # Writer = animation.writers['ffmpeg']
        writer = animation.FFMpegWriter(fps=30,codec='prores', extra_args=['-loglevel', 'verbose','-f','mov'])
    
        anim = animation.FuncAnimation(fig, animate_trans,
                                   frames=n_states, interval=1)
        anim.save(outdir + '/Trans.mov', writer=writer)

    def PCA(disps):
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

        # # I don't really know how to interpret the vectors...
        # vecs = disps.dot(eigvec[:,:n_class]/np.ones((n,1)).dot(np.reshape(np.sum(eigvec[:,:n_class], axis=0),(1,n_class))))

        # vecs = vecs/np.ones((3,1)).dot(np.reshape(la.norm(vecs, axis=0),(1,n_class)))
        
        # fig = plt.figure()
        # ax = fig.add_subplot(111, projection='3d')
        # maxXAxis = vecs.max()
        # ax.set_xlim([-maxXAxis, maxXAxis])
        # ax.set_ylim([-maxXAxis, maxXAxis])
        # ax.set_zlim([-maxXAxis, maxXAxis])
        # ax.set_aspect('equal')
        # print(np.shape(vecs), n_class)
        # ax.quiver(np.zeros((1,n_class)), np.zeros((1,n_class)), np.zeros((1,n_class)),  vecs[0,:], vecs[1,:], vecs[2,:])

        plt.figure()
        plt.plot(eigval,".")
        
        plt.figure()
        plt.semilogy(eigval,".")

        plt.figure()
        plt.plot(logdiffs,".-")
        
        return n_class

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    random = False
    
    A = read.poscar(fileA)
    B = read.poscar(fileB)

    if prim:
        A = primitive(A)
        B = primitive(B)
    
    print("POSCARS")
    print(A)
    print(B)

    mul = lcm(len(A),len(B))
    mulA = mul//len(A)
    mulB = mul//len(B)

    Acell = A.cell*float(A.scale)
    Bcell = B.cell*float(B.scale)

    if (abs(mulA*la.det(Acell)) < abs(mulB*la.det(Bcell))) != switch: # (is switched?) != switch
        print("Transition from %s to %s"%(fileA, fileB))
        tmp = deepcopy(B)
        tmpmul = mulB
        tmpcell = Bcell
        B = deepcopy(A)
        mulB = mulA
        Bcell = Acell
        A = tmp
        mulA = tmpmul
        Acell = tmpcell
    else:
        print("Transition from %s to %s"%(fileB, fileA))
    
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
    while (not found and rep < 1):
        rep += 1
        found = True
        
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

        Bpos = la.inv(tmat).dot(Bpos)
        Bcell = la.inv(tmat).dot(Bcell)
            
        assert all(mulA*atomsA == mulB*atomsB)
        atoms = mulA*atomsA
        
        if not use:
            Apos = np.asfortranarray(Apos)
            Bpos = np.asfortranarray(Bpos)
            tr.center(Bpos)
            oldBpos = np.asanyarray(deepcopy(Bpos))
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
            oldBpos = np.asanyarray(deepcopy(Bpos))
            Apos_map, Bpos, Bposst, n_map, natA , class_list, tmat, dmin = pickle.load(open(outdir+"/fastoptimization.dat","rb"))
            Bpos = np.asanyarray(Bpos)
            Apos = np.asanyarray(Apos)
        
        print("Total distance between structures:", dmin)
        
        class_list = class_list[:n_map]-1

        vec = Bposst[:,0:1] - tmat.dot(oldBpos)[:,0:1]
        vecv = vec.dot(np.ones((1,n_map)))

        Bpos = Bpos[:,:n_map]
        Bposst = Bposst[:,:n_map]
        Apos_map = Apos_map[:,:n_map]
        
        natB = n_map // np.sum(atoms)
        nat_map = n_map // np.sum(atoms)
        nat = np.shape(Apos)[1] // np.sum(atoms)
        print("NAT", nat)
        print("N_MAP",n_map)

        # disps = Apos_map - Bposst #TMP
        # class_list, vec_classes = classify(disps, tol = 1.e-2) #TMP

        foundcell, origin = find_cell(class_list, Bposst)

        try:
            foundcell, origin = find_cell(class_list, Bposst)
            if abs(abs(la.det(tmat)) - abs(mulA * la.det(Acell)/(mulB * la.det(Bcell)))) > tol_vol:
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
    ax.view_init(0,0) # TMP
    # ax.scatter(Apos.T[:,0],Apos.T[:,1])
    num_tot = 0


    # TMP
    # Stretching Matrix
    stMat = la.inv(tr.canonicalize(Bcell)).dot(tr.canonicalize(tmat.dot(Bcell)))
    
    stMat2 = Bcell.dot(stMat).dot(la.inv(Bcell))

    # Rotation Matrix
    rtMat = tmat.dot(Bcell).dot(la.inv(stMat)).dot(la.inv(Bcell))

    for i,num in enumerate(atoms):
        ax.scatter(Apos.T[num_tot*nat:num_tot*nat+natA*num+1,0],Apos.T[num_tot*nat:num_tot*nat+natA*num+1,1],Apos.T[num_tot*nat:num_tot*nat+natA*num+1,2], c=colorlist[2*i])
        ax.scatter(Apos.T[num_tot*nat+natA*num:(num_tot + num)*nat+1,0],Apos.T[num_tot*nat+natA*num:(num_tot + num)*nat+1,1],Apos.T[num_tot*nat+natA*num:(num_tot + num)*nat+1,2], c=colorlist[2*i], alpha=0.1)
        ax.scatter(Bpos.T[natB*num_tot:natB*(num_tot+num),0],Bpos.T[natB*num_tot:natB*(num_tot+num),1], Bpos.T[natB*num_tot:natB*(num_tot+num),2], c=colorlist[2*i+1])
        num_tot = num_tot + num
    
    centerofmassA = np.mean(Apos,axis=1)
    centerofmassB = np.mean(Bpos,axis=1)

    # TMP getting the rotation
    oldV1 = oldBpos[:,0]
    newV1 = (Bpos[:,0] - centerofmassA.reshape((3,1)))
    
    u = np.cross(oldV1,newV1)/(la.norm(oldV1)*la.norm(newV1))
    theta = np.arccos(oldV1.dot(newV1)/(la.norm(oldV1)*la.norm(newV1)))

    P = u.dot(u.T)
    Q = np.array([[0.0,-u[2,0],u[1,0]],[u[2,0],0.0,-u[0,0]],[-u[1,0],u[0,0],0.0]])

    R = P + (np.eye(3) - P)*np.cos(theta) + Q*np.sin(theta)

    print("RRRR", R)
    print("Bcell", R.dot(Bcell))
    
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
    disps_total = Apos_map - Bpos
    
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    ax.quiver(Bpos.T[:,0], Bpos.T[:,1], Bpos.T[:,2], disps_total.T[:,0], disps_total.T[:,1], disps_total.T[:,2])
    maxXAxis = np.max([Apos.max(), Bpos.max()]) + 1
    ax.set_xlim([-maxXAxis, maxXAxis])
    ax.set_ylim([-maxXAxis, maxXAxis])
    ax.set_zlim([-maxXAxis, maxXAxis])
    ax.set_aspect('equal')
    fig.savefig(outdir+'/DispLattice.svg')
    
    # Displacement with stretching
    disps = Apos_map - Bposst
    vec_classes = np.array([np.mean(disps[:,class_list==d_type], axis=1) for d_type in np.unique(class_list)])

    if pca:
        print("PCA found %d classes"%PCA(disps))

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
    ax.view_init(0,0) # TMP
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
            ax.scatter(Apos.T[num_tot*nat:num_tot*nat+natA*num+1,0],Apos.T[num_tot*nat:num_tot*nat+natA*num+1,1],Apos.T[num_tot*nat:num_tot*nat+natA*num+1,2], c=colorlist[2*i])
            ax.scatter(Apos.T[num_tot*nat+natA*num:(num_tot + num)*nat+1,0],Apos.T[num_tot*nat+natA*num:(num_tot + num)*nat+1,1],Apos.T[num_tot*nat+natA*num:(num_tot + num)*nat+1,2], c=colorlist[2*i], alpha=0.1)
            ax.scatter(Bposst.T[natB*num_tot:natB*(num_tot+num),0],Bposst.T[natB*num_tot:natB*(num_tot+num),1], Bposst.T[natB*num_tot:natB*(num_tot+num),2], c=colorlist[2*i+1])
            num_tot = num_tot + num

        for i in range(len(vec_classes)):
            disps_class = disps[:,class_list==i]
            Bposst_class = Bposst[:,class_list==i]
            ndisps = np.shape(disps_class)[1]
            ax.quiver(Bposst_class.T[:,0], Bposst_class.T[:,1], Bposst_class.T[:,2], disps_class.T[:,0], disps_class.T[:,1], disps_class.T[:,2], color=colorlist[i%10])
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
        ax.quiver(np.zeros((1,ndisps)), np.zeros((1,ndisps)), np.zeros((1,ndisps)), disps_class.T[:,0], disps_class.T[:,1], disps_class.T[:,2], color=colorlist[i%10])
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
    
    # plt.show()

    if not found:
        raise RuntimeError("Could not find good displacement cell. Increase system size")

    cell = foundcell
    
    pos_in_struc = Bposst - origin.dot(np.ones((1,np.shape(Bposst)[1])))
    pos_in_strucB = Bpos - origin.dot(np.ones((1,np.shape(Bposst)[1])))
    pos_in_strucA = Apos_map - origin.dot(np.ones((1,np.shape(Bposst)[1])))  

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
            i = np.argmin(la.norm(np.array([pos_in_struc[:,i] for i in idx]),axis=0))
            incell.append((i,pos_in_struc[:,i]))
            print("Could not find all disp. in first cell, certain displacement could be broken")
        
    for i, disp in incell:
        dispStruc.add_atom(*(tuple(disp)+(str(class_list[i]),)))
        stinitStruc.add_atom(*(tuple(disp)+(whattype(i, natB),)))

    if la.det(cell) < 0:
       cell[:,2] = -cell[:,2] 

    # Finds a squarer cell
    # cell = gruber(cell)

    dispStruc = supercell(dispStruc, cell)

    # Makes sure it is the primitive cell 
    dispStruc = primitive(dispStruc, tolerance = tol)
    stinitStruc = supercell(stinitStruc, dispStruc.cell)

    finalStruc = Structure(dispStruc.cell)
    for i,a in enumerate(dispStruc):
        finalStruc.add_atom(*(a.pos+vec_classes[int(a.type)]),stinitStruc[i].type)

    print("Is it a supercell?")
    print(la.inv(Acell)*dispStruc.cell)
    print(la.inv(Bcell)*dispStruc.cell)

    print("Number of A cell in dispCell:", la.det(dispStruc.cell)/(mulA*la.det(Acell)))
    print("Number of B cell in dispCell:", la.det(dispStruc.cell)/(mulB*la.det(tmat.dot(Bcell))))

    # Produce transition

    n_steps = 100

    os.makedirs(outdir+"/TransPOSCARS", exist_ok=True)

    spgList = []
    for i in range(n_steps):
        curMat = (la.inv(tmat)-np.eye(3))*i/n_steps + np.eye(3)
        curStruc = Structure(finalStruc.cell.dot(curMat))
        for j,a in enumerate(dispStruc):
            curDisp = vec_classes[int(a.type)]*i/n_steps
            curPos = (finalStruc[j].pos - curDisp).dot(curMat)
            curStruc.add_atom(*(curPos),finalStruc[j].type)
        write.poscar(curStruc, vasp5=True, file=outdir+"Trans_%d.POSCAR")
        spgList.append(space_group(curStruc))

    print(spgList)


    # Total displacement per unit volume a as metric
    Total_disp = 0
    for disp in dispStruc:
        Total_disp += la.norm(vec_classes[int(disp.type)])
    
    Total_disp = Total_disp / la.det(dispStruc.cell)
    
    cell = dispStruc.cell
    
    print("Displacement Structure")
    print(dispStruc)
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
        for i,disp in enumerate(dispStruc):
            ax.quiver(disp.pos[0], disp.pos[1], disp.pos[2], vec_classes[int(disp.type)][0],vec_classes[int(disp.type)][1], vec_classes[int(disp.type)][2], color=colorlist[i%10])
            ax.scatter(disp.pos[0], disp.pos[1], disp.pos[2], alpha = 0.5, s=10, color=colorlist[i%10])
            ax.scatter(finalStruc[i].pos[0], finalStruc[i].pos[1], finalStruc[i].pos[2], alpha = 1, s=10, color=colorlist[i%10])
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

    # Growth directions

    eigval, eigvec = la.eig(stMat2)

    stretch_dir = la.inv(Bcell).dot(eigvec)

    # stretch_dir = stretch_dir/np.ones((3,1)).dot(np.reshape(np.amax(abs(stretch_dir), axis=0),(1,3)))

    # # Failed attempt at findin uvw
    # for vec in stretch_dir.T:
    #     mul = np.zeros(3)
    #     for i in range(3):
    #         elema = np.round(1/tol_uvw*vec[i])
    #         print(elema)
    #         mul[i] = lcm(elema,1/tol_uvw)/elema
    #     print(mul)
    #     mul = lcm3(*mul)
    #     print(mul)
    #     print(mul*vec)

    # closest uvw for -10 to 10
    min_dist = np.ones(3)*1000
    min_uvw = np.zeros((3,3))
    for l,vec in enumerate(stretch_dir.T):
        for i in np.arange(-10,10):
            for j in np.arange(-10,10):
                for k in np.arange(-10,10):
                    if [i,j,k] != [0,0,0]:
                        unit_vec = Bcell.dot(np.array([i,j,k]))
                        unit_vec = unit_vec/la.norm(unit_vec)
                        cur_dist = la.norm(unit_vec-Bcell.dot(vec))
                        if cur_dist < min_dist[l]:
                            min_dist[l] = cur_dist
                            min_uvw[:,l] = np.array([i,j,k]).T
        gcd3 = gcd(min_uvw[0,l],gcd(min_uvw[1,l], min_uvw[2,l]))
        min_uvw[:,l] = min_uvw[:,l]/gcd3

    print("Exact directions of stretching in coord. of B")
    print(stretch_dir)
    print("Approx. UVW directions of stretching")
    print(min_uvw)
    print("distance between UVW and actual direction (in card. coord)")
    print(min_dist)
    print("Amount of stretching")
    print(eigval)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(Bpos.T[:,0], Bpos.T[:,1], Bpos.T[:,2], s=10, color = "C0")
    ax.quiver(np.zeros(3), np.zeros(3), np.zeros(3), eigvec[0,:], eigvec[1,:], eigvec[2,:], color = "red")
    maxXAxis = pos_in_struc.max() + 1
    ax.set_xlim([-maxXAxis, maxXAxis])
    ax.set_ylim([-maxXAxis, maxXAxis])
    ax.set_zlim([-maxXAxis, maxXAxis])
    ax.set_aspect('equal')

    color_array = np.array([a for i,num in enumerate(atoms) for a in [i]*(natB*num)])

    planes = [[1,0,0], [0,1,0], [0,0,1]]

    idx = [None]*len(planes)

    for l,pl in enumerate(planes):

        plane = dispStruc.cell.dot(pl)
        
        tickness = la.norm(plane)
        
        plane = plane/tickness
        
        print("THICKNESS", tickness)
        
        lattices = []
        
        # Cut along the current plane
        for i in np.arange(-3,3):
            for j in np.arange(-3,3):
                for k in np.arange(-3,3):
                    pos = dispStruc.cell.dot(np.array([i,j,k]))
                    if plane.dot(pos) < 0 and plane.dot(pos) > - tickness:
                        lattices.append([i,j,k])
        
                        
        dispSupercell = dispStruc.cell.dot(np.diag(np.max(abs(np.array(lattices)), axis = 0)))
        
        superDispStruc = supercell(dispStruc, dispSupercell)
        superstInitStruc = supercell(stinitStruc, dispSupercell)
        superFinalStruc = supercell(finalStruc, dispSupercell)
        
        write.poscar(stinitStruc, vasp5=True, file=outdir+"InitStruc_POSCAR_%s"%("-".join([str(p) for p in pl])))
        write.poscar(superstInitStruc, vasp5=True, file=outdir+"SuperstInitStruc_POSCAR_%s"%("-".join([str(p) for p in pl])))
        write.poscar(superFinalStruc, vasp5=True, file=outdir+"SuperFinalStruc_POSCAR_%s"%("-".join([str(p) for p in pl])))
        
        centerCell = np.mean([a.pos for a in superFinalStruc], axis=0) 
        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d') 
        for i,disp in enumerate(superDispStruc):
            color = np.where(atom_types == superFinalStruc[i].type)[0][0]
            ax.quiver(disp.pos[0] - centerCell[0], disp.pos[1] - centerCell[1], disp.pos[2] - centerCell[2], vec_classes[int(disp.type)][0],vec_classes[int(disp.type)][1], vec_classes[int(disp.type)][2], color=colorlist[2*color])
            ax.scatter(disp.pos[0] - centerCell[0], disp.pos[1] - centerCell[1], disp.pos[2] - centerCell[2], alpha = 0.5, s=10, color=colorlist[2*color])
            ax.scatter(superFinalStruc[i].pos[0] - centerCell[0], superFinalStruc[i].pos[1] - centerCell[1], superFinalStruc[i].pos[2] - centerCell[2], alpha = 0.5, s=10, color=colorlist[2*color + 1])
        ax.quiver(np.zeros(3), np.zeros(3), np.zeros(3), cell[0,:], cell[1,:], cell[2,:], color = "red", alpha = 0.3)
        
        ax.quiver(np.zeros(3), np.zeros(3), np.zeros(3), rtMat.dot(stretch_dir)[0,:]*10, rtMat.dot(stretch_dir)[1,:]*10, rtMat.dot(stretch_dir)[2,:]*10, color="blue")
        
        maxXAxis = abs(superDispStruc.cell).max() + 1
                     
        ax.set_xlim([-maxXAxis, maxXAxis])
        ax.set_ylim([-maxXAxis, maxXAxis])
        ax.set_zlim([-maxXAxis, maxXAxis])
        ax.set_aspect('equal')
        ax.view_init(*dir2angles(plane))
        fig.savefig(outdir+'/Displacement_structure_from_above.svg')
            

        # Animation

        if on_top == None:
            on_top = A[0].type 
        
        tickness = tickness*2

        offset = tickness/2

        tol_t = 1e-3

        plane_origin = Apos_map[:,np.argmax(plane.dot(Apos_map))]
        for i in range(n_map):
            if plane.dot(Apos_map[:,i]) + tol_t > offset and plane.dot(Apos_map[:,i]) < plane.dot(plane_origin) and atom_types[color_array[i]] == on_top:
                plane_origin = Apos_map[:,i]

        idx[l] = []
        for i in range(n_map):
            if plane.dot(Apos_map[:,i]) < plane.dot(plane_origin) + tol_t and plane.dot(Apos_map[:,i]) > plane.dot(plane_origin) + tol_t - tickness:
                idx[l].append(i)

    n_states = 120

    make_fig(0,n_states)
    make_fig(int(n_states/5),n_states)
    make_fig(int(2*n_states/5),n_states)
    make_fig(int(3*n_states/5),n_states)
    make_fig(int(4*n_states/5),n_states)
    make_fig(n_states,n_states)

    if display:
        plt.show()
    else:
        make_anim(n_states)


if __name__=='__main__':
    p2ptrans(*readOptions())
