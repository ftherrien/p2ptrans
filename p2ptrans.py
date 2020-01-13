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

colorlist=['#929591', 'r', 'k','b','#06470c','#ceb301', '#9e0168', '#26f7fd', '#f97306', '#c20078']
reccolor=['blue','green','red']

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
    parser.add_argument("-r","--noprim",dest="prim",action="store_false", default=True, help="Finds the primitive cell at the beginning") #TMP
    parser.add_argument("-a","--anim",dest="anim",action="store_true", default=False, help="Produce the animation") #TMP
    parser.add_argument("-v","--vol",dest="vol",action="store_true", default=False, help="Make the two (stochiometric) cells equal in volume")
    

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
    anim = options.anim
    vol = options.vol
    
    return fileA, fileB, ncell, filename, display, outdir, use, switch, prim, anim, vol

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
                            min_uvw[:,l] = np.array([i,j,k]).T

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

            # for j in range(len(pos.T)):
            #     if abs(pos[2,idx[j]]) < 1:
            #         print("---", pos[:,idx[j]].T)
            #     elif abs(pos[2,idx[j]]) > 13.1:
            #         print("+++", pos[:,idx[j]].T)
            #     else:
            #         print(pos[:,idx[j]].T)

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
                                    print("Found cell!")
                                    return newcell, origin
                    else:
                        continue
                    break
                else:
                    continue
                break
            print("WARNING: Could not find periodic cell using displacement %d. Increase sample size or use results with care."%i)
        print("WARNING: Could not find cell using shortest distances, trying random order") 
    raise RuntimeError("Could not find periodic cell for any displacement. Increase sample size.")                

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

def p2ptrans(fileA, fileB, ncell, filename, display, outdir, use, switch, prim, anim, vol):

    on_top = None

    os.makedirs(outdir, exist_ok=True)

    if not display:
        matplotlib.use('Agg')

    import matplotlib.pyplot as plt

    plt.rcParams["figure.figsize"] = [5, 5]
    
    import matplotlib.gridspec as gridspec
    from matplotlib.patches import Rectangle

    def add_panel(fig,g,plane, state, p, anchor):
        ax = fig.add_subplot(g, projection='3d', proj_type = 'ortho')
        ax.set_anchor(anchor)
        # ax.view_init(*angles)
        ax.view_init(azim=-90, elev=90) # x-y plane view
        maxXAxis = np.abs([c for a in Tpos for b in a for c in b.flatten()]).max() + 1
        ax.set_xlim([-maxXAxis, maxXAxis])
        ax.set_ylim([-maxXAxis, maxXAxis])
        ax.set_zlim([-maxXAxis, maxXAxis])
        
        print("ROTATION", set_view(plane))
        toplot = set_view(-plane).dot(Tpos[state][p])
        color_to_plot = np.array(color_array[state][p])
        # idxx = np.argsort(toplot.T.dot(transStruc[state].cell.dot(viewDirs[p]).reshape(3,1)), axis=0).T[0]
        idxx = np.argsort(toplot[2,:])
        print("IDXX", idxx)
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
        

    def all_panels(fig, gs, state, label):
        for p,pl in enumerate(viewDirs):
            if isinstance(pl[0],(list, np.ndarray)): 
                plane = normal(transStruc[state].cell).dot(pl[0] + state/n_steps*(pl[1] - pl[0]))
            else:
                plane = normal(transStruc[state].cell).dot(pl)
            plane = plane/la.norm(plane)
            #angles = dir2angles(plane)
            if p ==0:
                add_panel(fig,gs[0,0], plane, state, p, 'E')
            elif p == 1:
                add_panel(fig,gs[0,1], plane, state, p, 'W')
            else:
                add_panel(fig,gs[1,0], plane, state, p, 'E')
            
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
            print("ORIGIN", origin)
            ax.quiver(*(-origin), *plane, color=reccolor[p])
            # ax.quiver(*np.zeros(3), *plane, pivot='tip', color=reccolor[p])
            
        fig.suptitle("Space group: " + spgList[state], fontsize=16)
        for x in ax.get_children():
            if isinstance(x, matplotlib.legend.Legend):
                break
        else:
            fig.legend()

    def make_fig(state, save):
        fig = plt.figure(figsize=[7.2,7.2])
        # fig = plt.figure()
        gs = gridspec.GridSpec(2, 2)
        gs.update(wspace=0.03, hspace=0.03)
        all_panels(fig,gs, state, True)
        if save:
            fig.savefig(outdir+"/Trans_%d.svg"%state)
    
    def make_anim(n_states):
        fig = plt.figure(figsize=[12.8,7.2])
        gs = gridspec.GridSpec(2, 2)
        gs.update(wspace=0.03, hspace=0.03)
        def animate_trans(state):
            all_panels(fig,gs, state, state==0)

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

    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    A = read.poscar(fileA)
    B = read.poscar(fileB)

    if prim:
        print("Making the cells primitive.")
        A = primitive(A, tol)
        B = primitive(B, tol)
    
    print("POSCARS")
    print(A)
    print(B)

    mul = lcm(len(A),len(B))
    mulA = mul//len(A)
    mulB = mul//len(B)

    Acell = A.cell*float(A.scale)
    Bcell = B.cell*float(B.scale)
            
    if (abs(mulA*la.det(Acell)) < abs(mulB*la.det(Bcell))) != switch: # (is switched?) != switch
        print("Transition from %s to %s"%(fileB, fileA))
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
        print("Transition from %s to %s"%(fileA, fileB))

    if vol:
        normalize = (abs(mulA*la.det(Acell)) / abs(mulB*la.det(Bcell)))**(1./3.)
        Bcell = normalize*Bcell
        B.cell = Bcell
        for b in B:
            b.pos = normalize*b.pos

    initSpg = get_spacegroup(to_spglib(A), symprec=0.3, angle_tolerance=3.0)
    finalSpg = get_spacegroup(to_spglib(B), symprec=0.3, angle_tolerance=3.0)
            
    print("Initial SpaceGroup:", initSpg)
    print("Final SpaceGroup:", finalSpg)

    print("Number of Bcells in sphere:", mulB*ncell)
    print("Number of Acells in sphere:", mulA*ncell)
    print("Number of atoms in sphere:", mulA*ncell*len(A))
    
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
            Apos_map, Bpos, Bposst, n_map, natA, class_list, tmat, dmin, vec = tr.fastoptimization(Apos, Bpos, Acell, la.inv(Acell), mulA * la.det(Acell)/(mulB * la.det(Bcell)), atoms, filename)
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

        class_list = class_list[:n_map]-1

        print("Number of classes:", len(np.unique(class_list)))

        Bpos = Bpos[:,:n_map]
        Bposst = Bposst[:,:n_map]
        Apos_map = Apos_map[:,:n_map]

        for class_type in np.unique(class_list):
            print("Class:", class_type)
            print("n:", np.sum(class_list == class_type))
            id = np.where([class_list==class_type])[1][0]
            print("value:", Apos_map[:,id] - Bposst[:,id], la.norm(Apos_map[:,id] - Bposst[:,id])) 
            
        print("Number of mapped atoms:", n_map)
        print("Total distance between structures:", dmin)
        dmin = sum(np.sqrt(sum((Apos_map - Bpos)**2,0)))
        print("Total distance between structures in python:", dmin)

        natB = n_map // np.sum(atoms)
        nat_map = n_map // np.sum(atoms)
        nat = np.shape(Apos)[1] // np.sum(atoms)

        try:
            print("Looking for periodic cell...")        
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
    ax.view_init(0,90) # TMP
    # ax.scatter(Apos.T[:,0],Apos.T[:,1])
    num_tot = 0

    for i,num in enumerate(atoms):
        ax.scatter(Apos.T[num_tot*nat:num_tot*nat+natA*num+1,0],Apos.T[num_tot*nat:num_tot*nat+natA*num+1,1],Apos.T[num_tot*nat:num_tot*nat+natA*num+1,2], c=colorlist[2*i])
        ax.scatter(Apos.T[num_tot*nat+natA*num:(num_tot + num)*nat+1,0],Apos.T[num_tot*nat+natA*num:(num_tot + num)*nat+1,1],Apos.T[num_tot*nat+natA*num:(num_tot + num)*nat+1,2], c=colorlist[2*i], alpha=0.1)
        ax.scatter(Bpos.T[natB*num_tot:natB*(num_tot+num),0],Bpos.T[natB*num_tot:natB*(num_tot+num),1], Bpos.T[natB*num_tot:natB*(num_tot+num),2], c=colorlist[2*i+1])
        num_tot = num_tot + num
    
    centerofmassA = np.mean(Apos,axis=1)
    centerofmassB = np.mean(Bpos,axis=1)

    ax.scatter(centerofmassA[0], centerofmassA[1], centerofmassA[2], s=60, c='red')
    ax.scatter(centerofmassB[0], centerofmassB[1], centerofmassB[2], s=60, c='green')

    maxXAxis = np.max([Apos.max(), Bpos.max()]) + 1
    ax.set_xlim([-maxXAxis, maxXAxis])
    ax.set_ylim([-maxXAxis, maxXAxis])
    ax.set_zlim([-maxXAxis, maxXAxis])
        
    
    # Displacements without stretching (for plotting)
    disps_total = Apos_map - Bpos
    
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    ax.quiver(Bpos.T[:,0], Bpos.T[:,1], Bpos.T[:,2], disps_total.T[:,0], disps_total.T[:,1], disps_total.T[:,2])
    maxXAxis = np.max([Apos.max(), Bpos.max()]) + 1
    ax.set_xlim([-maxXAxis, maxXAxis])
    ax.set_ylim([-maxXAxis, maxXAxis])
    ax.set_zlim([-maxXAxis, maxXAxis])
    
    fig.savefig(outdir+'/DispLattice.svg')
    
    # Displacement with stretching
    disps = Apos_map - Bposst

    vec_classes = np.array([np.mean(disps[:,class_list==d_type], axis=1) for d_type in np.unique(class_list)])

    if pca:
        print("PCA found %d classes"%PCA(disps))
    
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.view_init(0,90) # TMP
    maxXAxis = np.max([Apos.max(), Bposst.max()]) + 1
    ax.set_xlim([-maxXAxis, maxXAxis])
    ax.set_ylim([-maxXAxis, maxXAxis])
    ax.set_zlim([-maxXAxis, maxXAxis])
    

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

    if anim:
        anim = animation.FuncAnimation(fig, animate, init_func=init_disps,
                                       frames=490, interval=30)
        anim.save(outdir+'/Crystal+Disps.gif', fps=30, codec='gif')
    else:
        init_disps()
    

    # Plotting just the displacements
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    maxXAxis = disps.max()
    ax.set_xlim([-maxXAxis, maxXAxis])
    ax.set_ylim([-maxXAxis, maxXAxis])
    ax.set_zlim([-maxXAxis, maxXAxis])
    
    for i in range(len(vec_classes)):
        disps_class = disps[:,class_list==i]
        ndisps = np.shape(disps_class)[1]
        ax.quiver(np.zeros((1,ndisps)), np.zeros((1,ndisps)), np.zeros((1,ndisps)), disps_class.T[:,0], disps_class.T[:,1], disps_class.T[:,2], color=colorlist[i%10])
    fig.savefig(outdir+'/DispOverlayed.svg')
    
    # Centers the position on the first atom    
    print("Volume stretching factor:", la.det(tmat))
    print("Cell volume ratio (should be exactly the same):", mulA * la.det(Acell)/(mulB * la.det(Bcell)))
        
    print("Showing")

    if display:
        plt.show()

    if not found:
        raise RuntimeError("Could not find good displacement cell. Increase system size")

    cell = foundcell
    
    pos_in_struc = Bposst - origin.dot(np.ones((1,np.shape(Bposst)[1])))

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

    print("LEN1", len(dispStruc))

    # Makes sure it is the primitive cell 
    dispStruc = primitive(dispStruc, tolerance = tol)

    print("LEN2", len(dispStruc))

    tmpStruc = Structure(dispStruc.cell)
    to_add = [np.mod(la.inv(dispStruc.cell).dot(a.pos)+tol,1)-tol for a in stinitStruc]
    for idx, pos in zip(*uniqueclose(np.array(to_add).T, tol)):
        tmpStruc.add_atom(*dispStruc.cell.dot(pos),stinitStruc[idx[0]].type)
    
    stinitStruc = tmpStruc

    write.poscar(stinitStruc, vasp5=True, file="POSCAR_init")
    
    assert len(stinitStruc) == len(dispStruc)

    cell = dispStruc.cell

    print("VOLUME", la.det(dispStruc.cell))
    
    finalStruc = Structure(dispStruc.cell)
    for i,a in enumerate(dispStruc):
        print("1", a)
        print("2", a.pos+vec_classes[int(a.type)])
        finalStruc.add_atom(*(a.pos+vec_classes[int(a.type)]),stinitStruc[i].type)

    print("Is it a supercell?")
    print(la.inv(Acell).dot(dispStruc.cell))
    print(la.inv(tmat.dot(Bcell)).dot(dispStruc.cell))

    print("Number of A cell in dispCell:", la.det(dispStruc.cell)/(mulA*la.det(Acell)))
    print("Number of B cell in dispCell:", la.det(dispStruc.cell)/(mulB*la.det(tmat.dot(Bcell))))

    print("dispCell in A frame:")
    print("A coord")
    print(la.inv(Acell).dot(dispStruc.cell))
    print(find_uvw(la.inv(Acell).dot(dispStruc.cell)))
    print("Cartesian Coord.")
    print(dispStruc.cell)
    print(find_uvw(dispStruc.cell))
    print("dispCell in B frame:")
    print("B coord.")
    print(la.inv(Bcell).dot(la.inv(tmat).dot(dispStruc.cell)))
    print(find_uvw(la.inv(Bcell).dot(la.inv(tmat).dot(dispStruc.cell))))
    print("Cartesian Coord.")
    print(la.inv(tmat).dot(dispStruc.cell))
    print(find_uvw(la.inv(tmat).dot(dispStruc.cell)))
    
    # Total displacement per unit volume a as metric
    Total_disp = 0
    for disp in dispStruc:
        Total_disp += la.norm(vec_classes[int(disp.type)])
    
    Total_disp = Total_disp / la.det(dispStruc.cell)
    
    print("Total displacement stretched cell:", Total_disp)

    # Growth directions

    stMat = rotate(tmat, np.eye(3)).dot(tmat)

    eigval, eigvec = la.eig(stMat)

    # closest uvw for -10 to 10

    print("Exact directions of stretching in B frame")
    print("B  coord.")
    print(la.inv(Bcell).dot(eigvec))
    print(find_uvw(la.inv(Bcell).dot(eigvec)))
    print("Cartesian  coord.")
    print(eigvec)         
    print(find_uvw(eigvec))
    print("Amount of stretching")
    print(eigval)    
    
    neweigval, P = la.eig(tmat.T.dot(tmat))
    neweigval = np.sqrt(neweigval)
    
    P = normal(P)

    U = P.dot(np.diag(neweigval)).dot(P.T)
    
    print("Deformation Gradient Method (FtF)")

    print("Exact directions of stretching in final frame (%s) coord."%finalSpg)
    print("B  coord.")
    print(la.inv(Bcell).dot(P))
    print(find_uvw(la.inv(Bcell).dot(P)))
    print("Cartesian  coord.")
    print(P)
    print(find_uvw(P))
    print("Amount of stretching")
    print(neweigval)

    neweigval, Q = la.eig(la.inv(tmat).T.dot(la.inv(tmat)))
    neweigval = np.sqrt(neweigval)
    idx = np.argsort(1/neweigval)
    #Q = Q[:,idx]
    
    Q = normal(Q)
    
    print("Exact directions of stretching in intial frame (%s) coord."%initSpg)
    print("A coord.")
    print(la.inv(Acell).dot(Q))
    print(find_uvw(la.inv(Acell).dot(Q)))
    print("cartesian coord.")
    print(Q)
    print(find_uvw(Q))
    print("Amount of stretching 1/lmabda")
    print(1/neweigval)

    #Plane (miller indices) and direction in final frame or B coord
    useBcoord = False
    recB = la.inv(Bcell).T
    planeDir = [1,1,0]
    orDir = [1,-1,-1]
    if useBcoord:
        planeDir = recB.dot(planeDir)
        orDir = Bcell.dot(orDir)

    print("R", tmat.dot(P.T.dot(np.diag(neweigval)).dot(P)))
    
    print("Orientation Relationship for thin film")
    resPlaneDir = Q.dot(P.T).dot(planeDir)
    print("A/B coord.")
    print("Planes:", la.inv(Bcell).dot(planeDir), la.inv(Acell).dot(resPlaneDir))
    print("Planes:", la.inv(Bcell).dot(planeDir), find_uvw(la.inv(Acell).dot(resPlaneDir.reshape((3,1)))).T[0])
    resOrDir = Q.dot(P.T).dot(orDir)
    print("Directions:", la.inv(Bcell).dot(orDir), la.inv(Acell).dot(resOrDir))
    print("Directions:", la.inv(Bcell).dot(orDir), find_uvw(la.inv(Acell).dot(resOrDir.reshape((3,1)))).T[0])
    print("Cartesian coord.")
    print("Planes:", planeDir, resPlaneDir)
    print("Planes:", planeDir, find_uvw(resPlaneDir.reshape((3,1))).T[0])
    resOrDir = Q.dot(P.T).dot(orDir)
    print("Directions:", orDir, resOrDir)
    print("Directions:", orDir, find_uvw(resOrDir.reshape((3,1))).T[0])

    twin = Q.dot(np.diag((1,-1,1))).dot(Q.T)
    
    print("Plus twinning")
    print("Plane:", la.inv(Acell).dot(twin.dot(resPlaneDir)))
    print("Directions:", la.inv(Acell).dot(twin.dot(resOrDir)))

    print("Habit Plane (Closest to)")
    idv = np.argmin(abs(eigval - 1))
    idh = np.arange(3)[np.arange(3)!=idv]
    ratio = np.sqrt(abs((eigval[idh[1]]**2 - eigval[idv]**2)/(eigval[idv]**2 - eigval[idh[0]]**2)))
    #ratio = np.sqrt(abs((eigval[idh[1]]**2 - 1)/(1 - eigval[idh[0]]**2)))

    # Only for cubic systems!
    
    planeHab = np.zeros((3,2))
    planeHab[:,0] = eigvec[:,idh[0]] + ratio*eigvec[:,idh[1]]
    planeHab[:,1] = eigvec[:,idh[0]] - ratio*eigvec[:,idh[1]]
    
    if useBcoord:
        print("Exact Habit plane:", la.inv(Bcell).dot(planeHab))
        print("Closest uvw:", find_uvw(la.inv(Bcell).dot(planeHab)))
    else:
        print("Exact Habit plane:", planeHab)
        print("Closest uvw:", find_uvw(planeHab))

    planeHab=planeHab[:,0]

    vhab = np.zeros((3,3))
    rhab = np.zeros((3,3))
    
    vhab[:,0] = ratio*eigvec[:,idh[0]] - eigvec[:,idh[1]]
    vhab[:,1] = eigvec[:,idv]
    vhab[:,2] = planeHab/la.norm(planeHab) 
    
    print("test 1", la.norm(vhab,axis=0), la.norm(tmat.dot(vhab), axis=0)) 
    
    print("Orientation relationship with Habit Plane")
    rhab[:,:2] = tmat.dot(vhab[:,:2])
    rhab[:,2] = np.cross(rhab[:,0], rhab[:,1])
    rhab[:,2] = rhab[:,2] / la.norm(rhab[:,2])

    r2hab = np.zeros((3,3))
    
    r2hab[:,:2] = U.dot(vhab[:,:2])
    r2hab[:,2] = np.cross(r2hab[:,0], r2hab[:,1])
    r2hab[:,2] = r2hab[:,2] / la.norm(r2hab[:,2])


    print("P", P)

    print("Q", Q)

    print("QQ'", Q.dot(Q.T))

    print("PP'", P.dot(P.T))

    print("QP'", Q.dot(P.T)) 

    print('RR', tmat.dot(la.inv(U)))
    
    Rtest = tmat.dot(la.inv(U)).dot(r2hab.dot(la.inv(vhab)))
    R = rhab.dot(la.inv(vhab))
    
    def find_R_RU(mat):
        if np.all(mat==np.eye(3)):
            return mat

        eigval, eigvec = la.eig(mat)
        idv = np.argmin(abs(eigval - 1))
        idh = np.arange(3)[np.arange(3)!=idv]
        ratio = np.sqrt(abs((eigval[idh[1]]**2 - eigval[idv]**2)/(eigval[idv]**2 - eigval[idh[0]]**2)))

        # Only for cubic systems!
        
        planeHab = np.zeros((3,2))
        planeHab = eigvec[:,idh[0]] + ratio*eigvec[:,idh[1]]
        
        vhab = np.zeros((3,3))
        rhab = np.zeros((3,3))
        
        vhab[:,0] = ratio*eigvec[:,idh[0]] - eigvec[:,idh[1]]
        vhab[:,1] = -eigvec[:,idv]
        vhab[:,2] = planeHab/la.norm(planeHab) 
        
        r2hab[:,:2] = mat.dot(vhab[:,:2])
        r2hab[:,2] = np.cross(r2hab[:,0], r2hab[:,1])
        r2hab[:,2] = r2hab[:,2] / la.norm(r2hab[:,2])

        return r2hab.dot(la.inv(vhab)).T

    R = Q.dot(P.T).dot(r2hab.dot(la.inv(vhab)))

    print('R', R)
    print('Rtest', Rtest)
    
    print("test habit", R.dot(vhab), la.inv(R).dot(vhab), tmat.dot(vhab))
    resPlaneDir = R.dot(planeDir)
    print("A/B coord.")
    print("Planes:", la.inv(Bcell).dot(planeDir), la.inv(Acell).dot(resPlaneDir))
    print("Planes:", la.inv(Bcell).dot(planeDir), find_uvw(la.inv(Acell).dot(resPlaneDir.reshape((3,1)))).T[0])
    resOrDir = R.dot(orDir)
    print("Directions:", la.inv(Bcell).dot(orDir), la.inv(Acell).dot(resOrDir))
    print("Directions:", la.inv(Bcell).dot(orDir), find_uvw(la.inv(Acell).dot(resOrDir.reshape((3,1)))).T[0])
    print("Cartesian coord.")
    print("Planes:", planeDir, resPlaneDir)
    print("Planes:", planeDir, find_uvw(resPlaneDir.reshape((3,1))).T[0])
    resOrDir = R.dot(orDir)
    print("Directions:", orDir, resOrDir)
    print("Directions:", orDir, find_uvw(resOrDir.reshape((3,1))).T[0])
    
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
        ax.quiver(np.zeros(3), np.zeros(3), np.zeros(3), Acell[0,:], Acell[1,:], Acell[2,:], color = "blue", alpha = 0.3)
        maxXAxis = abs(cell).max() + 1
        ax.set_xlim([-maxXAxis, maxXAxis])
        ax.set_ylim([-maxXAxis, maxXAxis])
        ax.set_zlim([-maxXAxis, maxXAxis])
        
        fig.savefig(outdir+'/Displacement_structure.svg')
        return fig,
    
    if anim:
        anim = animation.FuncAnimation(fig, animate, init_func=init_struc,
                                       frames=490, interval=30)
        anim.save(outdir+'/DispStruc.gif', fps=30, codec='gif')
    else:
        init_struc()

    # Produce transition

    print("Producing transition!")

    n_steps = 60

    PoscarDirName = "/TransPOSCARS_K-S"
    
    os.makedirs(outdir+PoscarDirName, exist_ok=True)

    itmat = rotate(la.inv(tmat).dot(finalStruc.cell), finalStruc.cell).dot(la.inv(tmat))

    print("ITMAT", itmat)
    print("U-1", Q.dot(P.T).dot(la.inv(U)).dot(P.dot(Q.T)))
    
    spgList = []
    transStruc = []
    # viewirs = [[la.inv(finalStruc.cell).dot([0,0,1]),
    #             la.inv(la.inv(tmat).dot(finalStruc.cell)).dot([1,0,0])],
    #             [1,0,0], [0,1,0]]

    # viewDirs = [[0, np.sqrt(3), 3], [1,0,0], [0, np.sqrt(3), -1]]

    # viewDirs = np.eye(3).tolist()

    viewDirs = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    
    print("Viewing directions in units of dispCell:", viewDirs)
    # print("Choosing the right plane", la.inv(tmat).dot(finalStruc.cell).dot(la.inv(finalStruc.cell).dot([0,0,1])))
    color_array = []
    size = 3
    Tpos = [] 
    for i in range(n_steps+1):
        # Pitsch
        curMat = (itmat-np.eye(3))*i/n_steps + np.eye(3)
        # K-S
        curMat = find_R_RU(curMat).dot(curMat) # K-S only remove this line for Pitsch
        curStruc = Structure(curMat.dot(finalStruc.cell))
        for j,a in enumerate(dispStruc):
            curDisp = vec_classes[int(a.type)]*i/n_steps
            curPos = curMat.dot((finalStruc[j].pos - curDisp).reshape((3,1)))
            curStruc.add_atom(*(curPos.T.tolist()[0]),finalStruc[j].type)
            # curStruc = supercell(curStruc, curStruc.cell)
        write.poscar(curStruc, vasp5=True, file=outdir+PoscarDirName+"/POSCAR_%03d"%i)
        spgList.append(get_spacegroup(to_spglib(curStruc), symprec=0.3, angle_tolerance=3.0))
        transStruc.append(curStruc)
        
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

    print("Spacegroups along the transition:", spgList)
    
    # Displays displacement with the disp cell overlayed
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.quiver(pos_in_struc.T[:,0], pos_in_struc.T[:,1], pos_in_struc.T[:,2], disps.T[:,0], disps.T[:,1], disps.T[:,2], color = "C0")
    ax.scatter(pos_in_struc.T[:,0], pos_in_struc.T[:,1], pos_in_struc.T[:,2], s=10, color = "C0")
    ax.quiver(np.zeros(3), np.zeros(3), np.zeros(3), foundcell[0,:], foundcell[1,:], foundcell[2,:], color = "red")
    ax.quiver(np.zeros(3), np.zeros(3), np.zeros(3), cell[0,:], cell[1,:], cell[2,:], color = "green")
    maxXAxis = pos_in_struc.max() + 1
    ax.set_xlim([-maxXAxis, maxXAxis])
    ax.set_ylim([-maxXAxis, maxXAxis])
    ax.set_zlim([-maxXAxis, maxXAxis])
    fig.savefig(outdir+'/DispLattice_stretched_cell_primitive.svg')

    plt.show() #TMP

        
    # Showing some of the steps
    make_fig(0, True)
    make_fig(int(n_steps/4), True)
    make_fig(int(2*n_steps/4), True)
    make_fig(int(3*n_steps/4), True)
    make_fig(n_steps, True)

    if display:
        print("Showing")
        plt.show()

    if anim:
        make_anim(n_steps)


if __name__=='__main__':
    p2ptrans(*readOptions())
