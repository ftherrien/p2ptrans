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
import re
from glob import iglob
import warnings
from format_spglib import from_spglib, to_spglib
from spglib import get_spacegroup

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
    parser.add_argument("-O","--orig",dest="C",type=str, default='./POSCAR_C', help="Original Structure")
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
    fileC = options.C
    ncell = options.ncell
    filename = options.filename
    display = options.display
    outdir = options.outdir
    use = options.use
    switch = options.switch
    prim = options.prim
    
    return fileA, fileB, fileC, ncell, filename, display, outdir, use, switch, prim
    
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

def p2ptrans(fileA, fileB, fileC, ncell, filename, display, outdir, use, switch, prim):

    on_top = None

    os.makedirs(outdir, exist_ok=True)

    if not display:
        matplotlib.use('Agg')

    import matplotlib.pyplot as plt

    import matplotlib.gridspec as gridspec
    from matplotlib.patches import Rectangle

    def add_panel(fig,g,angles, state, p, label, anchor):
        color_list=[]
        ax = fig.add_subplot(g, projection='3d', proj_type = 'ortho')
        ax.set_anchor(anchor)
        ax.view_init(*angles)
        maxXAxis = np.abs([c for a in Tpos for b in a for c in b.flatten()]).max() + 1
        ax.set_xlim([-maxXAxis, maxXAxis])
        ax.set_ylim([-maxXAxis, maxXAxis])
        ax.set_zlim([-maxXAxis, maxXAxis])
        ax.set_aspect('equal')
        toplot = Tpos[state][p]
        color_to_plot = np.array(color_array[state][p])
        idxx = np.argsort(toplot.T.dot(transStruc[state].cell.dot(planes[p]).reshape(3,1)), axis=0).T[0]
        toplot = toplot[:,idxx]
        color_to_plot = color_to_plot[idxx]
        ax.set_axis_off()
        ax.dist = 2
        axlims = [a for b in ax.get_position().get_points() for a in b]
        rec = Rectangle((axlims[0],axlims[1]),(axlims[2]-axlims[0]),(axlims[3]-axlims[1]), transform = fig.transFigure, fill=False,lw=1, color="k")
        fig.patches.append(rec)
        for i,point in enumerate(toplot.T):
            if any(color_list==color_to_plot[i]) or not label or p:
                ax.scatter(*point, c=colorlist[color_to_plot[i]], s=40, depthshade=False)
            else:
                color_list.append(color_to_plot[i])
                ax.scatter(*point, c=colorlist[color_to_plot[i]], s=40, depthshade=False, label=atom_types[color_to_plot[i]])

    def all_panels(fig, gs, state, label):
        for p,pl in enumerate(planes):
            plane = transStruc[state].cell.dot(pl)
            plane = plane/la.norm(plane)
            angles = dir2angles(plane)
            if p ==0:
                add_panel(fig,gs[0,0], angles, state, p, label, 'E')
            elif p == 1:
                add_panel(fig,gs[0,1], angles, state, p, label, 'W')
            else:
                add_panel(fig,gs[1,0], angles, state, p, label, 'E')
            
        ax = fig.add_subplot(gs[1,1], projection='3d', proj_type = 'ortho')
        ax.set_anchor('W')
        ax.set_axis_off()
        ax.view_init(*(np.array(dir2angles(transStruc[state].cell[:,1]))+np.array([0,0])))
        ax.dist = 4
        maxXAxis = abs(transStruc[state].cell).max() + 1
        ax.set_xlim([-maxXAxis, maxXAxis])
        ax.set_ylim([-maxXAxis, maxXAxis])
        ax.set_zlim([-maxXAxis, maxXAxis])
        ax.set_aspect('equal')
        axlims = [a for b in ax.get_position().get_points() for a in b]
        rec = Rectangle((axlims[0],axlims[1]),(axlims[2]-axlims[0]),(axlims[3]-axlims[1]), transform = fig.transFigure, fill=False,lw=1, color="k")
        fig.patches.append(rec)

        origin = np.sum(transStruc[state].cell, axis=1)/2

        for i in range(3):
            base = np.array([np.zeros(3), transStruc[state].cell[:,(i+1)%3],
                             transStruc[state].cell[:,(i+2)%3], 
                             transStruc[state].cell[:,(i+1)%3] + transStruc[state].cell[:,(i+2)%3]])
            vec = transStruc[state].cell[:,i:i+1].dot(np.ones((1,4)))
            ax.quiver(base[:,0]-origin[0], base[:,1]-origin[1], base[:,2]-origin[2], vec[0,:], vec[1,:], vec[2,:], arrow_length_ratio=0, color="k", alpha=0.5)
        for a in transStruc[state]:
            ax.scatter(a.pos[0]-origin[0], a.pos[1]-origin[1], a.pos[2]-origin[2], alpha = 1, s=400, color=colorlist[(np.where(atom_types == a.type)[0][0])%10])

        fig.suptitle("Space group: " + spgList[state], fontsize=16)
        fig.legend()

    def make_fig(state):
        fig = plt.figure(figsize=[12.8,7.2])
        gs = gridspec.GridSpec(2, 2)
        gs.update(wspace=0.01, hspace=0.01)
        all_panels(fig,gs, state, True)
    
    def make_anim(n_states):
        fig = plt.figure(figsize=[12.8,7.2])
        gs = gridspec.GridSpec(2, 2)
        gs.update(wspace=0.01, hspace=0.01)
        def animate_trans(state):
            all_panels(fig,gs, state, state==1)

        # animation.verbose.set_level('debug')
    
        plt.rcParams['animation.ffmpeg_path'] = '/home/felixt/projs/bin/ffmpeg'
        # Writer = animation.writers['ffmpeg']
        writer = animation.FFMpegWriter(fps=30,codec='prores', extra_args=['-loglevel', 'verbose','-f','mov'])
    
        anim = animation.FuncAnimation(fig, animate_trans,
                                   frames=n_states+1, interval=1)
        anim.save(outdir + '/Trans.mov', writer=writer)

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

    # Finding the right mapping
    transfiles = []
    idx = []
    for i,ff in enumerate(iglob(fileA)):
        idx.append(int(re.findall("([0-9]+)[^/]*$",ff)[-1]))
        transfiles.append(ff)

    idx = np.argsort(idx)
    if switch:
        idx = idx[::-1]
    transfiles = [transfiles[i] for i in idx]

    print("Reading A from %s"%(transfiles[0]))
    A = read.poscar(transfiles[0])
    print("Reading B from %s"%(transfiles[-1]))
    B = read.poscar(transfiles[-1])

    mapAdjust = np.zeros((len(A),3))
    
    for i,ff in enumerate(transfiles):
        S = read.poscar(ff)
        if i:
            iScell = la.inv(S.cell)
            iSprevcell = la.inv(Sprev.cell)
            for j,s in enumerate(S): 
                spos = iScell.dot(s.pos)
                sposprev = iSprevcell.dot(Sprev[j].pos)
                trans = -(np.round(spos)*2 - 1)
                cond = np.argmin([abs(spos - sposprev), abs(spos + trans - sposprev)],0)
                if any(cond):
                    print("FOUND CROSSING")
                    mapAdjust[j] += -cond*trans
        Sprev = S

    print("ADJUSTMENTS TO CELL:")
    print(mapAdjust)

    for i,a in enumerate(A):
        a.pos = a.pos + A.cell.dot(mapAdjust[i]) 

    C = read.poscar(fileC)

    print("POSCARS")
    print("A", A)
    print("B", B)
    print("C", C)
    
    # B should be a supercell of C
    
    Acell = A.cell*float(A.scale)
    Bcell = B.cell*float(B.scale)
    
    if prim:
        print("Making the original cell primitive.")
        C = primitive(C, tol)
    
    Ccell = C.cell*float(C.scale)

    ncell_match = 30

    mulC = abs(la.det(Bcell)/la.det(Ccell))
    print("The ratio between the final structure and the original one is:", mulC) 
    assert abs(mulC - round(mulC)) < 1e-5
    mulC = int(round(mulC))
    mulB = 1
    print("There will be %d cells of original structure"%(mulC*ncell))
    print("Initial SpaceGroup:", get_spacegroup(to_spglib(B), symprec=0.3, angle_tolerance=3.0))
    print("Final SpaceGroup:", get_spacegroup(to_spglib(A), symprec=0.3, angle_tolerance=3.0))
    
    tmat = np.eye(3)
    found = False
    rep = 0
    while (not found and rep < 1):
        rep += 1
        found = True

        # Adds atoms to C and B (for cell with different types of atoms)
        Cpos = []
        atom_types = np.array([], np.str)
        atomsC = np.array([], np.int)
        for a in C:
            if any(atom_types == a.type):
                idx = np.where(atom_types == a.type)[0][0]
                Cpos[idx] = np.concatenate((Cpos[idx], t.sphere(Ccell, mulC*ncell_match, a.pos*float(C.scale))), axis = 1) 

                # Order the atoms in terms of distance
                Cpos[idx] = Cpos[idx][:,np.argsort(la.norm(Cpos[idx],axis=0))] 
                atomsC[idx] += 1
            else:
                print(Ccell, mulC*ncell_match, a.pos*float(C.scale))
                Cpos.append(t.sphere(Ccell, mulC*ncell_match, a.pos*float(C.scale)))
                atom_types = np.append(atom_types, a.type)
                atomsC = np.append(atomsC,1)
        
        Cpos = np.concatenate(Cpos, axis=1)
        
        # Temporarly stretching Bcell, for tiling
        Bcell = tmat.dot(Bcell)

        Bpos = [None]*len(atom_types)
        atomsB = np.zeros(len(atom_types), np.int)
        for a in B:
            idx = np.where(atom_types == a.type)[0][0]
            if atomsB[idx] == 0:
                Bpos[idx] = t.sphere(Bcell, mulB*ncell_match, tmat.dot(a.pos)*float(B.scale))
            else:
                Bpos[idx] = np.concatenate((Bpos[idx], t.sphere(Bcell, mulB*ncell_match, tmat.dot(a.pos)*float(B.scale))), axis = 1) 
                # Order the atoms in terms of distance
                Bpos[idx] = Bpos[idx][:,np.argsort(la.norm(Bpos[idx],axis=0))]
            atomsB[idx] += 1
        
        Bpos = np.concatenate(Bpos, axis=1)

        Bpos = la.inv(tmat).dot(Bpos)
        Bcell = la.inv(tmat).dot(Bcell)

        assert all(mulC*atomsC == mulB*atomsB)
        atoms = mulC*atomsC

        oldcmC = np.mean(Cpos,axis=1)
        oldcmB = np.mean(Bpos,axis=1)

        print(oldcmC, oldcmB)

        if not use:
            Cpos = np.asfortranarray(Cpos)
            Bpos = np.asfortranarray(Bpos)
            tr.center(Bpos)
            oldBpos = np.asanyarray(deepcopy(Bpos))
            t_time = time.time()
            Cpos_map, Bpos, Bposst, n_map, natC, class_list, tmat, dmin, vec = tr.fastoptimization(Cpos, Bpos, Ccell, la.inv(Ccell), la.det(Ccell)/(la.det(Bcell)), atoms, "check.in")
            t_time = time.time() - t_time
            Bpos = np.asanyarray(Bpos)
            Cpos = np.asanyarray(Cpos)
            # --------

            pickle.dump((Cpos_map, Bpos, Bposst, n_map, natC, class_list, tmat, dmin, vec), open(outdir+"/check.dat","wb"))

        else:

            Cpos_map, Bpos, Bposst, n_map, natC, class_list, tmat, dmin, vec = pickle.load(open(outdir+"/check.dat","rb"))

        class_list = class_list[:n_map]-1

        print("Number of classes:", len(np.unique(class_list)))

        Bpos = Bpos[:,:n_map]
        Bposst = Bposst[:,:n_map]
        Cpos_map = Cpos_map[:,:n_map]

        natB = n_map // np.sum(atoms)
        nat_map = n_map // np.sum(atoms)
        nat = np.shape(Cpos)[1] // np.sum(atoms)

        if display:
        
            # Plotting the Cpos and Bpos overlayed
            fig = plt.figure(22)
            ax = fig.add_subplot(111, projection='3d')
            ax.view_init(0,90) # TMP
            # ax.scatter(Cpos.T[:,0],Cpos.T[:,1])
            num_tot = 0
            
            for i,num in enumerate(atoms):
                ax.scatter(Cpos.T[num_tot*nat:num_tot*nat+natC*num+1,0],Cpos.T[num_tot*nat:num_tot*nat+natC*num+1,1],Cpos.T[num_tot*nat:num_tot*nat+natC*num+1,2], c=colorlist[2*i])
                ax.scatter(Cpos.T[num_tot*nat+natC*num:(num_tot + num)*nat+1,0],Cpos.T[num_tot*nat+natC*num:(num_tot + num)*nat+1,1],Cpos.T[num_tot*nat+natC*num:(num_tot + num)*nat+1,2], c=colorlist[2*i], alpha=0.1)
                ax.scatter(Bpos.T[natB*num_tot:natB*(num_tot+num),0],Bpos.T[natB*num_tot:natB*(num_tot+num),1], Bpos.T[natB*num_tot:natB*(num_tot+num),2], c=colorlist[2*i+1])
                num_tot = num_tot + num
            
            centerofmassC = np.mean(Cpos,axis=1)
            centerofmassB = np.mean(Bpos,axis=1)
            
            ax.scatter(centerofmassC[0], centerofmassC[1], centerofmassC[2], s=60, c='red')
            ax.scatter(centerofmassB[0], centerofmassB[1], centerofmassB[2], s=60, c='green')
            
            maxXCxis = np.max([Cpos.max(), Bpos.max()]) + 1
            ax.set_xlim([-maxXCxis, maxXCxis])
            ax.set_ylim([-maxXCxis, maxXCxis])
            ax.set_zlim([-maxXCxis, maxXCxis])
            ax.set_aspect('equal')    
            
            # Displacements without stretching (for plotting)
            disps_total = Cpos_map - Bpos
            
            # fig = plt.figure()
            # ax = fig.add_subplot(111, projection='3d')
            ax.quiver(Bpos.T[:,0], Bpos.T[:,1], Bpos.T[:,2], disps_total.T[:,0], disps_total.T[:,1], disps_total.T[:,2])
            maxXCxis = np.max([Cpos.max(), Bpos.max()]) + 1
            ax.set_xlim([-maxXCxis, maxXCxis])
            ax.set_ylim([-maxXCxis, maxXCxis])
            ax.set_zlim([-maxXCxis, maxXCxis])
            ax.set_aspect('equal')
            fig.savefig(outdir+'/DispLattice.svg')
            
            plt.show()
        
        # --------
        Bcell = tmat.dot(Bcell)
        newBcell = Structure(Bcell)
        for a in B:
            newBcell.add_atom(*(tmat.dot(a.pos.reshape((3,1)) - oldcmB.reshape((3,1))) + oldcmC.reshape((3,1)) + vec.reshape((3,1))).T.tolist()[0], a.type)

        # # TMP not sure
        # for i,a in enumerate(newBcell):
        #     if not i:
        #         shift = Ccell.dot(np.round(la.inv(Ccell).dot(a.pos - tmat.dot(B[i].pos))))
        #     a.pos = a.pos - shift
            
        B = deepcopy(newBcell)

        tol_ch = 1e-3

        adjust = []
        # Into b with tol, adjust A in case B is actually outside the cell
        for b in B:
            print("Into cell:", la.inv(Bcell).dot(b.pos))
            adjust.append(-Acell.dot(np.floor(la.inv(Bcell).dot(b.pos) + tol_ch)))
            b.pos = Bcell.dot(np.mod(la.inv(Bcell).dot(b.pos) + tol_ch, 1) - tol_ch)
            
        Bnew = supercell(C,Ccell.dot(np.round(la.inv(Ccell).dot(B.cell),2)/float(C.scale)))
        sitepos = []
        for i, a in enumerate(deepcopy(Bnew)):
            apos_in = Bnew.cell.dot(np.mod(la.inv(Bnew.cell).dot(a.pos) + tol_ch, 1) - tol_ch)
            for j, b in enumerate(B):
                if a.type == b.type and np.allclose(b.pos, apos_in*float(C.scale), atol = tol_ch):
                    print("Found correspondance for atom %d"%i)
                    Bnew[j] = a
                    adjust[j] = Acell.dot(np.floor(la.inv(Bnew.cell).dot(a.pos) + tol_ch))
                    if a.site == 0:
                        sitepos.append((a.pos*float(C.scale)).reshape((3,1)))
                    break
            else:
                raise RuntimeError("No correspondance found for atom %d!!"%i)

        print("ADDITIONAL ADJUSTMENTS TO CELL")
        print(adjust)

        tol_ch = 1e-10
            
        B = deepcopy(Bnew)

        write.poscar(B,  vasp5=True, file="POSCAR_B")

        Bcell = B.cell*float(B.scale)

        sitepos = [s - sitepos[0] for s in sitepos]

        dminList = []

        print("SITES TO TRY")
        print(sitepos)
        
        for si, site in enumerate(sitepos):
        
            # Adds atoms to C and B (for cell with different types of atoms)
            Cpos = []
            Apos = []
            atom_types = np.array([], np.str)
            atoms = np.array([], np.int)
            for a in C:
                # Csph = np.array([a.pos*float(B.scale) for a in supercell(B,B.cell.dot(np.diag([2,4,2])))]).T
                Csph = t.sphere(Ccell, ncell, a.pos*float(C.scale))
                Bgrid = la.inv(Bcell).dot(Csph - site.dot(np.ones((1,np.shape(Csph)[1]))))
                Asph = []
                for pos in Bgrid.T:
                    for ib, b in enumerate(B):
                        if np.allclose(b.pos*float(B.scale), Bcell.dot(np.mod(pos+tol_ch,1)-tol_ch), atol = tol_ch):
                            Asph.append(Acell.dot(np.floor(pos+tol_ch)) + adjust[ib]  + A[ib].pos*float(A.scale))
                            break
                    else:
                        print(pos)
                        raise RuntimeError("Could not find all atoms")
            
                Asph = np.array(Asph).T
                if any(atom_types == a.type):
                    idx = np.where(atom_types == a.type)[0][0]
            
                    Apos[idx] = np.concatenate((Apos[idx], Asph), axis = 1)
                    Cpos[idx] = np.concatenate((Cpos[idx], Csph), axis = 1) 
            
                    # Order the atoms in terms of distance
                    idx_sort = np.argsort(la.norm(Cpos[idx],axis=0))
                    Cpos[idx] = Cpos[idx][:, idx_sort]
                    Apos[idx] = Apos[idx][:, idx_sort]
                    atoms[idx] += 1
                else:
                    Cpos.append(Csph)
                    Apos.append(Asph)
                    atom_types = np.append(atom_types, a.type)
                    atoms = np.append(atoms,1)
            
            Bpos = np.concatenate(Cpos, axis=1)
            Apos = np.concatenate(Apos, axis=1)
            
            if not False:

                cma = np.mean(Apos,1).reshape((3,1))
            
                Apos = np.asfortranarray(Apos)
                Bpos = np.asfortranarray(Bpos)
                tr.center(Bpos)
                oldBpos = np.asanyarray(deepcopy(Bpos))
                t_time = time.time()
                Apos_map, Bpos, Bposst, n_map, natA, class_list, tmat, dmin, vec = tr.fastoptimization(Apos, Bpos, Acell, la.inv(Acell), la.det(Acell)/(la.det(Bcell)), atoms, filename)
                t_time = time.time() - t_time
                Bpos = np.asanyarray(Bpos)
                Apos = np.asanyarray(Apos)
            
                print("Mapping time:", t_time)
                            
                pickle.dump((Apos, Apos_map, Bpos, Bposst, n_map, natA, class_list, tmat, dmin, vec), open(outdir+"/check_%d.dat"%si,"wb"))
            
            else:
                print("Using data from "+outdir+"/check_%d.dat"%si)
                Apos = np.asfortranarray(Apos)
                Bpos = np.asfortranarray(Bpos) 
                tr.center(Apos)
                tr.center(Bpos)
                oldBpos = np.asanyarray(deepcopy(Bpos))
                Apos, Apos_map, Bpos, Bposst, n_map, natA , class_list, tmat, dmin, vec = pickle.load(open(outdir+"/check_%d.dat"%si,"rb"))
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

            dminList.append(dmin)

        print("Getting the shortest mapping from "+outdir+"/check_%d.dat"%(np.argmin(dminList)))
        Apos = np.asfortranarray(Apos)
        Bpos = np.asfortranarray(Bpos) 
        tr.center(Apos)
        tr.center(Bpos)
        oldBpos = np.asanyarray(deepcopy(Bpos))
        Apos, Apos_map, Bpos, Bposst, n_map, natA , class_list, tmat, dmin, vec = pickle.load(open(outdir+"/check_%d.dat"%(np.argmin(dminList)),"rb"))
        Bpos = np.asanyarray(Bpos)
        Apos = np.asanyarray(Apos)

        print("----------------------------------------------")
        print(" The shortest distance is:", dmin)
        print("----------------------------------------------")
        
        class_list = class_list[:n_map]-1
            
        print("Number of classes:", len(np.unique(class_list)))
        print("Original number of classes:", len(B))
            
        Bpos = Bpos[:,:n_map]
        Bposst = Bposst[:,:n_map]
        Apos_map = Apos_map[:,:n_map]
        
        natB = n_map // np.sum(atoms)
        nat_map = n_map // np.sum(atoms)
        nat = np.shape(Apos)[1] // np.sum(atoms)

        try:
            print("Looking for periodic cell...")
            foundcell, origin = find_cell(class_list, Bposst)
            if abs(abs(la.det(tmat)) - abs(la.det(Acell)/(la.det(Bcell)))) > tol_vol:
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

    #ax.scatter(centerofmassA[0], centerofmassA[1], centerofmassA[2], s=60, c='red')
    #ax.scatter(centerofmassB[0], centerofmassB[1], centerofmassB[2], s=60, c='green')

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
            
    maxXBxis = np.max([Apos.max(), Bpos.max()]) + 1
    ax.set_xlim([-maxXBxis, maxXBxis])
    ax.set_ylim([-maxXBxis, maxXBxis])
    ax.set_zlim([-maxXBxis, maxXBxis])
    ax.set_aspect('equal')
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
        plane = Acell[:,1]
        plane = plane/la.norm(plane)
        angles = dir2angles(plane)
        ax.view_init(*angles)

        origin = cma

        ax.scatter(*(-origin), color='g')
            
        for i in range(3):
            base = np.array([np.zeros(3), Acell[:,(i+1)%3],
                             Acell[:,(i+2)%3], 
                             Acell[:,(i+1)%3] + Acell[:,(i+2)%3]])
            vec = Acell[:,i:i+1].dot(np.ones((1,4)))
            ax.quiver(base[:,0]-origin[0], base[:,1]-origin[1], base[:,2]-origin[2], vec[0,:], vec[1,:], vec[2,:], arrow_length_ratio=0, color="k", alpha=0.5)
    
        fig.savefig(outdir+'/DispLattice_stretched.svg')
        return fig,

    if not display:
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
    ax.set_aspect('equal')
    for i in range(len(vec_classes)):
        disps_class = disps[:,class_list==i]
        ndisps = np.shape(disps_class)[1]
        ax.quiver(np.zeros((1,ndisps)), np.zeros((1,ndisps)), np.zeros((1,ndisps)), disps_class.T[:,0], disps_class.T[:,1], disps_class.T[:,2], color=colorlist[i%10])
    fig.savefig(outdir+'/DispOverlayed.svg')
    
    # Centers the position on the first atom    
    print("Volume stretching factor:", la.det(tmat))
    print("Cell volume ratio (should be exactly the same):", la.det(Acell)/(la.det(Bcell)))
        
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

    print("Number of A cell in dispCell:", la.det(dispStruc.cell)/(la.det(Acell)))
    print("Number of B cell in dispCell:", la.det(dispStruc.cell)/(la.det(tmat.dot(Bcell))))

    # Total displacement per unit volume a as metric
    Total_disp = 0
    for disp in dispStruc:
        Total_disp += la.norm(vec_classes[int(disp.type)])
    
    Total_disp = Total_disp / la.det(dispStruc.cell)
    
    print("Total displacement stretched cell:", Total_disp)

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
        ax.set_aspect('equal')
        fig.savefig(outdir+'/Displacement_structure.svg')
        return fig,
    
    if not display:
        anim = animation.FuncAnimation(fig, animate, init_func=init_struc,
                                       frames=490, interval=30)
        anim.save(outdir+'/DispStruc.gif', fps=30, codec='gif')
    else:
        init_struc()

    # Produce transition

    print("Producing transition!")

    n_steps = 100

    os.makedirs(outdir+"/TransPOSCARS", exist_ok=True)

    itmat = rotate(la.inv(tmat).dot(finalStruc.cell), finalStruc.cell).dot(la.inv(tmat))
    
    print("tmat", la.det(tmat), tmat)
    print("itmat", la.det(itmat), itmat)

    spgList = []
    transStruc = []
    planes = [[1,0,0], [0,1,0], [0,0,1]]
    color_array = []
    size = 4
    Tpos = [] 
    for i in range(n_steps+1):
        curMat = (itmat-np.eye(3))*i/n_steps + np.eye(3)
        curStruc = Structure(curMat.dot(finalStruc.cell))
        for j,a in enumerate(dispStruc):
            curDisp = vec_classes[int(a.type)]*i/n_steps
            curPos = curMat.dot((finalStruc[j].pos - curDisp).reshape((3,1)))
            curStruc.add_atom(*(curPos.T.tolist()[0]),finalStruc[j].type)
            curStruc = supercell(curStruc, curStruc.cell)
        write.poscar(curStruc, vasp5=True, file=outdir+"/TransPOSCARS"+"/POSCAR_%03d"%i)
        spgList.append(get_spacegroup(to_spglib(curStruc), symprec=0.3, angle_tolerance=3.0))
        transStruc.append(curStruc)
        
        color_array.append([])
        Tpos.append([])
        for l,pl in enumerate(planes):
        
            plane = dispStruc.cell.dot(pl)
            tickness = 3 * la.norm(plane)
            plane = plane/tickness
            
            lattices = []
        
            Tpos_tmp = []
            types = []
            
            # Cut along the current plane
            for l in np.arange(-size,size):
                for j in np.arange(-size,size):
                    for k in np.arange(-size,size):
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
    ax.set_aspect('equal')
    fig.savefig(outdir+'/DispLattice_stretched_cell_primitive.svg')

    plt.show() #TMP

    # Growth directions

    stMat = rotate(tmat, np.eye(3)).dot(tmat)

    eigval, eigvec = la.eig(stMat)

    stretch_dir = la.inv(Bcell).dot(eigvec)

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


    # Showing some of the steps
    make_fig(0)
    make_fig(int(n_steps/5))
    make_fig(int(2*n_steps/5))
    make_fig(int(3*n_steps/5))
    make_fig(int(4*n_steps/5))
    make_fig(n_steps)

    if display:
        print("Showing")
        plt.show()
    else:
        print("make_anim(n_steps)")


if __name__=='__main__':
    p2ptrans(*readOptions())
