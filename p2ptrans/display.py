import matplotlib
from mpl_toolkits.mplot3d import Axes3D

from .config import *
from .utils import find_uvw, makeInitStruc, normal

plt=None
isset = False

def setplt(interactive):
    global isset
    global plt
    if not isset:
        if not interactive:
            matplotlib.use('Agg')
            
        import matplotlib.pyplot as plt
        
        plt.rcParams["figure.figsize"] = [5, 5]

        isset = True

def displayStats(stats, n_iter, peak_thetas, ttrans, dmin, n_peaks, sym, interactive, savedisplay, outdir):
    """ Displays the statistics of the interface minimization """

    setplt(interactive)

    angles = stats[:,0]
    angles= np.mod(angles,np.pi*2/sym)

    vol = stats[:,1]
        
    dists = stats[:,2]
        
    idx = np.argsort(angles)
    angles = 180*angles[idx]/np.pi
    vol = vol[idx]
    dists = dists[idx]

    size = int(0.02*n_iter)

    n_bins = n_iter//10
        
    fig = plt.figure("Stats", figsize=[15,7.5])

    gs = fig.add_gridspec(3, 6)

    ax = [None]*6
        
    ax[0] = fig.add_subplot(gs[0:2,0])
    ax[0].hist(dists,n_bins, orientation="horizontal")
    ax[0].set_xlim(ax[0].get_xlim()[::-1])
    ax[0].set_axis_off()
        
    ax[1] = fig.add_subplot(gs[0:2,1:3])
    ax[1].set_title("Distance vs angles")
    idvol = np.argsort(vol)[::-1]
    im = ax[1].scatter(angles[idvol], dists[idvol], c=vol[idvol])
        
    meanrun = np.zeros(n_iter-size)
    for i in range(n_iter-size):
        meanrun[i] = np.mean(dists[i:i+size])
            
    ax[1].plot(angles[size//2:-size//2],meanrun, 'gray')

    for i in range(n_peaks):
        ax[1].scatter(np.mod(180*peak_thetas[i]/np.pi,360/sym), dmin[i], marker="X", c="red")
            
    fig.colorbar(im, ax=ax[0:1], location='left')
        
    ax[2] = fig.add_subplot(gs[0:2,3:5])
    ax[2].set_title("Volume vs angles")
    iddist = np.argsort(dists)[::-1]
    im = ax[2].scatter(angles[iddist], vol[iddist], c=dists[iddist])

    for i in range(n_peaks):
        ax[2].scatter(np.mod(180*peak_thetas[i]/np.pi,360/sym), la.det(ttrans[i,:,:3]), marker="X", c="red")
        
    ax[3] = fig.add_subplot(gs[0:2,5])
    ax[3].hist(vol,n_bins, orientation="horizontal")
    ax[3].set_axis_off()
        
    fig.colorbar(im, ax=ax[3])

    ax[4] = fig.add_subplot(gs[2,1:3])
    ax[4].hist(angles,n_bins)
    ax[4].set_ylim(ax[4].get_ylim()[::-1])
    ax[4].set_axis_off()
        
    ax[5] = fig.add_subplot(gs[2,3:5])
    ax[5].hist(angles,n_bins)
    ax[5].set_ylim(ax[5].get_ylim()[::-1])
    ax[5].set_axis_off()

    if savedisplay:
        fig.savefig(outdir+"/stats.png")
        
    if interactive:
        plt.show()
        
def displayOptimalResult(Apos, Bpos, Bposst, disps_total, disps, class_list, vec_classes,
                         nat, natA, natB, atoms, outdir, savedisplay, interactive):

    """Displays or saves pictures of the optimal result"""
    
    setplt(interactive)
    
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

def displayTransCell(disps, dispStruc, foundcell,
                     pos_in_struc, vec_classes, outdir, interactive, savedisplay):   

    initStruc = makeInitStruc(dispStruc, vec_classes)
    
    cell = dispStruc.cell
    
    # Displays only the cell and the displacements in it
    fig = plt.figure("Transformation Cell", figsize = [10,5])
    
    ax = fig.add_subplot(121, projection='3d')
    ax.set_title('Transformation cell with displacements\n and atomic positions')
    for i,disp in enumerate(dispStruc):
        ax.quiver(disp.pos[0], disp.pos[1], disp.pos[2], vec_classes[int(disp.type)][0],vec_classes[int(disp.type)][1], vec_classes[int(disp.type)][2], color=colorlist[i%10])
        ax.scatter(disp.pos[0], disp.pos[1], disp.pos[2], alpha = 0.5, s=10, color=colorlist[i%10])
        ax.scatter(initStruc[i].pos[0], initStruc[i].pos[1], initStruc[i].pos[2], alpha = 1, s=10, color=colorlist[i%10])
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

def set_view(p,angle=0):
    p=p / la.norm(p)
    v1 = np.array([0,p[2],-p[1]])
    v1 = v1 / la.norm(v1)
    v2 = np.cross(p,v1)
    v2 = v2 / la.norm(v2)
    return la.inv(np.array([v1*np.cos(angle) + v2*np.sin(angle), -v1*np.sin(angle) + v2*np.cos(angle),
                     p]).T)

def dir2angles(plane):
    plane = plane/la.norm(plane)
    angles=np.zeros(2)
    a0 = np.arccos(plane[2])
    angles[0] = np.pi/2 - a0
    angles[1] = np.arctan2(plane[1], plane[0])
    return angles*180/np.pi

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

    apos1 = transStruc[state][np.argmin([la.norm(a.pos) for a in transStruc[state]])].pos
    a_list = []
    for a in transStruc[state]:
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
        
def make_fig(state, Tpos, color_array, transStruc, atom_types, spgList, outdir, savedisplay,
             interactive):
    setplt(interactive)
    fig = plt.figure(figsize=[7.2,7.2])
    gs = matplotlib.gridspec.GridSpec(2, 2)
    gs.update(wspace=0.03, hspace=0.03)
    all_panels(fig,gs, state, True, Tpos, color_array, transStruc, atom_types, spgList)
    if savedisplay:
        fig.savefig(outdir+"/Trans_%d.png"%state)
    if interactive:
        plt.show
        
def make_anim(n_states, Tpos, color_array, transStruc, atom_types, spgList, outdir):
    setplt(interactive)
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

def printMatAndDir(A,ccell):
    print("--------Matrix--------|-----Closest uvw------")
    print("    v1    v2    v3    |    d1    d2    d3    ")
    B = find_uvw(A, basis=ccell)
    A = la.inv(ccell).dot(A)
    for i, rowA in enumerate(A):
          rowB = B[i]
          print(' '.join(["% 5.3f"%(val) for val in rowA]), " |", ' '.join(["%5d"%(val) for val in rowB]))

