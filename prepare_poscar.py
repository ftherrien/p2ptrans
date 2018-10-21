from pylada.crystal import read, supercell, write, Structure
import pickle
import matplotlib.pyplot as plt
import numpy as np
from p2ptrans import tiling as t
from copy import deepcopy

la = np.linalg

spacing = 2 # Angstrom

tickness = 10

Alabel, Blabel, rtrans, ttrans, cell, centerA, centerB = pickle.load(open("transformation.dat","rb")) 

cell[:2,:2] = 2*cell[:2,:2]

centerA = np.append(centerA,0)
centerB = np.append(centerB,0)

A = read.poscar('/home/felixt/Recherche/p2ptrans/Surfaces/surfaces/YSZ/poscars/POSCAR_111_2lay_10vac')
B = read.poscar('/home/felixt/Recherche/p2ptrans/Surfaces/surfaces/ni/poscars/POSCAR_111_3lay_10vac')

oAcell = A.cell # original A.cell
oBcell = B.cell

for a in A:
    if a.type == Alabel:
        break
else:
    A = tmp
    A = B
    B = tmp
    for a in A:
        if a.type == Alabel:
            break
    else:
        raise RuntimeError('Wrong poscars')

for b in B:
    if b.type == Blabel:
        break
else:
    raise RuntimeError('Wrong poscars')

shiftA = -np.ones(3)*1e10
for a in A:
    if a.type == Alabel:
        if a.pos[2] > shiftA[2]:
            shiftA = a.pos

shiftB = np.ones(3)*1e10
for b in B:
    if b.type == Blabel:
        if b.pos[2] < shiftB[2]:
            shiftB = b.pos

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_aspect('equal')

for a in supercell(A, np.array([[10,0,0],[0,10,0],[0,0,1]]).dot(A.cell):
    if abs(a.pos[2] - shiftA) < 1e-5:
        ax.scatter(a.pos[0], a.pos[1], color="red", alpha='0.5')

for b in supercell(B, np.array([[10,0,0],[0,10,0],[0,0,1]]).dot(B.cell):
    if abs(b.pos[2] - shiftB) < 1e-5:
        ax.scatter(b.pos[0], b.pos[1], color="blue", alpha='0.5')

A.cell[:,2] = -A.cell[:,2] # The substrate is at the bottom
newA = Structure(A.cell)
for a in A:
    a.pos = a.pos - shiftA - centerA
    if a.pos[2] < 1e-10:
        newA.add_atom(*(tuple(a.pos)+(a.type,)))

A = deepcopy(newA) #TODO: remove would be more elegant, does not work for some reason

Asupercell = cell
Asupercell[:,2] = A.cell[:,2]

# Visually checks that the common cell is a supercell of A
ASC = t.circle(A.cell[:2,:2], 200)
BSC = t.circle(Asupercell[:2,:2], 10)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.scatter(ASC[0,:], ASC[1,:], c='k')
maxXAxis = np.max([ASC[0,:].max(),ASC[1,:].max()])
minXAxis = np.min([ASC[0,:].min(),ASC[1,:].min()])
ax.set_xlim([minXAxis-1, maxXAxis+1])
ax.set_ylim([minXAxis-1, maxXAxis+1])
ax.set_aspect('equal')
ax.scatter(BSC[0,:], BSC[1,:], facecolor='w', edgecolor='k')

Asupercell_prev = deepcopy(Asupercell)

Asupercell = A.cell.dot(np.round(np.linalg.inv(A.cell).dot(Asupercell)))

A = supercell(A, Asupercell)

B.cell = ttrans[:,:3].dot(B.cell)
for b in B:
    b.pos = (ttrans[:,:3].dot((b.pos - shiftB - centerB).reshape(3,1)) + ttrans[:,3:4]).T[0]

B = supercell(B, B.cell)
    
# B = supercell(B, np.concatenate((2*B.cell[:,:2],B.cell[:,2:3]), axis=1))
    
Bsupercell = cell
Bsupercell[:,2] = B.cell[:,2]

# Check that the common cell is a supercell of A
ASC = t.circle(B.cell[:2,:2], 200)
BSC = t.circle(Bsupercell[:2,:2], 10)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.scatter(ASC[0,:], ASC[1,:], c='k')
maxXAxis = np.max([ASC[0,:].max(),ASC[1,:].max()])
minXAxis = np.min([ASC[0,:].min(),ASC[1,:].min()])
ax.set_xlim([minXAxis-1, maxXAxis+1])
ax.set_ylim([minXAxis-1, maxXAxis+1])
ax.set_aspect('equal')
ax.scatter(BSC[0,:], BSC[1,:], facecolor='w', edgecolor='k')

Bsupercell_prev = deepcopy(Bsupercell)

Bsupercell = B.cell.dot(np.round(np.linalg.inv(B.cell).dot(Bsupercell)))

print("Bsupercell diff", Bsupercell - Bsupercell_prev)

B = supercell(B, Bsupercell)

# Angles between main planes

niv1 = ttrans[:,:3].dot(oBcell).dot([1,-1,0])
zrv1 = oAcell.dot([1,-1,0])

a1 = np.arccos(niv1.dot(zrv1)/(la.norm(niv1)*la.norm(zrv1)))
print("The angle between Ni[-1,1,0] and YSZ[1-10] is:", 360*a1/(2*np.pi))

niv2 = ttrans[:,:3].dot(oBcell).dot([1,1,0])
zrv2 = oAcell.dot([1,1,0])

a2 = np.arccos(niv2.dot(zrv2)/(la.norm(niv2)*la.norm(zrv2)))
print("The angle between Ni[1,1,-2] and YSZ[11-2] is:", 360*a1/(2*np.pi))

ABcell = cell
ABcell[:,2] = B.cell[:,2]*2

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_aspect('equal')

once=False
for a in A:
    if abs(a.pos[2]) < 1e-5:
        ax.scatter(a.pos[0], a.pos[1], color="red", alpha='0.5')
        if not once:
            once = True
            ax.quiver(a.pos[0], a.pos[1], zrv1[0], zrv1[1], scale_units='xy', scale=1)
            ax.quiver(a.pos[0], a.pos[1], zrv2[0], zrv2[1], scale_units='xy', scale=1)

once=False
for b in B:
    if abs(b.pos[2]) < 1e-5:
        ax.scatter(b.pos[0], b.pos[1], color="blue", alpha='0.5')
        if not once:
            once = True
            ax.quiver(b.pos[0], b.pos[1], niv1[0], niv1[1], scale_units='xy', scale=1)
            ax.quiver(b.pos[0], b.pos[1], niv2[0], niv2[1], scale_units='xy', scale=1)

plt.show()
            
AB = Structure(ABcell)
for a in A:
    AB.add_atom(*(tuple(a.pos+B.cell[:,2])+(a.type,)))
for a in B:
    AB.add_atom(*(tuple(a.pos+B.cell[:,2]+spacing)+(a.type,)))

A = supercell(A,A.cell)
B = supercell(B,B.cell)
AB = supercell(AB,AB.cell)
    
write.poscar(A, vasp5=True, file="POSCAR_A")
write.poscar(B, vasp5=True, file="POSCAR_B")
write.poscar(AB,vasp5=True, file="POSCAR_AB")

# Cutting surfaces (Assuming A is the substrate)

planes = np.array([[-1,1,0], [-1,-1,0]]).T

ABcell = np.array([[4,0,0],[0,4,0],[0,0,1]]).dot(AB.cell)
cm = np.sum(ABcell/2, axis=1)

tol = 1e-5

for i,p in enumerate(ttrans[:,:3].dot(oBcell).dot(planes).T):

    p = p / np.linalg.norm(p)

    struc = Structure(ABcell)

    for a in supercell(AB, ABcell):
        if p.dot(a.pos) < p.dot(cm) + tol and p.dot(a.pos) > p.dot(cm) - tickness:
            struc.add_atom(*(tuple(a.pos)+(a.type,)))

    print("Writing POSCAR for plane %s"%"-".join([str(e) for e in planes[:,i]]))
    print("Oui?")
    write.poscar(struc, vasp5=True, file="POSCAR_cut%s"%"-".join([str(e) for e in planes[:,i]]))
    print("Non?")

    print("plotting from top", len(struc))
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_aspect('equal')
    for a in struc:
        if a.type == Alabel:
            ax.scatter(a.pos[0], a.pos[1], color="red", alpha='0.5')
        if a.type == Blabel:
            ax.scatter(a.pos[0], a.pos[1], color="green", alpha='0.5')

    print("plotting from the side")
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_aspect('equal')
    for a in struc:
        if a.type == Alabel:
            ax.scatter(a.pos.dot(np.cross(p,np.array([0,0,1]))), a.pos[2], color="red", alpha='0.5')
        if a.type == Blabel:
            ax.scatter(a.pos.dot(np.cross(p,np.array([0,0,1]))), a.pos[2], color="green", alpha='0.5')
    
plt.show()
