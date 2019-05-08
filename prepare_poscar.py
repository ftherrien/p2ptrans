from pylada.crystal import read, supercell, write, Structure, primitive
import pickle
import matplotlib
# matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from p2ptrans import tiling as t
from copy import deepcopy
import sys

la = np.linalg

cut = False
spacing = 2 # Angstrom
tickness = 15
vacuum = 10 # on each side
substrate = "Zr"

Alabel, Blabel, rtrans, ttrans, cell, centerA, centerB = pickle.load(open("transformation_%s.dat"%(sys.argv[1]),"rb"))

if Blabel == substrate:
    tmp = Alabel
    Alabel = Blabel
    Blabel = tmp
    ttrans[:,:3] = la.inv(ttrans[:,:3])
    ttrans[:,3] = -ttrans[:,:3].dot(ttrans[:,3])
    cell = ttrans[:,:3].dot(cell)
    tmp = centerA
    centerA = centerB
    centerB = tmp

cell[:2,:2] = cell[:2,:2]

centerA = np.append(centerA,0)
centerB = np.append(centerB,0)

A = read.poscar('POSCAR_111_2lay_10vac')
B = read.poscar('POSCAR_111_3lay_10vac')

for a in A:
    if a.type == Alabel:
        break
else:
    tmp = A
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

oAcell = A.cell # original A.cell
oBcell = B.cell

B = primitive(B)
halfed_cell = oBcell
halfed_cell[:2,:2] = 0.5*halfed_cell[:2,:2]
B.cell = B.cell.dot(np.round(np.linalg.inv(B.cell).dot(halfed_cell)))
B = supercell(B, B.cell)

print("Alabel",Alabel,A)
print("Blabel",Blabel,B)

shiftA = -np.ones(3)*1e10
minA = 1e10
for a in A:
    if a.type == Alabel:
        if a.pos[2] > shiftA[2]:
            shiftA = a.pos
    if a.pos[2] < minA:
        minA = a.pos[2]

shiftB = np.ones(3)*1e10
maxB = -1e10
for b in B:
    if b.type == Blabel:
        if b.pos[2] < shiftB[2]:
            shiftB = b.pos
    if b.pos[2] > maxB:
        maxB = b.pos[2]

# fig = plt.figure()
# ax = fig.add_subplot(111)
# ax.set_aspect('equal')

# for a in supercell(A, np.array([[10,0,0],[0,10,0],[0,0,1]]).dot(A.cell):
#     if abs(a.pos[2] - shiftA) < 1e-5:
#         ax.scatter(a.pos[0], a.pos[1], color="red", alpha='0.5')

# for b in supercell(B, np.array([[10,0,0],[0,10,0],[0,0,1]]).dot(B.cell):
#     if abs(b.pos[2] - shiftB) < 1e-5:
#         ax.scatter(b.pos[0], b.pos[1], color="blue", alpha='0.5')

A.cell[2,2] = -(shiftA[2] - minA + vacuum) # The substrate is at the bottom
newA = Structure(A.cell)
for a in A:
    a.pos = a.pos - shiftA - centerA
    if a.pos[2] < 1e-10:
        newA.add_atom(*(tuple(a.pos)+(a.type,)))

A = deepcopy(newA) #TODO: remove would be more elegant, does not work for some reason

print("After shift (A):", len(A)) 

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

print("Should be", la.det(Asupercell) / la.det(A.cell))

A = supercell(A, Asupercell)

print("After supercell (A):", len(A))

B.cell = ttrans[:,:3].dot(B.cell)
print(maxB, shiftB[2], vacuum)
B.cell[2,2] = maxB - shiftB[2] + vacuum
for b in B:
    b.pos = (ttrans[:,:3].dot((b.pos - shiftB - centerB).reshape(3,1)) + ttrans[:,3:4]).T[0]

B = supercell(B, B.cell)

print("After shift (B):", len(B))

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

print("Should be", la.det(Bsupercell) / la.det(B.cell))

B = supercell(B, Bsupercell)

print("After supercell (B):", len(B))

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
ABcell[2,2] = B.cell[2,2] - A.cell[2,2] + spacing

print("ABcell", ABcell)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_aspect('equal')

once=False
for i,a in enumerate(A):
    if abs(a.pos[2]) < 1e-5:
        if not once:
            once = True
            ax.quiver(a.pos[0], a.pos[1], zrv1[0], zrv1[1], scale_units='xy', scale=1)
            ax.quiver(a.pos[0], a.pos[1], zrv2[0], zrv2[1], scale_units='xy', scale=1)

        if not np.mod(i,4):
            ax.scatter(a.pos[0], a.pos[1], color="red", alpha='0.5')
        else:
            ax.scatter(a.pos[0], a.pos[1], color="green", alpha='0.5')

once=False
for i, b in enumerate(B):
    if abs(b.pos[2]) < 1e-5:

        if not once:
            once = True
            center = b.pos
            ax.quiver(b.pos[0], b.pos[1], niv1[0], niv1[1], scale_units='xy', scale=1)
            ax.quiver(b.pos[0], b.pos[1], niv2[0], niv2[1], scale_units='xy', scale=1)

        print(abs(la.inv(oBcell).dot(b.pos - center) - np.round(la.inv(2*oBcell).dot(b.pos - center))))
        
        ax.scatter(b.pos[0], b.pos[1], color="blue", alpha='0.5')
            
AB = Structure(ABcell)
for a in A:
    AB.add_atom(*(tuple(a.pos-A.cell[:,2])+(a.type,)))
for a in B:
    AB.add_atom(*(tuple(a.pos-A.cell[:,2]+spacing)+(a.type,)))

# print("Finding primitive cell")
# A = primitive(supercell(A,A.cell), 1e-3)
# B = primitive(supercell(B,B.cell), 1e-3)
# AB = primitive(supercell(AB,AB.cell), 1e-3)

print("Writting poscars!")
print("A:", len(A))
print("B:", len(B))
print("AB:", len(AB))

write.poscar(A, vasp5=True, file="POSCAR_A_%s"%(sys.argv[1]))
write.poscar(B, vasp5=True, file="POSCAR_B_%s"%(sys.argv[1]))
write.poscar(AB,vasp5=True, file="relaxations/poscars/POSCAR_AB_%s"%(sys.argv[1]))

plt.show()

# Cutting surfaces (Assuming A is the substrate)

if cut:
    
    planes = np.array([[-1,1,0], [-1,-1,0]]).T
    
    ABcell = np.array([[2,0,0],[0,2,0],[0,0,1]]).dot(AB.cell)
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
