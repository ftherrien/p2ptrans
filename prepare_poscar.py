from pylada.crystal import read, supercell, write, Structure
import pickle
import matplotlib.pyplot as plt
import numpy as np
from p2ptrans import tiling as t
from copy import deepcopy

nlayers = 1

Alabel, Blabel, rtrans, ttrans, cell, centerA, centerB = pickle.load(open("transformation.dat","rb")) 

cell[:2,:2] = 2*cell[:2,:2]

centerA = np.append(centerA,0)
centerB = np.append(centerB,0)

A = read.poscar('/home/felixt/Recherche/p2ptrans/Surfaces/surfaces/ysz/poscars/POSCAR_111_2lay_10vac')
B = read.poscar('/home/felixt/Recherche/p2ptrans/Surfaces/surfaces/ni/poscars/POSCAR_111_3lay_10vac')

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

shiftA = np.ones(3)*1e10
for a in A:
    if a.type == Alabel:
        if a.pos[2] < shiftA[2]:
            shiftA = a.pos

for a in A:
    a.pos = a.pos - shiftA - centerA
    
Asupercell = cell
Asupercell[:,2] = A.cell[:,2]*nlayers

fig = plt.figure()
ax = fig.add_subplot(111)
ax.quiver(np.ones(2), np.ones(2), Asupercell[:2,0], Asupercell[:2,1], scale_units='xy', scale=1)
ax.quiver(np.ones(2), np.ones(2), A.cell[:2,0], A.cell[:2,1], scale_units='xy', scale=1, color = 'red')
maxXAxis = abs(Asupercell[:2,:2]).max() + 1
ax.set_xlim([-maxXAxis+1, maxXAxis+1])
ax.set_ylim([-maxXAxis+1, maxXAxis+1])
ax.set_aspect('equal')

ASC = t.circle(A.cell[:2,:2], 200)
BSC = t.circle(Asupercell[:2,:2], 10)

# Plot gamma points of each A cell
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

B = supercell(B, np.concatenate((B.cell[:,:2],-B.cell[:,2:3]), axis=1)) 

shiftB = np.ones(3)*-1e10
for b in B:
    if b.type == Blabel:
        if b.pos[2] > shiftB[2]:
            shiftB = b.pos

B.cell = ttrans[:,:3].dot(B.cell)
for b in B:
    b.pos = (ttrans[:,:3].dot((b.pos - shiftB - centerB).reshape(3,1)) + ttrans[:,3:4]).T[0]

B = supercell(B, B.cell)
    
B = supercell(B, np.concatenate((2*B.cell[:,:2],B.cell[:,2:3]), axis=1))
    
Bsupercell = cell
Bsupercell[:,2] = B.cell[:,2]*nlayers

Bsupercell_prev = deepcopy(Bsupercell)

Bsupercell = B.cell.dot(np.round(np.linalg.inv(B.cell).dot(Bsupercell)))

print("Bsupercell diff", Bsupercell - Bsupercell_prev)

B = supercell(B, Bsupercell)

ABcell = cell
ABcell[:,2] = A.cell[:,2]*2

fig = plt.figure()
ax = fig.add_subplot(111)
for a in A:
    if abs(a.pos[2]) < 1e-5:
        ax.scatter(a.pos[0], a.pos[1], color="red", alpha='0.5')

for b in B:
    if abs(b.pos[2]) < 1e-5:
        ax.scatter(b.pos[0], b.pos[1], color="blue", alpha='0.5')

plt.show()
        

AB = Structure(ABcell)
for a in A:
    AB.add_atom(*(tuple(a.pos-B.cell[:,2])+(a.type,)))
for a in B:
    AB.add_atom(*(tuple(a.pos-B.cell[:,2])+(a.type,)))

print(AB)
    
write.poscar(A, vasp5=True, file="POSCAR_A")
write.poscar(B, vasp5=True, file="POSCAR_B")
write.poscar(AB,vasp5=True, file="POSCAR_AB")

    


    
