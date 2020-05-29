#!/usr/bin/env python
from os import path
import p2ptrans as p2p
from pylada.crystal import read
import itertools
import numpy as np

bcc = read.poscar(path.join('ref','BCC_Fe_POSCAR'))
fcc = read.poscar(path.join('ref','FCC_Fe_POSCAR'))

outdir = 'results'

# n_steps, (tmat, dispStruc, vec_classes), outdir, display
tmat, dispStruc, vec_classes = p2p.findMatching(bcc, fcc, 100, outdir=outdir)
p2p.produceTransition(60, tmat, dispStruc, vec_classes, outdir, False)

def poscars_in_dir(main_dir):
    poscars = []
    for i in itertools.count():
        poscar_path = path.join(main_dir, 'POSCAR_{:03}'.format(i))
        if path.exists(poscar_path):
            poscars.append(poscar_path)
        else:
            return poscars

def structs_in_dir(main_dir):
    return [read.poscar(poscar) for poscar in poscars_in_dir(main_dir)]

def compare_structs(struct_a, struct_b, tol=1e-4):
    cells = compare_cells(struct_a, struct_b)
    num_atoms = len(struct_a) == len(struct_b)
    atom_poses = all([np.allclose(a.pos, b.pos, atol=tol) for a, b in zip(struct_a, struct_b)])
    atom_types = all([a.type == b.type for a, b in zip(struct_a, struct_b)])
    return cells and num_atoms and atom_poses and atom_types
    
def compare_cells(struct_a, struct_b, tol=1e-4):
    cell_vols = np.isclose(struct_a.volume, struct_b.volume, atol=tol)
    cell_vectors_a_lens = np.linalg.norm(struct_a.cell, axis=0)
    cell_vectors_b_lens = np.linalg.norm(struct_b.cell, axis=0)
    cell_vectors_a_lens = np.sort(cell_vectors_a_lens)
    cell_vectors_b_lens = np.sort(cell_vectors_b_lens)
    return np.allclose(cell_vectors_a_lens, cell_vectors_b_lens, atol=tol) and cell_vols 

# then verifiy

ref_structs = structs_in_dir(path.join('ref', 'TransPOSCARS'))
produced_structs = structs_in_dir(path.join('results', 'TransPOSCARS'))

for ref_struct, res_struct in zip(ref_structs, produced_structs):
    if not compare_structs(ref_struct, res_struct, tol=0.1):
        print(ref_struct, 'DNE', res_struct)
        exit()

print('Succeded!')




