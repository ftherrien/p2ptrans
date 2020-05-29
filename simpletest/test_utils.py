from os import path
from pylada.crystal import read, Structure
import numpy as np
import itertools

def poscars_in_dir(main_dir):
    '''returns a sorted list of all poscars in a dir in the format POSCAR_XXX sequentially'''
    poscars = []
    for i in itertools.count():
        poscar_path = path.join(main_dir, 'POSCAR_{:03}'.format(i))
        if path.exists(poscar_path):
            poscars.append(poscar_path)
        else:
            return poscars

def structs_in_dir(main_dir):
    '''returns a sorted list of pylada structures of poscars in a dir of the format POSCAR_XXX'''
    return [read.poscar(poscar) for poscar in poscars_in_dir(main_dir)]

def compare_structs(struct_a, struct_b, tol=1e-4):
    '''NOT INTENDED TO BE CORRECT, just a quick check
       Will produce false positives (cell) and false negitives (supercells, transposed atoms)'''
    cells = compare_cells(struct_a, struct_b)
    num_atoms = len(struct_a) == len(struct_b)
    atom_poses = all([np.allclose(a.pos, b.pos, atol=tol) for a, b in zip(struct_a, struct_b)])
    atom_types = all([a.type == b.type for a, b in zip(struct_a, struct_b)])
    return cells and num_atoms and atom_poses and atom_types

def compare_cells(struct_a, struct_b, tol=1e-4):
    '''NOT INTENDED TO BE CORRECT, just a quick check
       Will produce false positives (doesn't check angles)'''
    cell_vols = np.isclose(struct_a.volume, struct_b.volume, atol=tol)
    cell_vectors_a_lens = np.linalg.norm(struct_a.cell, axis=0)
    cell_vectors_b_lens = np.linalg.norm(struct_b.cell, axis=0)
    cell_vectors_a_lens = np.sort(cell_vectors_a_lens)
    cell_vectors_b_lens = np.sort(cell_vectors_b_lens)
    return np.allclose(cell_vectors_a_lens, cell_vectors_b_lens, atol=tol) and cell_vols


def compare_tmats(tmat_a, tmat_b, tol=1e-4):
    '''just a quick check'''
    norma, normb = np.linalg.norm(tmat_a, axis=1), np.linalg.norm(tmat_b, axis=1)
    norma, normb = np.sort(norma), np.sort(normb)
    return np.allclose(norma, normb, atol=tol)




