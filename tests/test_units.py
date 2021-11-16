from p2ptrans import read, analysis, core
from pylada.crystal import Structure
import numpy as np
import pytest
import os, glob, shutil


BCC_file = './BCC_POSCAR'
FCC_file = './FCC_POSCAR'
tol = 0.0001

def assert_structs_approx_eq(A, B, tol=0.001):    
    assert len(A) == len(B)
    np.testing.assert_allclose(A.cell, B.cell, tol)
    assert A.scale == pytest.approx(B.scale)
    for aAtom, bAtom in zip(A, B):
        assert aAtom.type == bAtom.type
        assert aAtom.pos == pytest.approx(bAtom.pos, tol)

def read_FCC_BCC():
    BCC = read.poscar(BCC_file)
    FCC = read.poscar(FCC_file)
    return BCC, FCC    

def cleanup():
    for dat in glob.iglob('*dat'):
        os.remove(dat)
    if os.path.exists('progress.txt'):
        os.remove('progress.txt')
    if os.path.exists('TransPOSCARS'):
        shutil.rmtree('TransPOSCARS')

def test_read_poscars():
    '''Test that pylada reads test files correctly'''
    (BCC, FCC) = read_FCC_BCC()
    test_BCC = Structure([[0.5,-0.5, 0.5],
                        [ 0.5, 0.5,-0.5],
                        [-0.5, 0.5, 0.5]], scale=2.87)\
               .add_atom(0., 0., 0., 'Fe')
    test_FCC = Structure([[0.5, 0.5, 0.0],
                        [ 0.5, 0.0, 0.5],
                        [ 0.0, 0.5, 0.5]], scale=3.57)\
               .add_atom(0., 0., 0., 'Fe')

    assert_structs_approx_eq(BCC, test_BCC)
    assert_structs_approx_eq(FCC, test_FCC)

def test_read_cryst_defaults():
    '''Test that the cryst param defaults are the unit matrix'''
    ccell1, ccell2, planehkl, diruvw = analysis.readCrystParam('./FILE_DNE')
    assert (ccell1 == [[1., 0., 0.],
                       [0., 1., 0.],
                       [0., 0., 1.]]).all
    assert (ccell2 == [[1., 0., 0.],
                       [0., 1., 0.],
                       [0., 0., 1.]]).all
    assert planehkl == [1, 0, 0]
    assert diruvw == [0, 1, 0]

def test_read_cryst():
    '''Tests that the cryst param file is read correctly'''
    ccell1, ccell2, planehkl, diruvw = analysis.readCrystParam('./cryst_test_file')
    assert ccell1 == [[1., 2., 3.],[4., 5., 6.],[7., 8., 9.]]
    assert ccell2 == [[9., 8., 7.],[6., 5., 4.],[3., 2., 1.]]
    assert planehkl == [1, 2, 3]
    assert diruvw == [3, 2, 1]

@pytest.mark.skip(reason="Replaced with smaller and larger tests...")
def test_optimize():
    cleanup()
    BCC, FCC = read_FCC_BCC()
    BCC_cell = BCC.cell * 2.87
    FCC_cell = FCC.cell * 3.57
    # optimization(A, Acell, mulA, B, Bcell, mulB, ncell, filename, outdir, max_cell_size)
    result = core.optimization(BCC, BCC_cell, 1, FCC, FCC_cell, 1, 300, './p2p.in', '.', 1000)
    Apos, Apos_map, Bpos, Bposst, n_map, natA, class_list, tmat, dmin, atoms, atom_types, foundcell, vec = result
    print(Apos, Apos_map, Bpos, Bposst, n_map, natA, class_list, tmat, dmin, atoms, atom_types, foundcell, vec)
    cleanup()
