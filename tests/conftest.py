from pylada.crystal import Structure
from p2ptrans import read
import os, glob, shutil
import numpy as np
import pytest

BCC_filename = './BCC_POSCAR'
FCC_filename = './FCC_POSCAR'
testmat_filename = './testmat.txt'
testmat_map_filename = './testmat_map.txt'

_bcc = Structure([[0.5,-0.5, 0.5],
                  [ 0.5, 0.5,-0.5],
                  [-0.5, 0.5, 0.5]], scale=2.87, name='BCC structure')\
           .add_atom(0., 0., 0., 'Fe')

_fcc = Structure([[0.5, 0.5, 0.0],
                  [ 0.5, 0.0, 0.5],
                  [ 0.0, 0.5, 0.5]], scale=3.57, name='FCC structure')\
           .add_atom(0., 0., 0., 'Fe')

@pytest.fixture(scope="session")
def bcc_fcc_s():
    '''Fixture that returns the test fcc and bcc structures (session scope)'''
    return _bcc, _fcc

@pytest.fixture()
def bcc_fcc():
    '''Fixture that returns the test fcc and bcc structures'''
    return _bcc, _fcc

@pytest.fixture()
def bcc_fcc_filenames():
    '''Fixture that returns the fcc and bcc filenames'''
    return BCC_filename, FCC_filename

@pytest.fixture(scope='session')
def bcc_fcc_filenames_s():
    '''Fixture that returns the fcc and bcc filenames (session scope)'''
    return BCC_filename, FCC_filename



@pytest.fixture(scope="session")
def p2p_default_options_s():
    '''Fixture that returns the defualt options from readOptions in the runner'''
    # (fileA, fileB, ncell, filename, interactive, savedisplay, outdir,
    # use, switch, prim, anim, vol, minimize, test, crystfile, n_steps,
    # showversion, map_ncell)
    return ('./POSCAR_A', './POSCAR_B', 300, './p2p.in', False, False, '.',
        False, False, True, False, False, False, False, './cryst.in', 60,
        False, None)

def p2p_int_default_options():
    '''Fixture that returns the defualt p2pint options'''
    # (fileA, fileB, planeA, planeB,
    # ruleA, ruleB, ncell, n_iter, filename,
    # interactive, savedisplay, outdir, use, minimize, test, n_iter, sym,
    # vacuum, layers, surface, max_thickness, showversion)
    return [['./POSCAR_A'], ['./POSCAR_B'], [0, 0, 1], [1, 1, 0],
            {'1': {'Si'}}, {'1': {'Si'}}, 100, 1000, './p2p.in',
            False, False, '.', False, False, False, 1000, 1,
            10, 1, False, None, False]


@pytest.fixture()
def lap_test_matrix():
    return np.loadtxt(testmat_filename)

@pytest.fixture()
def lap_test_map():
    return np.loadtxt(testmat_map_filename)



def assert_structs_approx_eq(A, B, tol=0.001):
    '''Assert that the two structures have the same cell, scales, and atoms of the same type in the same positions'''
    assert len(A) == len(B)
    np.testing.assert_allclose(A.cell, B.cell, tol)
    assert A.scale == pytest.approx(B.scale)
    for aAtom, bAtom in zip(A, B):
        assert aAtom.type == bAtom.type
        assert aAtom.pos == pytest.approx(bAtom.pos, tol)


def cleanup():
    '''Removes dat files, the progress file, and tranPOSCARS from the current directory'''
    for dat in glob.iglob('*dat'):
        os.remove(dat)
    if os.path.exists('progress.txt'):
        os.remove('progress.txt')
    if os.path.exists('out.txt'):
        os.remove('out.txt')
    if os.path.exists('TransPOSCARS'):
        shutil.rmtree('TransPOSCARS')
    if os.path.exists('DC_C_POSCAR-DC_Si_POSCAR'):
        shutil.rmtree('DC_C_POSCAR-DC_Si_POSCAR')
    

@pytest.fixture()
def double_cleanup():
    '''Fixture that cleans up both before and after'''
    cleanup()
    yield
    cleanup()

@pytest.fixture(scope="session")
def double_cleanup_s():
    '''Fixture that cleans up both before and after (session scope)'''
    cleanup()
    yield
    cleanup()
