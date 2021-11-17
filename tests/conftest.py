from pylada.crystal import Structure
from p2ptrans import read
import os, glob, shutil
import numpy as np
import pytest

BCC_file = './BCC_POSCAR'
FCC_file = './FCC_POSCAR'

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
    return _bcc, _fcc

@pytest.fixture()
def bcc_fcc():
    return _bcc, _fcc

@pytest.fixture()
def bcc_fcc_filenames():
    return BCC_file, FCC_file

@pytest.fixture(scope='session')
def bcc_fcc_filenames_s():
    return BCC_file, FCC_file



@pytest.fixture(scope="session")
def default_options_s():
    # (fileA, fileB, ncell, filename, interactive, savedisplay, outdir,
    # use, switch, prim, anim, vol, minimize, test, crystfile, n_steps,
    # showversion, map_ncell)
    return ('./POSCAR_A', './POSCAR_B', 300, './p2p.in', False, False, '.',
        False, False, True, False, False, False, False, './cryst.in', 60,
        False, None)

def assert_structs_approx_eq(A, B, tol=0.001):    
    assert len(A) == len(B)
    np.testing.assert_allclose(A.cell, B.cell, tol)
    assert A.scale == pytest.approx(B.scale)
    for aAtom, bAtom in zip(A, B):
        assert aAtom.type == bAtom.type
        assert aAtom.pos == pytest.approx(bAtom.pos, tol)


def cleanup():
    for dat in glob.iglob('*dat'):
        os.remove(dat)
    if os.path.exists('progress.txt'):
        os.remove('progress.txt')
    if os.path.exists('TransPOSCARS'):
        shutil.rmtree('TransPOSCARS')

@pytest.fixture()
def double_cleanup():
    cleanup()
    yield
    cleanup()

@pytest.fixture(scope="session")
def double_cleanup_s():
    cleanup()
    yield
    cleanup()
