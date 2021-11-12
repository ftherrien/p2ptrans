from p2ptrans import *
import os


BCC_file = './BCC_POSCAR'
FCC_file = './FCC_POSCAR'

def default_options():
    # (fileA, fileB, ncell, filename, interactive, savedisplay, outdir,
    # use, switch, prim, anim, vol, minimize, test, crystfile, n_steps,
    # showversion, map_ncell)
    return ('./POSCAR_A', './POSCAR_B', 300, './p2p.in', False, False, '.',
     False, False, True, False, False, False, False, './cryst.in', 60,
     False, None)

def read_FCC_BCC():
    BCC = read.poscar(BCC_file)
    FCC = read.poscar(FCC_file)
    return BCC, FCC

def test_read_poscars():
    '''Test that pylada reads test files correctly'''
    (BCC, FCC) = read_FCC_BCC()

    assert BCC.scale == 2.87
    assert (BCC.cell == [[0.5,-0.5, 0.5],
                        [ 0.5, 0.5,-0.5],
                        [-0.5, 0.5, 0.5]]).all
    assert (BCC[0].pos == [0, 0,-0]).all
    
    assert FCC.scale == 3.57
    assert (FCC.cell == [[0.5, 0.5, 0.0],
                        [ 0.5, 0.0, 0.5],
                        [ 0.0, 0.5, 0.5]]).all
    assert (BCC[0].pos == [0, 0, 0]).all

    assert BCC[0].type == FCC[0].type

def test_read_cryst_defaults():
    ccell1, ccell2, planehkl, diruvw = analysis.readCrystParam('./FILE_DNE')
    assert (ccell1 == [[1., 0., 0.],
                       [0., 1., 0.],
                       [0., 0., 1.]]).all
    assert (ccell2 == [[1., 0., 0.],
                       [0., 1., 0.],
                       [0., 0., 1.]]).all
    assert planehkl == [1, 0, 0]
    assert diruvw == [0, 1, 0]
    
