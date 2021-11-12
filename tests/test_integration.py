from _pytest.python_api import approx
from p2ptrans import *
import pytest
import os
import glob
import numpy as np


BCC_file = './BCC_POSCAR'
FCC_file = './FCC_POSCAR'
atol = 0.0001

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

#TODO test actual read_cryst    

def cleanup():
    for dat in glob.iglob('*dat'):
        os.remove(dat)
    if os.path.exists('progress.txt'):
        os.remove('progress.txt')

def test_matching():
    '''Runs a full matching, then tests that the tmats, dispCell, and first two dmins are the same'''
    cleanup()
    (_, _, ncell, filename, interactive, savedisplay, outdir,
        use, switch, prim, anim, vol, minimize, test, crystfile, n_steps,
        showversion, map_ncell) = default_options()
    
    BCC, FCC = read_FCC_BCC()

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    try:
        with open(filename, "r") as f:
            filecontent = f.readlines()
    except FileNotFoundError:
        filecontent = ""

    ccell1, ccell2, planehkl, diruvw = analysis.readCrystParam('./FILE_DNE')

    tmat, dispStruc, vec_classes, dmin = findMatching(BCC, FCC, ncell, fileA=BCC_file, fileB=FCC_file,
                                            ccellA=ccell1, ccellB=ccell2,
                                            filename=filename, interactive=interactive,
                                            savedisplay=savedisplay, outdir=outdir,
                                            switch=switch, prim=prim, vol=vol,
                                            minimize=minimize, test=test, map_ncell=map_ncell)
    

    # matricies can be in many diffrent forms...
    assert (abs(tmat[0]) == pytest.approx([8.03921569e-01, 8.03921569e-01, 0], atol) or
           abs(tmat[0]) == pytest.approx([0,  0, 8.03921569e-01]), atol)
    assert (abs(tmat[1]) == pytest.approx([8.03921569e-01, 8.03921569e-01, 0], atol) or
           abs(tmat[1]) == pytest.approx([0,  0, 8.03921569e-01]), atol)
    assert (abs(tmat[2]) == pytest.approx([8.03921569e-01, 8.03921569e-01, 0], atol) or
           abs(tmat[2]) == pytest.approx([0,  0, 8.03921569e-01]), atol)
    assert np.linalg.det(tmat) == pytest.approx(1.0391327633537175, atol)
    
    assert abs(dispStruc.cell).flatten() == pytest.approx([1.435, 1.435, 1.435,
                                                           1.435, 1.435, 1.435,
                                                           1.435, 1.435, 1.435], atol)
    assert np.linalg.det(tmat) == pytest.approx(1.0391327619090698, atol)

    #TODO: what is this pattern?
    #assert vec_classes == pytest.approx([1.86087616e-08, -6.60150734e-09, 1.37320999e-08])
    # 1.15223042e-06,  8.17389458e-07, -8.06306018e-08
    # 1.86087616e-08, -6.60150734e-09, 1.37320999e-08
    # 8.92814178e-09,  3.42579574e-08, -2.57275090e-10
    # 7.52845120e-09, -4.41008163e-09,  2.90160900e-08

    assert dmin[0] == pytest.approx(6.65931827e+01, 0.1)
    assert dmin[1] == pytest.approx(1.67732607e-01, 0.1)
    #TODO: what is the third dmin???
    # [ 6.65931827e+01  1.67732607e-01 -1.17534531e-02]
    # [ 6.65349152e+01  1.67163017e-01 -9.84830562e-03]
    # [ 6.65388928e+01  1.67598176e-01 -1.17542355e-02]
    cleanup()

