from p2ptrans import analysis, findMatching
import pytest
import os
import numpy as np
from test_units import read_FCC_BCC, BCC_file, FCC_file, tol, cleanup



def default_options():
    # (fileA, fileB, ncell, filename, interactive, savedisplay, outdir,
    # use, switch, prim, anim, vol, minimize, test, crystfile, n_steps,
    # showversion, map_ncell)
    return ('./POSCAR_A', './POSCAR_B', 300, './p2p.in', False, False, '.',
     False, False, True, False, False, False, False, './cryst.in', 60,
     False, None)

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
    
    print(tmat)
    print()

    # matricies can be in many diffrent forms...
    assert abs(tmat[0]) == pytest.approx([8.03921569e-01, 8.03921569e-01, 0], tol) or\
           abs(tmat[0]) == pytest.approx([0,  0, 8.03921569e-01]), tol
    assert abs(tmat[1]) == pytest.approx([8.03921569e-01, 8.03921569e-01, 0], tol) or\
           abs(tmat[1]) == pytest.approx([0,  0, 8.03921569e-01]), tol
    assert abs(tmat[2]) == pytest.approx([8.03921569e-01, 8.03921569e-01, 0], tol) or\
           abs(tmat[2]) == pytest.approx([0,  0, 8.03921569e-01]), tol
    assert np.linalg.det(tmat) == pytest.approx(1.0391327633537175, tol)
    
    assert abs(dispStruc.cell).flatten() == pytest.approx([1.435, 1.435, 1.435,
                                                           1.435, 1.435, 1.435,
                                                           1.435, 1.435, 1.435], tol)
    assert np.linalg.det(tmat) == pytest.approx(1.0391327619090698, tol)

    #TODO: what is this pattern?
    #assert vec_classes == pytest.approx([1.86087616e-08, -6.60150734e-09, 1.37320999e-08])
    # 1.15223042e-06,  8.17389458e-07, -8.06306018e-08
    # 1.86087616e-08, -6.60150734e-09, 1.37320999e-08
    # 8.92814178e-09,  3.42579574e-08, -2.57275090e-10
    # 7.52845120e-09, -4.41008163e-09,  2.90160900e-08

    assert dmin[0] == pytest.approx(6.65931827e+01, 0.1)
    assert dmin[1] == pytest.approx(1.67732607e-01, 0.1)
    assert dmin[2] == pytest.approx(-1.17534531e-02, 0.1) or \
        dmin[2] == pytest.approx(-9.84830562e-03, 0.1) # why is it like this?
    cleanup()

def test_crystallography():
    test_tmat = [[-8.03921569e-01, -8.03921569e-01, 0],
                 [ 8.03921569e-01, -8.03921569e-01, 0],
                 [ 0, 0, 8.03921569e-01]]
    ccell1, ccell2, planehkl, diruvw = analysis.readCrystParam('./FILE_DNE')
    BCC, FCC = read_FCC_BCC()
    eigval, U, P, Q, planeHab = analysis.crystallography(np.linalg.inv(test_tmat), FCC, BCC, ccell2, ccell1, planehkl,
                                                            diruvw, fileA=BCC_file, fileB=FCC_file)
    
    print(eigval)
    print(U)
    print(P)
    print(Q)
    print(planeHab)
    assert eigval == pytest.approx([0.87957185, 0.87957185, 1.24390244], tol)
    np.testing.assert_allclose(U,[[0.88275147, 0.05278822, 0.        ],
                [0.05278822, 0.87639223, 0.        ],
                [0.,         0.,         1.24390244]])
    np.testing.assert_allclose( P, [[1.,         0.06012458, 0.,       ],
                [0. ,        0.99819088, 0.,        ],
                [0.,         0.,         1.,        ]])
    np.testing.assert_allclose( Q, [[0.43187329, 1.,         0.,        ],
                [0.90193429, 0.,         0.,        ],
                [0.,         0.,         1.,        ]])
    np.testing.assert_allclose( planeHab, [[ 0., -0.],
    [ 0., -0.],
    [ 1., -1.]])

if __name__ == '__main__':
    test_crystallography()

   

    
