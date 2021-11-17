from p2ptrans import analysis, findMatching, produceTransition
import pytest
import os
import numpy as np
from pylada.crystal import Structure
from conftest import assert_structs_approx_eq

@pytest.fixture(scope="session")
def find_matching(bcc_fcc_s, bcc_fcc_filenames_s, double_cleanup_s, default_options_s):
    '''Runs a full matching, then cleanup the files'''
    (_, _, ncell, filename, interactive, savedisplay, outdir,
        use, switch, prim, anim, vol, minimize, test, crystfile, n_steps,
        showversion, map_ncell) = default_options_s
    bcc, fcc = bcc_fcc_s
    BCC_file, FCC_file = bcc_fcc_filenames_s
    
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    try:
        with open(filename, "r") as f:
            filecontent = f.readlines()
    except FileNotFoundError:
        filecontent = ""

    ccell1, ccell2, planehkl, diruvw = np.eye(3), np.eye(3), [1, 0, 0], [0, 1, 0]

    tmat, dispStruc, vec_classes, dmin = findMatching(bcc, fcc, ncell, fileA=BCC_file, fileB=FCC_file,
                                            ccellA=ccell1, ccellB=ccell2,
                                            filename=filename, interactive=interactive,
                                            savedisplay=savedisplay, outdir=outdir,
                                            switch=switch, prim=prim, vol=vol,
                                            minimize=minimize, test=test, map_ncell=map_ncell)
    
    return tmat, dispStruc, vec_classes, dmin

def test_tmat(find_matching):
    '''Test that the tmats are similar'''
    tmat, dispStruc, vec_classes, dmin = find_matching

    # matricies can be in many different forms...
    # TODO: row reduction or something?
    weird_tmat = abs(tmat)[abs(tmat)[:, 0].argsort()]
    print(tmat)
    print(weird_tmat)

    np.testing.assert_allclose(weird_tmat, [[0, 0, 0.8039216],
                                            [0.8039216, 0.8039216, 0],
                                            [0.8039216, 0.8039216, 0]], 0.001, atol=1e-14)
    #        [[0, 0.8039216, 0], # rare possiblity...
    #         [0.8039216, 0, 0.8039216],
    #         [0.8039216, 0, 0.8039216]]
    assert np.linalg.det(tmat) == pytest.approx(1.0391327633537175, 0.0001)

def test_dispCell(find_matching):
    '''Test that the dispCells are the same'''
    tmat, dispStruc, vec_classes, dmin = find_matching
    print(dispStruc)
    np.testing.assert_allclose(abs(dispStruc.cell), [[1.435, 1.435, 1.435],
                                                     [1.435, 1.435, 1.435],
                                                     [1.435, 1.435, 1.435]], 0.0001)
    assert np.linalg.det(tmat) == pytest.approx(1.0391327619090698, 0.0001)    

@pytest.mark.skip(reason="No idea how to test this...")
def test_vec_classes(find_matching):
    tmat, dispStruc, vec_classes, dmin = find_matching 
    print(vec_classes)
    #TODO: what is this pattern?
    #assert vec_classes == pytest.approx([1.86087616e-08, -6.60150734e-09, 1.37320999e-08])
    # 1.15223042e-06,  8.17389458e-07, -8.06306018e-08
    # 1.86087616e-08, -6.60150734e-09, 1.37320999e-08
    # 8.92814178e-09,  3.42579574e-08, -2.57275090e-10
    # 7.52845120e-09, -4.41008163e-09,  2.90160900e-08   

def test_dmin_12(find_matching):
    '''Test that the first 2 dmins match expected'''
    tmat, dispStruc, vec_classes, dmin = find_matching  

    assert dmin[0] == pytest.approx(6.65931827e+01, 0.1)
    assert dmin[1] == pytest.approx(1.67732607e-01, 0.1)

@pytest.mark.skip(reason="What does this mean?")
def test_dmin_3(find_matching):
    '''Test that the third dmin is in the list'''
    tmat, dispStruc, vec_classes, dmin = find_matching 
    assert dmin[2] == pytest.approx(-1.17534531e-02, 0.1) or \
            dmin[2] == pytest.approx(-9.84830562e-03, 0.1) or\
            dmin[2] == pytest.approx( 0.00861095521, 0.1) or\
            dmin[2] == pytest.approx(-0.00826468673, 0.1) or\
            dmin[2] == pytest.approx(-0.01332142476, 0.1) # why is it like this?


def test_crystallography(bcc_fcc, bcc_fcc_filenames):
    '''Test crystallography '''
    test_tmat = [[-8.03921569e-01, -8.03921569e-01, 0],
                 [ 8.03921569e-01, -8.03921569e-01, 0],
                 [ 0, 0, 8.03921569e-01]]
    bcc, fcc = bcc_fcc
    BCC_file, FCC_file = bcc_fcc_filenames
    ccell1, ccell2, planehkl, diruvw = analysis.readCrystParam('./FILE_DNE')
    eigval, U, P, Q, planeHab = analysis.crystallography(np.linalg.inv(test_tmat), fcc, bcc, ccell2, ccell1, planehkl,
                                                            diruvw, fileA=BCC_file, fileB=FCC_file)
    
    assert eigval == pytest.approx([0.87957185, 0.87957185, 1.24390244], 0.0001)
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

@pytest.mark.usefixtures("double_cleanup")
def test_transisiton():
    test_tmat = [[-8.03921569e-01, -8.03921569e-01, 0],
                 [ 8.03921569e-01, -8.03921569e-01, 0],
                 [ 0, 0, 8.03921569e-01]]
    test_disp = Structure([[ 1.435, -1.435,  1.435],
                           [ 1.435, -1.435, -1.435],
                           [ 1.435,  1.435, -1.435]])\
    .add_atom(0.09, -2.74, 0.014, '0', atom='Fe', site=0)
    vec_classes = [np.array([ 9.40440525e-09, -3.27675309e-09,  4.11677981e-08])]
    result = produceTransition(5, test_tmat, test_disp, vec_classes,
                               '.', False,
                               a_name='BCC', b_name='FCC')
    path = result[0]
    spacegroups = result[1]

    assert len(path) == 6
    assert (spacegroups == ['Im-3m (229)', 'I4/mmm (139)', 'I4/mmm (139)', 'I4/mmm (139)', 'I4/mmm (139)', 'Fm-3m (225)'])

    res0 = Structure([[ 1.435, -1.435,  1.435],
                      [ 1.435, -1.435, -1.435],
                      [ 1.435,  1.435, -1.435]])\
        .add_atom(0.09000, -2.74000, 0.014000, 'Fe')
    res2 = Structure([[ 1.318706, -1.413042,  1.410834],
                      [ 1.413042, -1.318706, -1.410834],
                      [ 1.571877,  1.571877, -1.505171]])\
        .add_atom(0.083026, -2.60537, 0.081111, 'Fe')
    res5 = Structure([[ 1.144264, -1.380106,  1.374586],
                      [ 1.380106, -1.144265, -1.374586],
                      [ 1.777193,  1.77719283, -1.610428]])\
        .add_atom(0.072567, -2.40343, 0.18178, 'Fe')

    assert_structs_approx_eq(path[0], res0)
    assert_structs_approx_eq(path[2], res2)
    assert_structs_approx_eq(path[5], res5)


if __name__ == '__main__':
    test_crystallography()

   

    
