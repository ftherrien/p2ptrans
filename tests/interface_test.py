from p2ptrans.interfaces import findMatchingInterfaces
from conftest import assert_structs_approx_eq
from pylada.crystal import Structure
import numpy as np
import pytest
import os


@pytest.fixture(scope="module")
def find_matching_interface(double_cleanup_s):
    (A, B, ncell, n_iter, sym, filename, interactive,
    savedisplay, term_outdir, minimize, test, A3D, B3D) = (\
    Structure(np.array([[3.57371   , 0.        , 0.        ],
                        [0.        , 2.52699457, 0.        ],
                        [0.        , 0.        , 2.52699457]]), name='C8 (1 1 0) 0')\
    .add_atom(0.0, 0.0, 0.0, '1', site=0)\
    .add_atom(2.6802824999999992, 1.2634972874970873, 6.661338147750939e-16, '1', site=1),
    Structure(np.array([[ 3.86697465, -6.89049062e-16,  0.        ],
                        [ 0.,          3.86697465e+00,  0.        ],
                        [ 0.,          0.,              5.46872800]]), name='Si8 (0 0 1) 0')\
    .add_atom(0.0, 0.0, 0.0, '1', site=0),
    100, 1000, 1, './p2p.in', False, False, './DC_C_POSCAR-DC_Si_POSCAR/term_000-000', False, False,
    np.array([[ 0., -0.70710678,  0.70710678],
            [ 0.,  0.70710678,  0.70710678],
            [-1.,  0.,          0.,       ]]),
    np.array([[-0.70710678,  0.70710678,  0.],
            [-0.70710678, -0.70710678,  0.],
            [ 0.,          0.,          1.]]))

    ttrans, dispStruc, vec_classes, dmin = findMatchingInterfaces(A, B, ncell, 10000,
                                sym=sym, filename=filename,
                                interactive=interactive,
                                savedisplay=savedisplay,
                                outdir=term_outdir,
                                minimize=minimize, test=test,
                                A3D=A3D, B3D=B3D)
    print('ASH')
    print(ttrans, dispStruc, vec_classes, dmin)
    return (ttrans, dispStruc, vec_classes, dmin)

@pytest.mark.skip(reason="Too inconsistant...")
def test_interface_ttrans(find_matching_interface):
    #TODO: find pattern
    #[[[ 0.924162, -0.23104 ,  0.      , -0.893427], another possiblity...
    #  [ 0.326741,  0.980222,  0.      ,  1.895085],
    #  [ 0.      ,  0.      ,  1.      ,  2.247815]]]
    ttrans, dispStruc, vec_classes, dmin = find_matching_interface    
    test_ttrans = [[[ 0.92416173, -0.23104043, 0., -0.89342744],
                    [ 0.32674051, 0.98022154,  0.,  1.89508484],
                    [ 0.,         0.,          1.,  2.24781482]]]
    np.testing.assert_allclose(ttrans, test_ttrans)

@pytest.mark.skip(reason="Too inconsistant...")
def test_interface_dispStruct(find_matching_interface):
    #TODO: need to fix cell negetives
    ttrans, dispStruc, vec_classes, dmin = find_matching_interface
    test_dispStruc = Structure(np.array([[ -7.14742   ,   3.57371   ,   0.        ],
                                        [ -2.52699457, -15.16196742,   0.        ],
                                        [   0.       ,  -0.        ,   7.99572257]]))\
    .add_atom(-3.5737104474246926, -17.6888747705878, 2.247814822377832, '0', atom='1', site=0)\
    .add_atom(-0.893427947424688, -12.6348856305878, 2.247814822377832, '1', atom='1', site=1)\
    .add_atom(-4.467137947424693, -13.898382915587801, 2.247814822377832, '2', atom='1', site=2)\
    .add_atom(-4.474246906305268e-07, -16.4253774855878, 2.247814822377832, '3', atom='1', site=3)\
    .add_atom(-6.253992947424691, -6.317399205587799, 2.247814822377832, '0', atom='1', site=4)\
    .add_atom(-5.360565447424691, -10.107891060587797, 2.247814822377832, '2', atom='1', site=5)\
    .add_atom(-2.680282947424688, -5.0539019205878, 2.247814822377832, '3', atom='1', site=6)\
    .add_atom(-1.78685544742469, -8.844393775587802, 2.247814822377832, '1', atom='1', site=7)
    assert_structs_approx_eq(dispStruc[0], test_dispStruc)

def test_interface_vec_classes(find_matching_interface):
    ttrans, dispStruc, vec_classes, dmin = find_matching_interface
    vec_classes = np.array(vec_classes[0])
    #sorted_vec_classes = vec_classes[np.lexsort(np.transpose(np.round(vec_classes, 2)))]
    for vec_class in vec_classes:
        np.testing.assert_allclose(abs(vec_class), [ 0.893501,  0.631446, 2.247841], 0.001)

@pytest.mark.skip(reason="Too inconsistant...")
def test_interface_dmin(find_matching_interface):
    #TODO: what is the reasonable tolerance on this?  
    ttrans, dispStruc, vec_classes, dmin = find_matching_interface
    test_dmin = -44.17527158660921
    assert dmin == pytest.approx(test_dmin)



