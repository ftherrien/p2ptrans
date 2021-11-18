from p2ptrans import analysis, findMatching, produceTransition
from p2ptrans.interfaces import findMatchingInterfaces
from conftest import assert_structs_approx_eq
from pylada.crystal import Structure
import numpy as np
import pytest
import os


class TestTransitions:
    @pytest.fixture(scope="class")
    def find_matching(bcc_fcc, bcc_fcc_filenames_s, double_cleanup_s, p2p_default_options_s):
        '''Runs a full matching, then cleanup the files'''
        (_, _, ncell, filename, interactive, savedisplay, outdir,
            use, switch, prim, anim, vol, minimize, test, crystfile, n_steps,
            showversion, map_ncell) = p2p_default_options_s
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
        '''Test that the tmat is similar to expected'''
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
        '''Test that the dispCells match expected'''
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
        '''Test that the eigenval, U, P, Q, and plane habit match expected'''
        test_tmat = [[-8.03921569e-01, -8.03921569e-01, 0],
                    [ 8.03921569e-01, -8.03921569e-01, 0],
                    [ 0, 0, 8.03921569e-01]]
        bcc, fcc = bcc_fcc
        BCC_file, FCC_file = bcc_fcc_filenames
        ccell1, ccell2, planehkl, diruvw = analysis.readCrystParam('./FILE_DNE')
        eigval, U, P, Q, planeHab = analysis.crystallography(np.linalg.inv(test_tmat), fcc, bcc, ccell2, ccell1, planehkl,
                                                                diruvw, fileA=BCC_file, fileB=FCC_file)
        
        assert eigval == pytest.approx([0.87957185, 0.87957185, 1.24390244], 0.0001)
        np.testing.assert_allclose(U, [[0.88275147, 0.05278822, 0.        ],
                                    [0.05278822, 0.87639223, 0.        ],
                                    [0.,         0.,         1.24390244]])
        np.testing.assert_allclose(P, [[1.,         0.06012458, 0.,       ],
                                    [0. ,        0.99819088, 0.,        ],
                                    [0.,         0.,         1.,        ]])
        np.testing.assert_allclose(Q, [[0.43187329, 1.,         0.,        ],
                                    [0.90193429, 0.,         0.,        ],
                                    [0.,         0.,         1.,        ]])
        np.testing.assert_allclose(planeHab, [[ 0., -0.],
                                            [ 0., -0.],
                                            [ 1., -1.]])

    def test_transition(double_cleanup):
        '''Test that the transition pathway produced matches expected'''
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
        assert_structs_approx_eq(path[5], res5)`




class TestInterfaces:
    @pytest.fixture(scope="class")
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
        return ttrans, dispStruc, vec_classes, dmin

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
        #TODO: needs to be sorted or something
        ttrans, dispStruc, vec_classes, dmin = find_matching_interface
        test_vec_classes = [[np.array([-0.89342756, -0.63158755, -2.24781482]),
                            np.array([ 0.89342744,  0.63190973, -2.24781482]),
                            np.array([ 0.89342744, -0.63158756, -2.24781482]),
                            np.array([-0.89342756,  0.63190973, -2.24781482])]]   
        np.testing.assert_allclose(vec_classes, test_vec_classes, 0.0001)


    def test_interface_dmin(find_matching_interface):
        #TODO: what is the reasoable tolerance on this?  
        ttrans, dispStruc, vec_classes, dmin = find_matching_interface
        test_dmin = -44.17527158660921
        assert dmin == pytest.approx(test_dmin)
    

    
