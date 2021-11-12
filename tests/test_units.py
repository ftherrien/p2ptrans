from p2ptrans import read, analysis


BCC_file = './BCC_POSCAR'
FCC_file = './FCC_POSCAR'
tol = 0.0001


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

def test_read_cryst():
    '''Tests that the cryst param file is read correctly'''
    ccell1, ccell2, planehkl, diruvw = analysis.readCrystParam('./cryst_test_file')
    assert ccell1 == [[1., 2., 3.],[4., 5., 6.],[7., 8., 9.]]
    assert ccell2 == [[9., 8., 7.],[6., 5., 4.],[3., 2., 1.]]
    assert planehkl == [1, 2, 3]
    assert diruvw == [3, 2, 1]
