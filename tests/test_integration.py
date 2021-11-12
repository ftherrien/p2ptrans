from p2ptrans import *

def default_options():
    # (fileA, fileB, ncell, filename, interactive, savedisplay, outdir,
    # use, switch, prim, anim, vol, minimize, test, crystfile, n_steps,
    # showversion, map_ncell)
    return ('./POSCAR_A', './POSCAR_B', 300, './p2p.in', False, False, '.',
     False, False, True, False, False, False, False, './cryst.in', 60,
     False, None)

def read_FCC_BCC():
    BCC = read.poscar('./BCC_POSCAR')
    FCC = read.poscar('./FCC_POSCAR')
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
