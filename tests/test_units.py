from p2ptrans import read, analysis, core, format_spglib
from p2ptrans.fmodules import hungarian as lap
from conftest import assert_structs_approx_eq
from spglib import get_spacegroup
import numpy as np
import pytest


def test_read_poscars(bcc_fcc, bcc_fcc_filenames):
    '''Test that pylada reads test files correctly'''
    bcc, fcc = bcc_fcc
    bcc_file, fcc_file = bcc_fcc_filenames
    read_bcc = read.poscar(bcc_file)
    read_fcc = read.poscar(fcc_file)

    assert_structs_approx_eq(read_bcc, bcc)
    assert_structs_approx_eq(read_fcc, fcc)

def test_spacegroup(bcc_fcc):
    '''Tests that spglib reports the correct spacegroups'''
    bcc, fcc = bcc_fcc
    assert get_spacegroup(format_spglib.to_spglib(bcc), symprec=0.3, angle_tolerance=3.0) == 'Im-3m (229)'
    assert get_spacegroup(format_spglib.to_spglib(fcc), symprec=0.3, angle_tolerance=3.0) == 'Fm-3m (225)'

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

def test_p2pint_extra_params():
    assert False

@pytest.mark.skip(reason="Replaced with smaller and larger tests...")
def test_optimize(double_cleanup, bcc_fcc):
    BCC, FCC = bcc_fcc
    BCC_cell = BCC.cell * 2.87
    FCC_cell = FCC.cell * 3.57
    # optimization(A, Acell, mulA, B, Bcell, mulB, ncell, filename, outdir, max_cell_size)
    result = core.optimization(BCC, BCC_cell, 1, FCC, FCC_cell, 1, 300, './p2p.in', '.', 1000)
    Apos, Apos_map, Bpos, Bposst, n_map, natA, class_list, tmat, dmin, atoms, atom_types, foundcell, vec = result
    print(Apos, Apos_map, Bpos, Bposst, n_map, natA, class_list, tmat, dmin, atoms, atom_types, foundcell, vec)


def test_LAP(lap_test_matrix, lap_test_map):
    test_mat = lap_test_matrix
    test_map = lap_test_map
    dist, mapping = lap.munkres(test_mat, )
    assert dist == pytest.approx(94.08380269604208, 0.0001)
    np.testing.assert_allclose(mapping, test_map)
    
'''!ASHHH ! way to print given matrix
    integer :: L

    Print *, "start Matrix"
    DO i = 1, N        ! print row-by-row
      DO j = 1, N
        print *, C(i, j)
      END DO
      !WRITE(*,*) (C(i,j), j=1, L)
    END DO
    print *, "end Matrix"
    stop
    !END ASHHHH'''
