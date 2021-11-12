from p2ptrans import analysis, findMatching
import pytest
import os, glob
import numpy as np
from test_units import read_FCC_BCC, BCC_file, FCC_file, tol, cleanup
import cProfile


def default_options():
    # (fileA, fileB, ncell, filename, interactive, savedisplay, outdir,
    # use, switch, prim, anim, vol, minimize, test, crystfile, n_steps,
    # showversion, map_ncell)
    return ('./POSCAR_A', './POSCAR_B', 300, './p2p.in', False, False, '.',
     False, False, True, False, False, False, False, './cryst.in', 60,
     False, None)

if __name__ == "__main__":
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

    matchCode = 'findMatching(BCC, FCC, ncell, fileA=BCC_file, fileB=FCC_file,'+\
                    'ccellA=ccell1, ccellB=ccell2,'+\
                    'filename=filename, interactive=interactive,'+\
                    'savedisplay=savedisplay, outdir=outdir,'+\
                    'switch=switch, prim=prim, vol=vol,'+\
                    'minimize=minimize, test=test, map_ncell=map_ncell)'

    cProfile.run(matchCode, 'stats')
    

    cleanup()
