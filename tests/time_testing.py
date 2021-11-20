from p2ptrans import analysis, findMatching
from p2ptrans.core import optimization
import pytest
import os, glob
import numpy as np
from conftest import _bcc, _fcc, BCC_filename, FCC_filename, cleanup
import cProfile, pstats

def default_options():
    # (fileA, fileB, ncell, filename, interactive, savedisplay, outdir,
    # use, switch, prim, anim, vol, minimize, test, crystfile, n_steps,
    # showversion, map_ncell)
    return ('./POSCAR_A', './POSCAR_B', 300, './p2p.in', False, False, '.',
     False, False, True, False, False, False, False, './cryst.in', 60,
     False, None)

def profile_match_code():
    cleanup()
    (_, _, ncell, filename, interactive, savedisplay, outdir,
        use, switch, prim, anim, vol, minimize, test, crystfile, n_steps,
        showversion, map_ncell) = default_options()
    
    BCC, FCC = _bcc, _fcc

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    try:
        with open(filename, "r") as f:
            filecontent = f.readlines()
    except FileNotFoundError:
        filecontent = ""

    ccell1, ccell2, planehkl, diruvw = analysis.readCrystParam('./FILE_DNE')

    matchCode = 'findMatching(BCC, FCC, ncell, fileA=BCC_filename, fileB=FCC_filename,'+\
                    'ccellA=ccell1, ccellB=ccell2,'+\
                    'filename=filename, interactive=interactive,'+\
                    'savedisplay=savedisplay, outdir=outdir,'+\
                    'switch=switch, prim=prim, vol=vol,'+\
                    'minimize=minimize, test=test, map_ncell=map_ncell)'

    cProfile.runctx(matchCode, globals(), locals(), 'stats', sort=-1)   

    cleanup()


def profile_optimization():
    cleanup()
    optCode = '''optimization(_bcc,
    np.array([[ 1.435, -1.435,  1.435],
              [ 1.435,  1.435, -1.435],
              [-1.435,  1.435,  1.435]]),
    1, _fcc,
    np.array([[1.785, 1.785, 0.   ],
              [1.785, 0.,    1.785],
              [0.,    1.785, 1.785]]),
    1, 300, './p2p.in', '.', '1000')'''
    cProfile.runctx(optCode, globals(), locals(), 'stats')
    cleanup()
    

if __name__ == "__main__":
    profile_optimization()

    # note: tottime is time spent just in function
    # cumtime includes time spent in subfunctions
    #profile = pstats.Stats('./stats')
    #profile.strip_dirs().sort_stats('cumulative')
    #profile.print_stats(15)
