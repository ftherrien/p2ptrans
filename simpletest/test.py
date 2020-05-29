#!/usr/bin/env python
from os import path, makedirs
import p2ptrans as p2p
from pylada.crystal import read
from test_utils import *
import pickle

bcc = read.poscar(path.join('ref','BCC_Fe_POSCAR'))
fcc = read.poscar(path.join('ref','FCC_Fe_POSCAR'))

outdir = 'results'
if not path.exists(outdir):
   makedirs(outdir) 

# n_steps, (tmat, dispStruc, vec_classes), outdir, display
tmat, dispStruc, vec_classes = p2p.findMatching(bcc, fcc, 100, outdir=outdir)
p2p.produceTransition(60, tmat, dispStruc, vec_classes, outdir, False)

# check that tmats are close
pickle.dump( tmat, open(path.join(outdir,'tmat.pickle'), 'wb') )
ref_tmat = pickle.load(open(path.join('ref','tmat.pickle'), 'rb'))
if not compare_tmats(ref_tmat, tmat, tol=1e-4):
    print('FAILED TMATS!!!!')
    exit()


# then verifiy that produced poscars are close to reference
ref_structs = structs_in_dir(path.join('ref', 'TransPOSCARS'))
produced_structs = structs_in_dir(path.join('results', 'TransPOSCARS'))

for ref_struct, res_struct in zip(ref_structs, produced_structs):
    if not compare_structs(ref_struct, res_struct, tol=1e-4):
        print('FAILED!!!!')
        print()
        print(ref_struct, 'DNE', res_struct)
        exit()

print('Succeded!')



