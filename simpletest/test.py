#!/usr/bin/env python
import os
import p2ptrans as p2p
from pylada.crystal import read

bcc = read.poscar(os.path.join('ref','BCC_Fe_POSCAR'))
fcc = read.poscar(os.path.join('ref','FCC_Fe_POSCAR'))

outdir = 'results'

# n_steps, (tmat, dispStruc, vec_classes), outdir, display
p2p.produceTransition(60, *p2p.findMatching(bcc, fcc, 100, outdir=outdir), outdir, False)

# then verifiy


