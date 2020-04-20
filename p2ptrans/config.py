from pylada.crystal import Structure, primitive, gruber, read, write, supercell, space_group
import numpy as np
import numpy.linalg as la
import pickle

supercell_original = supercell

def supercell(*args):
    if len(getattr(args[0], 'name', '')) > 0:
        tmpname = args[0].name
    result = supercell_original(*args)
    if len(getattr(result, 'name', '')) > 0:
        result.name = tmpname
    return result

# Default color styles
colorlist=['#929591', 'r', 'k','b','#06470c','#ceb301', '#9e0168', '#26f7fd', '#f97306', '#c20078']
reccolor=['blue','green','red']

# Default extra parameters +++++++
tol = 1e-5
tol_vol = 2*1e-3
tol_uvw = 1e-6
pca = False
nrep = 1
gif = False
lay = 1
vac = 10
# ++++++++++++++++++++++++++++++++

# Steps
PoscarDirName = "/TransPOSCARS"
viewDirs = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
size = 3
habit = False
n_frames = 5

