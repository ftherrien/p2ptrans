from pylada.crystal import Structure, primitive, gruber, read, write, supercell, space_group
import numpy as np
import numpy.linalg as la
import pickle

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
# ++++++++++++++++++++++++++++++++

# Steps
PoscarDirName = "/TransPOSCARS"
n_steps = 60
viewDirs = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
size = 3
habit = False
n_frames = 5
