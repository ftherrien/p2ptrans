Using as a Python Module
========================

Script Example for Interfaces
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This script will create the interface POSCARs for the first termination in Si and SiC with 3 layers and no
vacuum.

.. code-block:: python

   from p2ptrans.interfaces import read, readSurface, findMatchingInterfaces, createPoscar

   Si = read.poscar("POSCAR_Si") # Create a pylada Structure object for Si    
   SiC = read.poscar("POSCAR_SiC") # Create a pylada Structure object for SiC
   
   termSi = readSurface(Si, [1,1,0], {"1":{"Si"}}) # Create a list of possible terminations for Si   
   termSiC = readSurface(SiC, [0,0,1], {"1":{"Si"}}) # Create a list of possible termination for SiC
   
   A, equivTermA = termSi[0]    
   B, equivTermB = termSiC[0]

   ttrans, dispStruc, vec_classes, dmin = findMatchingInterfaces(A, B, 130, 1000) # Run the minimization algorithm with 130 cells and 1000*10 random initializations
   
   createPoscar(A, B, equivTermA[0], equivTermB[0], ttrans[0,:,:], dispStruc[0], layers=3, vacuum=0) # Creates the interface POSCAR for the first variant of termination of A and B
   
Script Example for Transformations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This script was used to create a transforming martensite plate sandwiched between two austenite slabs for figure 3 of `Therrien and Vladan Stevanović. "Minimization of atomic displacements as a guiding principle of the martensitic phase transformation", Physical Review Letters (Accepted) <https://arxiv.org/abs/1912.11915>`_.

.. code-block:: python

   from p2ptrans.interfaces import read, readSurface, findMatchingInterfaces, createPoscar

   import p2ptrans as p2p
   import numpy as np
   from pylada.crystal import write, read, supercell, Structure 
   import numpy.linalg as la
   from copy import deepcopy

   layBCC = 6 # Number of BCC (martensite) layers in the final state 
   layFCC = 4 # Number of FCC (austenite) layers in the final state
   dirname = "outdir"
   tol = 1e-3  

   # Read the position files (POSCAR) using the read pylada function
   A = read.poscar("POSCAR_FCC") 
   B = read.poscar("POSCAR_BCC")

   # Find best matching between A and B (from A to B)
   # Use minimize=False to avoid reoptimizing if there is already a result in "opt_dir"
   tmat, dispStruc, vec_classes = p2p.findMatching(A,B, 600, outdir="opt_dir", minimize=False)

   # Find the habit plane (planeHab) using the crystallography function
   # tmat transforms B into A, so we use la.inv(tmat) to get strain directions in coordinates of A
   eigval, U, P, Q, planeHab = p2p.analysis.crystallography(la.inv(tmat), B, A)

   # Find the closest hkl in the basis of the convential cell (cubic)
   planeHab = p2p.utils.find_uvw(planeHab[:,0:1], np.eye(3)).reshape(3)
   
   # Using a function from the interfaces submodule find a unit cell of the transformation structure (dispStruc)
   # where the cell vectors a and b are in the habit plane 
   cell2D, cell3D = p2p.interfaces.find_basis(planeHab, dispStruc.cell)

   # Replace the atoms in the structure
   dispStruc = supercell(dispStruc, cell3D)

   # Produce the transition structures with 15 steps
   result = p2p.produceTransition(15, tmat, dispStruc, vec_classes, dirname, False, habit=0)
   transStruc, spgList, Tpos, color_array, atom_types = result
       

   # The FCC (austenite) structure is the initial state [0] 
   FCC = transStruc[0]

   # Create a supercell with the specified number of layers
   FCC = supercell(FCC, FCC.cell.dot(np.array([[1,0,0],[0,1,0],[0,0,layFCC//2]])))

   
   for i, S in enumerate(transStruc):
       # Create a supercell with the specified number of layer for the transforming part of the cell
       S = supercell(S, S.cell.dot(np.array([[1,0,0],[0,1,0],[0,0,layBCC]])))

       # Create cell with inplane dimensions of FCC and z vector of S
       cell = deepcopy(FCC.cell)
       cell[:,2] = S.cell[:,2]

       # Create a structure to contain all the layers 
       slab = Structure(FCC.cell)
       slab.cell[:,2] = cell[:,2] + FCC.cell[:,2]*2

       print("Slab thickness is at step %d is:"%i, np.cross(cell[:,0],cell[:,1]).dot(cell[:,2])/la.norm(np.cross(cell[:,0],cell[:,1])))

       # Add the FCC atoms on both sides
       for a in FCC:
           slab.add_atom(*a.pos, a.type)
           slab.add_atom(*(a.pos + cell[:,2] + FCC.cell[:,2]), a.type)

       # Add the S atoms in the middle
       for a in S:
           if 1 - la.inv(S.cell).dot(a.pos)[2] < tol:
               a.pos = a.pos - S.cell[:,2]
           slab.add_atom(*(cell.dot(la.inv(S.cell).dot(a.pos)) + FCC.cell[:,2]), a.type)

       slab = supercell(slab, slab.cell) # Make sure all atoms are inside the cell

       write.poscar(slab, vasp5=True, file= dirname + "/POSCAR_%03d"%i) # Write resulting structure in POSCAR
