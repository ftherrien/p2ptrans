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
   
