p2ptrans
========

Command line
^^^^^^^^^^^^

Usage:

.. code-block:: console

   p2ptrans [-h] [-I INITIAL_STRUCTURE] [-F FINAL_STRUCTURE] [-n N_CELLS] [-i] [-d] [-p FILENAME] [-c CRYSTFILE]
            [-o OUTDIR] [-u USEDIR] [-m] [-s] [-r] [-a] [-v] [-t] [-f N_STEPS] [--version]


Optional arguments:
  -h, --help            Show help message and exit
  -I, --initial INITIAL_STRUCTURE
                        Initial crystal structure. It must be a `POSCAR <https://www.vasp.at/wiki/index.php/Input>`_.
  -F, --final FINAL_STRUCTURE
                        Final Structure structure. It must also be a `POSCAR <https://www.vasp.at/wiki/index.php/Input>`_.
  -n, --ncell N_CELL
                        Minimum number of cells to tile *DEFAULT: 300*
  -cn, N_CELL
                        If specified, number of cells used for the final mapping after tmat has been adjusted to be commensurate with both structures. If not specified, the mapping stays the same.
  -i, --interactive     Enable interactive display
  -d, --disp            Save figures
  -p, --param FILENAME
                        Name of the parameter file (see detail in `Extra parameters`_). *DEFAULT: p2p.in* **Caution**: if a file named *p2p.in* exists in the current directory, p2ptrans will read it.
  -c, --crystal CRYSTFILE
                        Name of the parameter file for crystallography analysis *DEFAULT: cryst.in*
  -o, --outdir OUTDIR
                        Output directory *DEFAULT: current directory*
  -u, --use USEDIR      Use previously calculated data in *USEDIR*
  -m, --minimize        Force new optimization even if data is available
  -s, --switch          Map the larger cell on the smaller cell instead of the opposite
  -r, --noprim          Do not try to find the primitive cell of the input structures
  -a, --anim            Produce an animation of the transition (slow)
  -v, --vol             Make the two (stoichiometric) cells equal in volume (not tested)
  -t, --test            Test the input file and prepare the run. You can continue this run
                        with the -u USEDIR -r option
  -f N_STEPS, --frames N_STEPS
                        Number of intermediate structures to create
  --version             Show installed version of p2p



Extra Parameters
^^^^^^^^^^^^^^^^

More parameters can be specified to p2ptrans via a parameter file. This file is in the Fortran namelist format. By default any file named *p2p.in* in the current directory will be read. The name of the parameter file can be changed using the ``-p`` option. The syntax of this file is the following:

.. code-block:: fortran

   &input
   init_class = 1.0d0
   tol = 1.0d-4
   tol_std = 1.0d-7 ! If not specified: tol_std = tol*1.0d-3
   tol_class = 1.0d-3
   fracA = 0.09d0
   fracB = 0.3d0
   n_iter = 1000
   n_ana = 300
   n_conv = 5
   n_class = 30
   n_adjust = 10
   max_vol = 4.0d0
   free = .false.
   check = .false.
   savebest = "best.dat" ! If not specified: savebest = trim(outdir)//"/best2d.dat"
   usebest = .false.
   remap = .true.
   /

This is file contains all the default parameters, if an entry is not specified, it will take the value shown above.

  Init_class
               Initial separation tolerance for displacement classes. At the initial classification step, if the norm of the difference between two vectors is larger than *init_class* they will be classified in different groups.
  tol
               Convergence criterion for the gradient descent
  tol_std
               Convergence criterion for the std minimization
  tol_class
               Convergence criterion for the classification *abs(std - previous std)*
  fracA
               Fraction of the mapped structure that constitutes core atoms
  fracB
               Fraction of the mapping structure that constitutes mapping atoms
  n_iter
               Number of random starts
  n_ana
               Maximum number of iterations in the gradient descent
  n_conv
               Maximum number of remappings per minimization
  n_class
               Maximum number of classification iterations
  n_adjsut
               Maximum number of unconstrained post-processing minimization iterations
  max_vol
               Maximum volume of the random starting *tmat* when using the unrestricted minimization (*free = .true.*)
  free
               Use unrestricted minimization. Not limited to rigid rotations. The optimal result is the one for which the sum of the unstrained distance and the strained (unrestricted) distance is minimal so that overly stretched results, where skipping is likely to occur, are penalized.
  check
               Use the mapping given by the order of the input structures directly (do not map atoms).
  savebest
               Name of the file to save the optimal result to at the end of the minimization, before the post-processing steps.
  remap
               If true, allows remapping during the post-processing steps.
	       
	  






