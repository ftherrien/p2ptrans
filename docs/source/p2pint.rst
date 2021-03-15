p2pint
======

Command Line
^^^^^^^^^^^^

Usage:

.. code-block:: console

   p2pint [-h] [-B BOTTOM_STRUCTURE] [-T TOP_STRUCTURE] [-n N_CELL] [-N N_ITER] [-i] [-d]
              [-p FILENAME] [-o OUTDIR] [-u USEDIR] [-m] [-t] [-b] [-s SYM] [-l LAYERS]
              [-v VACUUM] [-w WIDTH] [--version]

Optional arguments:
  -h, --help            Show help message and exit
  -B, --bottom BOTTOM_STRUCTURE
                        Structure at the bottom (substrate), miller indices and chemistry rule. This structure will not be strained. Format: *BOTTOM_STRUCTURE [h,k,l] 1 EL11 [EL12 ...] [2 EL21 [EL22 ...] ...]*  Example: *POSCAR_A [0,0,1] 1 Si*
  -T, --top TOP_STRUCTURE
                        Structure at the top, miller indices and chemistry rule. This structure will be strained.  Format: *TOP_STRUCTURE [h,k,l] 1 EL11 [EL12 ...] [2 EL21 [EL22 ...] ...]* Ex: POSCAR_B [1,1,0] 1 C
  -n, --ncell NCELL
                        Number of cells to tile *DEFAULT: 100*
  -cn N_CELL
                        Number of cells used for the final mapping after tmat has been adjusted to be commensurate with both structures. *DEFAULT: N_CELL*
  -N, --niter N_ITER
                        Number of random initial deformation matrices. *DEFAULT: 1000*
  -i, --interactive     Enable interactive display
  -d, --disp            Save figures
  -p, --param FILENAME
                        Name of the parameter file (see detail in `Extra parameters`_). *DEFAULT: p2p.in* **Caution**: if a file named *p2p.in* exists in the current directory, p2ptrans will read it.
  -o, --outdir OUTDIR
                        Output directory *DEFAULT: current directory*
  -u, --use USEDIR         Use previously calculated data in *USEDIR*
  -m, --minimize        Force new optimization even if data is available
  -t, --test            Test the input file and prepares the run, you can continue this run
                        with the -u USEDIR option
  -b, --surface-only    Only use the first termination from the bottom and do not find primitive cells (compatible with `cutsurfaces <https://www.github.com/ftherrien/cutsurfaces>`_)
  -s, --sym SYM         Maximum rotational symmetry around the axis perpendicular to the surface planes on both side of the interface. This will reduce the number of random starts (*N_ITER*) necessary to find the global minimum.
  -l, --layers LAYERS
                        Desired number of layers in the output interface structures
  -v, --vacuum VACUUM
                        Desired thickness of a vacuum in the output interface structures
  -w, --width WIDTH
                        Maximum width of corrugated surface to optimize. This defines a tolerance to account for atoms that are not exactly in-plane but are within *WIDTH* of the plane.
  --version             Show installed version of p2p

Extra Parameters
^^^^^^^^^^^^^^^^

More parameters can be specified to p2pint via a parameter file. This file is in the Fortran namelist format. By default any file named *p2p.in* in the current directory will be read. The name of the parameter file can be changed using the ``-p`` option. The syntax of this file is the following:

.. code-block:: fortran

   &input2d
   tol = 1.0d-6
   tol_std = tol*1.0d-3
   tol_class = 1.0d-3
   init_class = 1.0d0
   fracA = 0.0d0
   fracB = 0.25d0
   n_ana = 300
   n_conv = 5
   n_class = 30
   n_adjust = 10
   findpeaks = .false.
   free = .true.
   max_vol = 0.08d0
   savebest = trim(outdir)//"/best2d.dat"
   usebest = .false.
   remap = .true.
   pot = "LJ"
   param = 2.5d0
   check = .false.
   vecrep = 10
   min_prom = 0.6d0
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
               Fraction of the mapped structure that constitutes core atoms. For interfaces matching *fracA=0* by default so that one-to-one mapping is **not** enforced. 
  fracB
               Fraction of the mapping structure that constitutes mapping atoms
  n_ana
               Maximum number of iterations in the gradient descent
  n_conv
               Maximum number of remappings per minimization
  n_class
               Maximum number of classification iterations
  n_adjsut
               Maximum number of unconstrained post-processing minimization iterations
  max_vol
               Maximum in-plane strain or relative change in area *during the minimization* ( when *free = .true.*). Note that the final post-processed result may not meet this criterion.  
  free
               Use unrestricted minimization. Not limited to rigid rotations.
  savebest
               Name of the file to save the optimal result to at the end of the minimization, before the post-processing steps.
  remap
               If true, allows remapping during the post-processing steps.
  vecrep
               For each deformation matrix (*tmat*) try *vecrep* random initial translations. The total number of random starts is ``n_iter * vecrep``.
  findpeaks
               Find multiple local minima in distance, not only the absolute minimum. p2pint will find peaks in the distance vs. angle plot to determine the local minima. It will select the peaks that have a prominence greater than *min_prom*. Post-processing steps will be applied to all the selected minima.
  min_prom
               Minimum prominence of the peaks in the distance vs. angle plot to be selected for post-processing
  pot
               Type of potential to minimize. Currently, the choices are:
	       
	       :"LJ":          Lennard-Jones potential
	       :"Euclidean":   Euclidean distance		
  param
               Equilibrium length for the Lennard-Jones potential

