Interfaces
==========

This tutorial will show how to obtain the results presented in this paper: `Therrien, Félix, Peter Graf, and Vladan Stevanović. "Matching crystal structures atom-to-atom." The Journal of Chemical Physics 152.7 (2020): 074106.
<https://aip.scitation.org/doi/abs/10.1063/1.5131527>`_.

This assumes that you have installed p2ptrans (see :ref:`install`).

Top and Bottom Structures
^^^^^^^^^^^^^^^^^^^^^^^^^

In this tutorial, we will find the orientation between Diamond Si (110) and 6H-SiC (001).

First let's create a directory and move our working directory to it:

.. code-block:: console

   mkdir tutorial2
   cd tutorial2

The initial and final structure must be given in the `POSCAR format
<https://www.vasp.at/wiki/index.php/Input>`_.

Let's create the input file for Si and let's name it POSCAR_Si:

.. code-block:: console

   Si
      5.43
    1.0    0.0     0.0
    0.0    1.0     0.0
    0.0    0.0     1.0
    Si
    8
   Direct
    -0.125 -0.125 -0.125
    -0.125  0.375  0.375
     0.375 -0.125  0.375
     0.375  0.375 -0.125
     0.125  0.125  0.125
     0.125  0.625  0.625
     0.625  0.125  0.625
     0.625  0.625  0.125

Similarly for SiC, let's create the file POSCAR_SiC:

.. code-block:: console

   SiC
   1.000000
      1.54064500000000  -2.66847541642695   0.00000000000000
      1.54064500000000   2.66847541642695   0.00000000000000
      0.00000000000000   0.00000000000000  15.11975999999338
   C Si
   6 6
   Direct
      0.00000000000000   0.00000000000000   0.12540000000000    
      0.00000000000000   0.00000000000000   0.62540000000000    
      0.33333333333333   0.66666666666667   0.29215000000000    
      0.66666666666667   0.33333333333333   0.79215000000000    
      0.33333333333333   0.66666666666667  -0.04150000000000    
      0.66666666666667   0.33333333333333   0.45850000000000    
      0.00000000000000   0.00000000000000   0.00000000000000    
      0.00000000000000   0.00000000000000   0.50000000000000    
      0.33333333333333   0.66666666666667   0.16675000000000    
      0.66666666666667   0.33333333333333   0.66675000000000    
      0.33333333333333   0.66666666666667   0.83350000000000    
      0.66666666666667   0.33333333333333   1.33350000000000 

Running the algorithm
^^^^^^^^^^^^^^^^^^^^^

Now that we have the bottom (substrate) and top structures, we need to specify the miller indices of the interfacial plane in both structures along with the bonding rule. Let's run the following command:

.. code-block:: console

   p2pint -T POSCAR_Si [1,1,0] 1 Si -B POSCAR_SiC [0,0,1] 1 Si

Here we defined Si as the top structure (`-T`) and its corresponding h,k,l as (110). We defined SiC as the bottom structure (substrate) (`-B`) with the (001) interfacial plane. The rest of the command determines the bonding rule we set it as: `1 Si` which means "the first group is composed of Si". This will associate, and therefore, bond any atom labeled Si. To bond carbon with silicon instead, the command would have been:

.. code-block:: console

   p2pint -T POSCAR_Si [1,1,0] 1 Si -B POSCAR_SiC [0,0,1] 1 C

.. note:: If a termination had two different types of atoms (which is not the case here) then the rules could be: `-B POSCAR_AB [h,k,l] 1 A 2 B` and `-T POSCAR_CD [h,k,l] 1 C 2 D` which would bond A with C and B with D at the interface between AB and CD. Similarly, we could bond A and B to C with `-B POSCAR_AB [h,k,l] 1 A B` and `-T POSCAR_CD [h,k,l] 1 C`.    
   
If you would like to store the output files in a subdirectory (e.g. `outputdir`) just add `-o outputdir`:

.. code-block:: console

   p2pint -T POSCAR_Si [1,1,0] 1 Si -B POSCAR_SiC [0,0,1] 1 Si -o outputdir

This should take about 10 min to run on a laptop. p2pint will automatically use all threads on your computer so the computation time will depend on the number of cores on your computer.

.. note:: If you do not want p2pint to use all the available threads on the computer, limit the number of threads woth:
	  
   .. code-block:: console

   OMP_NUM_THREADS=1 p2pint -T POSCAR_Si [1,1,0] 1 Si -B POSCAR_SiC [0,0,1] 1 Si

Analyzing the output
^^^^^^^^^^^^^^^^^^^^

Let's analyze the standard output of p2pint:

.. code-block:: console

   Total number of atoms in each disk: 200

This is the number of atoms that will be mapped together, i.e it is the size of the cost matrix, this number has a very strong influence on the computational cost of running p2pint.

.. code-block:: console

   Check progress in ./POSCAR_SiC-POSCAR_Si/term_000-000/progress.txt

*progress.txt* contains a list of the random initializations minimizations that have been started and completed. 

.. code-block:: console

   Looking for periodic cell for peak 0...
   Found cell!

Contrary to p2ptrans, p2pint can look for the best results instead of only looking the absolute minimal distance. Each potential result represents a peak in the distance vs. angle plot. By default, however, this functionality is turned off and p2pint will only give one peak, corresponding to the optimal result. `Found cell!` indicates that p2pint has found the cell of correspondence (Interface Cell) between the two structures.

.. code-block:: console

   -------OPTIMIZATION RESULTS FOR PEAK 0--------
   
   Number of classes: 8
   Number of mapped atoms: 64
   Total distance between structures: -59.99274481751023
   Optimal angle between structures: 270.0201544091763
   Volume stretching factor (det(T)): 0.9858178881599341
   Cell volume ratio (initial cell volume)/(final cell volume): 0.7887531650921741

This block summarizes the result of the optimization. The number of classes is the number of types of connections i.e. the number of different "bonds" that were found. The total distance between the structures is the measure of how compatible they are with this choice of potential (lennard-Jones by default). The volume stretching factor indicates how much strain there is in the first layer of the top structure. The cell volume ratio indicates the ratio of specific areas of the two interfacial planes. Note that since this is a semi-coherent interface the specific area of the two planes is very different, i.e the lattices *do not match* in the conventional sense of lattice matching.

.. code-block:: console

   -----------PERIODIC CELL-----------
   
   Number of bonds in Interface Cell (IC): 8
   Number of SiC (0 0 1) 1 cells in IC: 9.998746698318255
   Number of Si (1 1 0) 0 cells in IC: 3.999999999999999

This block gives details about the Interface Cell. The number of SiC cells is not integer because of the level of precision of the classification algorithm (1e-3 by default).

.. code-block:: console

   Creating POSCARS for peak 0, bottom term. 0, top term 0
   Creating POSCARS for peak 0, bottom term. 1, top term 0
   Creating POSCARS for peak 0, bottom term. 2, top term 0
   Creating POSCARS for peak 0, bottom term. 3, top term 0
   Creating POSCARS for peak 0, bottom term. 4, top term 0
   Creating POSCARS for peak 0, bottom term. 5, top term 0
   
Once the interface cell is found, p2pint will create interface structures for each combination of possible terminations. In this case Si (110) has 1 possible termination with 4 variants that are all equivalent under translation, and SiC (001) also has 1 termination with 6 variants (i.e. the terminating plane is the same, but the rest of the structure is different).

For each termination three POSCARs are created: (1,2) Representation of Si and SiC with a common cell in the plane specified at the beginning, (3) the interface between Si and SiC. For example, if you have a POSCAR viewing software like VESTA you can run:

.. code-block:: console

   VESTA POSCAR_SiC-POSCAR_Si/term_000-000/peak_000/var_000-000/POSCAR_interface

You can adjust the number of layers of materials on each side of the interface with the `-l` option and you can adjust the amount of vacuum with the `-v` option.

At this point your output directory should have the following structure:

.. code-block:: console

   outputdir
   ├── out.txt
   ├── param.dat
   └── POSCAR_SiC-POSCAR_Si
       └── term_000-000
           ├── best2d.dat
           ├── intoptimization.dat
           ├── peak_000
           │   ├── var_000-000
           │   │   ├── POSCAR_Bottom
           │   │   ├── POSCAR_bottom
           │   │   ├── POSCAR_interface
           │   │   ├── POSCAR_Top
           │   │   └── POSCAR_top
           │   ├── var_001-000
           │   │   ├── POSCAR_Bottom
           │   │   ├── POSCAR_bottom
           │   │   ├── POSCAR_interface
           │   │   ├── POSCAR_Top
           │   │   └── POSCAR_top
           │   ├── var_002-000
           │   │   ├── POSCAR_Bottom
           │   │   ├── POSCAR_bottom
           │   │   ├── POSCAR_interface
           │   │   ├── POSCAR_Top
           │   │   └── POSCAR_top
           │   ├── var_003-000
           │   │   ├── POSCAR_Bottom
           │   │   ├── POSCAR_bottom
           │   │   ├── POSCAR_interface
           │   │   ├── POSCAR_Top
           │   │   └── POSCAR_top
           │   ├── var_004-000
           │   │   ├── POSCAR_Bottom
           │   │   ├── POSCAR_bottom
           │   │   ├── POSCAR_interface
           │   │   ├── POSCAR_Top
           │   │   └── POSCAR_top
           │   └── var_005-000
           │       ├── POSCAR_Bottom
           │       ├── POSCAR_bottom
           │       ├── POSCAR_interface
           │       ├── POSCAR_Top
           │       └── POSCAR_top
           └── progress.txt


   
Visualizing the result
^^^^^^^^^^^^^^^^^^^^^^

When running p2pint, the result is saved in different files in the output directory. p2ptrans can be rerun without having to reoptimize the result. To run p2ptrans in interactive mode (`-i`) and use the previous result (`-u`) simply run:

.. code-block:: console

   p2pint -i -u .

The period indicates that the output is in the current directory (.), if you specified a different directory with the `-o` option you must provide the path to that directory. To save the images instead of displaying them:

.. code-block:: console

   p2ptrans -d -u .

Those two options can be used simultaneously and they can be used without the -u option.

Running the algorithm on larger systems
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Let's now increase the size of the disks (number of atoms used during the minimization) in order to obtain the result presented in the paper.

.. tip:: I like to make sure all the parameters are ok before I truly run the code. For that you can use the ``--test`` option.

	  .. code-block:: console

	     p2pint -T POSCAR_Si [1,1,0] 1 Si -B POSCAR_SiC [0,0,1] 1 Si -o outputdir2 --test

	  That will tell you how many atoms will be in each disk which will give you an idea of how big the calculations will be--this is not always trivial when inputting two non-primitive structures of different sizes. It will also create the output directory and save the parameters of the run.

We are now ready to run the calculation:

.. code-block:: console

   p2pint -T POSCAR_Si [1,1,0] 1 Si -B POSCAR_SiC [0,0,1] 1 Si -n 130

.. note:: If you do not want to re-enter the same parameters you can also do: 

	  .. code-block:: console

	     p2pint -u newoutdir -m

	  The -m option used in concert with the -u option will use (`-u`) the parameters found in ``newoutdir`` and run the distance minimization (`-m`) on them. This will yield exactly the same results as the previous command.

The calculation should less than hour on a modern computer (9 min on 4-core Intel Core i7). If you are on a cluster, you can simply put the previous line in a submission script. p2ptrans is parallelized with OpenMP; it will automatically use all the cores in one node but cannot use multiple nodes.

.. tip:: I like to monitor the progress of the calculation using

	  .. code-block:: console

	     grep "Opt dist" progress.txt | wc -l

	  This will tell you how many initial random steps have completed, by default p2pint will do 10000 initial random steps.

**At the end of this calculation you should obtain the result presented in the article.**
