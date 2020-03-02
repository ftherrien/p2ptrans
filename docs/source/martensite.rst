Phase Transformations
=====================

This tutorial will show how to obtain the results presented in this paper: `F. Therrien and V. Stevanovic, arXiv:1912.11915 (2020)
<https://arxiv.org/abs/1912.11915>`_.

This assumes that you have installed p2ptrans (see :ref:`install`).

Initial and Final Structures
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this tutorial, we will find the transformation mechanism between the face-centered cubic and body-centered cubic phases of iron.

First let's create a directory and move our working directory to it:

.. code-block:: console

   mkdir tutorial
   cd tutorial

The initial and final structure must be given in the `POSCAR format
<https://www.vasp.at/wiki/index.php/Input>`_.

Let's create the input file for FCC and let's name it POSCAR_FCC:

.. code-block:: console

   FCC structure
     3.57
     0.5  0.5  0.0
     0.5  0.0  0.5
     0.0  0.5  0.5
     Fe
     1
   Direct
     0.0  0.0  0.0

Similarly for BCC, let's create the file POSCAR_BCC:

.. code-block:: console

   BCC structure
     2.87
     0.5  0.5 -0.5
     -0.5  0.5  0.5
     0.5 -0.5  0.5
     Fe
     1
   Direct
     0.0  0.0  0.0

Running the algorithm
^^^^^^^^^^^^^^^^^^^^^

Now that we have the initial and final structures we can simply run:

.. code-block:: console

   p2ptrans -I POSCAR_FCC -F POSCAR_BCC

If you would like to store the output files in a subdirectory (e.g. `outputdir`), run:

.. code-block:: console

   p2ptrans -I POSCAR_FCC -F POSCAR_BCC -o outputdir

The -I flag indicates the initial structure and the -F indicates the final structure. p2ptrans will return the transformation in that specific order, in this case, from FCC to BCC. This will automatically run on all the threads on your computer or node. On a relatively recent computer, this example should take less than 5 min.

Analyzing the output
^^^^^^^^^^^^^^^^^^^^

Let's analyze the standard output of p2ptrans:

.. code-block:: console

   The optimization will be performed in the reversed direction.

This is because the "mapped structure" (A) is always the one with the largest specific volume. In this case, BCC has the largest specific volume. p2ptrans will still give you the result in the order initially specified.

.. code-block:: console

   Check progress in ./progress.txt

*progress.txt* contains a list of the initial random starts that have been started and completed. 

.. code-block:: console

   Found cell!

The program found a periodic cell in the mapping. This is usually a good sign.

.. code-block:: console

   Number of classes: 1

There is only one type of connection. This means that the transformation is a distortion without any displacement. This is not consistent with the result in the paper. We will explain why below.

.. code-block:: console

   Volume stretching factor (det(T)): 1.0391327642627355
   Cell volume ratio (initial cell volume)/(final cell volume): -1.0391327619090702

The determinant of the transformation matrix is equal to the ratio in specific volumes. This should always be the case.

.. code-block:: console

   Size of the transformation cell (TC): 1

There is only one atom in the transformation cell. This is consistent with the fact that the transformation is fully distortive.

.. code-block:: console

   TC in FCC structure (../p2ptrans/examples/BCC2FCC/POSCAR_A) coordinates:
   --------Matrix--------|-----Closest uvw------
       v1    v2    v3    |    d1    d2    d3    
    1.785 -1.785 -0.000  |     1    -1     0
    0.000 -0.000 -1.785  |     0     0     1
    1.785  1.785 -1.785  |     1     1     1

   TC in BCC structure (../p2ptrans/examples/BCC2FCC/POSCAR_B) coordinates:
   --------Matrix--------|-----Closest uvw------
       v1    v2    v3    |    d1    d2    d3    
   -1.435  1.435 -1.435  |    -1     1     1
   -1.435  1.435  1.435  |    -1     1    -1
    1.435  1.435 -1.435  |     1     1     1

The coordinates of the transformation cells or those of the primitive cell since the transformation is fully distortive.

.. code-block:: console

   ----------CRYSTALLOGRAPHY----------

   Strain Directions in FCC structure (../p2ptrans/examples/BCC2FCC/POSCAR_A) coordinates:
       d1    d2    d3    
   --------Matrix--------|-----Closest uvw------
       v1    v2    v3    |    d1    d2    d3    
    0.000 -1.000  0.011  |     0     1     0
    0.000  0.011  1.000  |     0     0     1
    1.000  0.000 -0.000  |     1     0     0

   Strain Directions in BCC structure (../p2ptrans/examples/BCC2FCC/POSCAR_B) coordinates:
       d1    d2    d3    
   --------Matrix--------|-----Closest uvw------
       v1    v2    v3    |    d1    d2    d3    
   -0.000 -0.715 -0.699  |     0     1    -1
   -0.000 -0.699  0.715  |     0     1     1
    1.000 -0.000  0.000  |     1     0     0

   Strains + 1 (eigenvalues)
       e1    e2    e3    
    0.804  1.137  1.137

Those are the Bain strains and directions! *p2ptrans found the Bain correspondence* By default, p2ptrans will create two sets of atoms (spheres) of a size of *300* primitive cells, we call this parameter ``ncell``. In this case, ``ncell`` is too small to retrieve the slipping process presented in the paper. p2ptrans is not wrong, for a system of that size, the Bain path is actually the path of minimal distance; to find the actual path of minimal distance--for an infinite system, we have to make ``ncell`` as large as computationally possible.

Notice that p2ptrans created a folder named TransPOSCAR, this folder contains 60 POSCAR files that describe the evolution of the structure during the transition. If you wish to change the number of frames, specify it using the -f option. Beware, if you reduce the number of frames p2ptrans will not erase the already existing extra-frames.

Visualizing the result
^^^^^^^^^^^^^^^^^^^^^^

Before we increase ``ncell`` let's take a look at the result.

When running p2ptrans, the result is saved in different files in the output directory. p2ptrans can be rerun without having to reoptimize the result. To run p2ptrans in interactive mode (-i) and use the previous result (-u) simply run:

.. code-block:: console

   p2ptrans -i -u .

The period indicates that the output is in the current directory (.), if you specified a different directory with the -o option you must provide the path to that directory. To save the images instead of displaying them:

.. code-block:: console

   p2ptrans -d -u .

Those two options can be used simultaneously and they can be used without the -u option.

Running the algorithm on larger systems
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Let's now increase ``ncell`` to a larger number in order to obtain the result presented in the paper.

.. tip:: I like to make sure all the parameters are ok before I truly run the code. For that you can use the ``--test`` option.

	  .. code-block:: console

	     p2ptrans -I POSCAR_FCC -F POSCAR_BCC -o newoutdir -n 600 --test

	  That will tell you how many atoms will be in each sphere which will give you an idea of how big the calculations will be--this is not always trivial when inputting two non-primitive structures of different sizes. It will also create the output directory and save the parameters of the run.

We are now ready to run the calculation:

.. code-block:: console

   p2ptrans -I POSCAR_FCC -F POSCAR_BCC -o newoutdir -n 600

.. note:: If you do not want to re-enter the same parameters you can also do: 

	  .. code-block:: console

	     p2ptrans -u newoutdir -m

	  The -m option used in concert with the -u option will use (-u) the parameters found in ``newoutdir`` and run the distance minimization (-m) on them. This will yield exactly the same results as the previous command.

The calculation should take a couple of hours on a modern computer. If you are on a cluster, you can simply put the previous line in a submission script. p2ptrans is parallelized with OpenMP; it will automatically use all the cores in one node but cannot use multiple nodes.

.. tip:: I like to monitor the progress of the calculation using

	  .. code-block:: console

	     grep "Opt dist" progress.txt | wc -l

	  This will tell you how many initial random steps have completed, by default p2ptrans will do 1000 initial random steps.

**At the end of this calculation you should obtain the result presented in the article.**

Crystallography
^^^^^^^^^^^^^^^

Coming soon!


Fine-tuning the optimization
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Coming soon!

   
