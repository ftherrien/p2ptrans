# p2ptrans - A Structure Matching Algorithm

<img src="https://github.com/ftherrien/p2ptrans/blob/master/WelcomeImage.gif" width="300" height="300">

p2ptrans allows you to find the best matching between two crystal structures.

## Latest Updates

Current version: 2.1.0 (07.15.2022)

**03.14.2021**: VERSION 2.0 is out! See the section below for more details.  
**10.04.2020**: Added dmin as a 4th output to the findMatching() function in p2ptrans.  
**8.28.2020**: More detailed documentation is now available for [p2ptrans](https://p2ptrans.readthedocs.io/en/latest/p2ptrans.html) and [p2pint](https://p2ptrans.readthedocs.io/en/latest/p2pint.html).  

## What's new in 2.0
Version 2.0 includes a new method to find the periodicity, which means p2ptrans and p2pint become much more reliable at finding the transformation or interface structures.

In 3D, in many cases, this means that p2ptrans can be run with a smaller "-n" and still return a relevant result. Moreover, even if the optimal matching does not lead to the correct ratio specific volumes (a necessary condition for a one-to-one mapping), p2ptrans will return a transformation path that involves the creation of vacancies in one of the structures.

In 2D, the printed information has been clarified and several bugs have been fixed (see commits for details).

## Features
p2ptrans (and p2pint) can be used directly as a command-line interface (cli) or as a python package. It can be used for two main aspects:

### 1. Phase Transformations:
p2ptrans can find the optimal mechanism of transformation between any two structures. It can provide the following information:
* The transformation cell
* The evolution of the structures during the transformation in the [POSCAR](https://www.vasp.at/wiki/index.php/Input) format (ex: 60 steps)
* The total distance traveled by all the atoms during the transformation
* The principal strains and directions
* The uniformly strained plane (Habit Plane)
* The orientation relationship (constrained and unconstrained)
* An animation of the transformation from different points of view

### 2. Interfaces
Given the interfacial planes, p2pint finds the optimal matching between two interfaces. It can provide the following information:
* The distance between the two structures (how well they match)
* The cell of correspondence between the two structures (Interface Cell)
* The amount of strain at the interface
* A POSCAR file representing the interface for each termination

## Installation
    pip install git+https://github.com/ftherrien/p2ptrans

Note: If you do not have [pylada](https://github.com/pylada/pylada-light), you will need to install the py module first:
    pip install py


### Possible Errors
1. On certain systems, the pylada installation fails with `error: ‘v’ does not name a type`. If you encounter this error retry the installation with:
```
CXXFLAGS="-std=c++11" pip install git+https://github.com/ftherrien/p2ptrans
```
 
## Documentation & Tutorials

Please visit the [documentation for p2ptrans](https://p2ptrans.readthedocs.io)

To run the *transformation* finder:
    
    p2ptrans -I POSCAR_INITIAL -F POSCAR_FINAL
    
to get help:
    
    p2ptrans --help

To run the *interface* finder:

    p2pint -T POSCAR_TOP -B POSCAR_BOTTOM

to get help:
    
    p2pint --help

## Contribution
Any contribution including [raising issues](https://github.com/ftherrien/p2ptrans/issues) is greatly appreciated.

## Citation
If you use p2ptrans in your research, please cite:

[Therrien, F., Graf, P., & Stevanović, V. (2020). Matching crystal structures atom-to-atom. *The Journal of Chemical Physics, 152(7)*, 074106.](https://aip.scitation.org/doi/full/10.1063/1.5131527) [ArXiv](https://arxiv.org/abs/1909.12965)
