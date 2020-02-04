# P2PTRANS - A Structure Matching Algorithm

![A cool gif](https://github.com/ftherrien/p2ptrans/blob/master/WelcomeImage.gif)

p2ptrans allows you to find the best matching between two crystal structures.

## Features
p2ptrans can be used directly as a command line interface (cli) or as a python package. It can be used for two main aspects:

### 1. Phase Transformations:
p2ptrans can find the optimal mechanism of transformation between any two structures. It can provides the following information:
* The transfromation cell
* The evolution of the structures during the transformation in the [POSCAR](https://www.vasp.at/wiki/index.php/Input) format (ex: 60 steps)
* The total distance travelled by all the atoms during the transformation
* The principal strains and directions
* The uniformly strained plane (Habit Plane)
* The orientation relationship (constrained and unconstrained)
* An animation of the transfromation from different points of vue

### 2. Interfaces
The interface analysis (2D) is currently under developpement in the [2D](https://github.com/ftherrien/p2ptrans/tree/2D) branch. **It will be available soon.**

## Installation
    pip install git+https://github.com/ftherrien/p2ptrans
    
## Documentation

To run:
    
    p2ptrans -I POSCAR_INITIAL -F POSCAR_FINAL
    
to get help:
    
    p2ptrans --help

**More exhaustive documentation will be available soon.**

## Contribution
Any conribution including [raising issues](https://github.com/ftherrien/p2ptrans/issues) is greatly appreciated.

## Citation
If you use p2ptrans in your research, please cite:

> Therrien, F, Graf, G, Stevanovic, V, (2020). Matching Crystal Structures atom-by-atom. *The Journal of Chemical Physics* (Accepted), [ArXiv](https://arxiv.org/abs/1909.12965)
