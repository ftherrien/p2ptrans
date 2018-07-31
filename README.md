# P2PTRANS
## A Lattice Independant Structure Mapping Algorithm

p2ptrans returns (1) the transformation lattice between the 2 structures (from the smallest volume to the largest) and (2) the transformation matrix applied to the smallest lattice  

Installation:

    git clone https://github.com/ftherrien/p2ptrans.git
    cd p2ptrans
    make

To execute:

   python p2ptrans.py -I POSCAR_A -F POSCAR_B

Other options:

optional arguments:
`-`h, --help                  Show a help message and exit
`-`I A, --initial A           Initial Structure
`-`F B, --final B             Final Structure
`-`n NCELL, --ncell NCELL     Number of cells to tile
`-`a FRACA, --fracA FRACA     Fraction of the biggest structure to force to use in mapping (fracA < fracB)
`-`b FRACB, --fracB FRACB     Fraction of the smallest structure to force to use in mapping
`-`r NITER, --niter NITER     Number of (r)andom starts
`-`g NANA, --nana NANA        Number of iteration in (g)radient descent
`-`m REMAP, --remap REMAP     Number of re-(m)apping
`-`s ADJUST, --adjust ADJUST  Number of (s)cale adjusting steps
`-`d, --display               Unable interactive display