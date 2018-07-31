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


`-h`, `--help`    Show a help message and exit

`-I`, `--initial` Initial Structure

`-F`, `--final`   Final Structure

`-n`, `--ncell`   Number of cells to tile

`-a`, `--fracA`   Fraction of the biggest structure to force to use in mapping (fracA < fracB)

`-b`, `--fracB`   Fraction of the smallest structure to force to use in mapping

`-r`, `--niter`   Number of (r)andom starts

`-g`, `--nana`    Number of iteration in (g)radient descent

`-m`, `--remap`   Number of re-(m)apping

`-s`, `--adjust`  Number of (s)cale adjusting steps

`-d`, `--display` Unable interactive display
