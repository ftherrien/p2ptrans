COMP = gfortran
COMP2 = python3 -m numpy.f2py

files = source/lap.f90 source/utils.f90 source/tiling.f90 source/potential.f90 source/transform.f90
functions = munkres free_trans rot_mat center eye norm split det sort sphere circle distance derivative closest fixed_tmat intoptimization fastoptimization
flags = --f90exec=gfortran --f90flags="-g -Og -fbacktrace -fopenmp" -lgomp 

#-I/usr/lib64/openmpi/lib/ -L/usr/lib64/openmpi/lib/ -lmpi

all: fmodules

fmodules: $(files)
	$(COMP2) -c $(flags) -m p2ptrans.$@ $(files) only: $(functions)

gfortran: $(files) source/lap.f90
	$(COMP) -c -fPIC -fopenmp -Wall -g -Og source/lap.f90 $(files)
clean:
	rm -f *.mod *.so
