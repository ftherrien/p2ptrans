COMP = gfortran
COMP2 = python3 -m numpy.f2py

files = source/tiling.f90 source/utils.f90 source/potential.f90 source/transform.f90
functions = sphere circle center fastoptimization intoptimization
flags = --f90exec=gfortran --f90flags="-g -fbacktrace -fopenmp" -lgomp 

#-I/usr/lib64/openmpi/lib/ -L/usr/lib64/openmpi/lib/ -lmpi

all: fmodules

fmodules: $(files) source/lap.f90
	$(COMP) -c -fPIC source/lap.f90	
	$(COMP2) -c $(flags) -I. lap.o -m p2ptrans.$@ $(files) only: $(functions)

gfortran: $(files) source/lap.f90
	$(COMP) -c -fPIC -fopenmp -Wall source/lap.f90 $(files)
clean:
	rm -f *.mod *.so
