COMP = gfortran
COMP2 = python3 -m numpy.f2py

files = tiling.f90 transform.f90
functions = sphere parallelepiped circle trans center fastoptimization canonicalize
flags = --f90exec=gfortran --f90flags="-fopenmp" --opt=-O3 -lgomp 

#-I/usr/lib64/openmpi/lib/ -L/usr/lib64/openmpi/lib/ -lmpi

all: fmodules

fmodules: $(files) lap.f90
	$(COMP) -c -fPIC lap.f90	
	$(COMP2) -c $(flags) -I. lap.o -m $@ $(files) only: $(functions)

clean:
	rm -f *.mod *.so
