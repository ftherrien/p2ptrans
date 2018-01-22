COMP = gfortran
COMP2 = python3 -m numpy.f2py

files = transform.f90 tiling.f90
functions = trans center mapping fastmapping sphere parallelepiped
flags = --f90exec=gfortran

#-I/usr/lib64/openmpi/lib/ -L/usr/lib64/openmpi/lib/ -lmpi

all: p2ptrans

p2ptrans: $(files) lap.f90
	$(COMP) -c -fPIC lap.f90
	$(COMP2) -c $(flags) -I. lap.o -m $@ $(files) only: $(functions)

clean:
	rm -f *.mod *.so
