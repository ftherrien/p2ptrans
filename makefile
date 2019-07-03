COMP = gfortran
COMP2 = python3 -m numpy.f2py

files = tiling.f90 transform.f90
functions = circle trans center fastoptimization
flags = --f90exec=gfortran --f90flags="-fopenmp" -lgomp 

#-I/usr/lib64/openmpi/lib/ -L/usr/lib64/openmpi/lib/ -lmpi

all: p2ptrans

p2ptrans: $(files) lap.f90
	$(COMP) -c -fPIC lap.f90	
	$(COMP2) -c $(flags) --opt=-O3 -I. lap.o -m $@ $(files) only: $(functions)

debug: $(files) lap.f90
	$(COMP) -c -fPIC lap.f90		
	$(COMP2) -c $(flags) --opt=-Og -I. lap.o -m p2ptrans $(files) only: $(functions)

clean:
	rm -f *.mod *.so
