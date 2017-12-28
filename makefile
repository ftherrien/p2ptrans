COMP = gfortran
COMP2 = f2py

mod_files = transform.f90 main.f90 
flags = --f90exec=gfortran

#-I/usr/lib64/openmpi/lib/ -L/usr/lib64/openmpi/lib/ -lmpi

all: transform main

main: $(mod_files) 
	$(COMP) $^ -o $@

transform: lap.f90 transform.f90
	$(COMP) -c -fPIC lap.f90
	$(COMP2) -c $(flags) -I. lap.o -m $@ transform.f90 only: trans center mapping 

clean:
	rm -f *.mod *.so

