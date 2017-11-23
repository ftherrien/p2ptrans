COMP = gfortran
COMP2 = f2py

mod_files = transform.f90 main.f90 
flags = --f90exec=mpif90

#-I/usr/lib64/openmpi/lib/ -L/usr/lib64/openmpi/lib/ -lmpi

all: main

main: $(mod_files) 
	$(COMP) $^ -o $@

f90_pmpaths: transform.f90
	$(COMP2) -c $(flags) $^ -m $@ only: final_fix_gruber ucell_surface vec2alpha stats_to_value

clean:
	rm -f *.mod *.so

