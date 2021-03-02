FC = gfortran
#FCFLAGS = -fopenmp
PROJDIR := $(realpath $(CURDIR))
SOURCEDIR := $(PROJDIR)/src

#Objects
read_input.o: $(SOURCEDIR)/read_input.f90
	$(FC) $(FCFLAGS) -c $(SOURCEDIR)/read_input.f90

main.o: $(SOURCEDIR)/main.f90
	$(FC) $(FCFLAGS) -c $(SOURCEDIR)/main.f90

all: read_input.o main.o
	$(FC) $(FCFLAGS) read_input.o main.o -o convert_cp2k_to_deepmd

clean:
	rm -rf *.mod *.o

realclean:
	rm -rf *.mod *.o convert_cp2k_to_deepmd