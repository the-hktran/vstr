#Generic makefile for fortran90
FC := gfortran
# serial
#LIBS := -qmkl=sequential
#FFLAGS := -r8 -O3
# parallel  
LIBS := #-mkl=sequential
FFLAGS := -fdefault-real-8 -fdefault-double-8 -O3 -fopenmp -fPIC

OBJ := \
constants.o\
model_mod.o\
ttm2f_mod.o\
ttm3f_mod.o\
nasa_mod.o\
math.o\
smear.o\
nasa.o\
potential_mod.o\
bemsa2b.o\
bemsa3b.o\
bemsa4b.o\
pes_shell.o\

OBJPOT := \
calc_h2o_n_pot.o

%.o : %.f90
	$(FC) -c $(FFLAGS) $< 

all : $(OBJ) $(OBJPOT)
	$(FC) -c -fPIC $(FFLAGS) calc_h2o_n_pot.f90
	python3 -m numpy.f2py -lgomp -c $(OBJ) -m h2o_n_pot calc_h2o_n_pot.f90

clean :
	rm *.o *.mod 
