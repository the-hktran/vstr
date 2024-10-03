#Generic makefile for fortran90

LIBDIR := /insomnia001/home/hkt2112/Code/vstr/examples/dms/h2o_n/lib

FC := gfortran
LIBS :=  -L$(LIBDIR) -ldms2b  
FFLAGS := -fdefault-real-8 -fdefault-double-8 -fopenmp -O3 -fPIC

vpath %.f90 src

OBJ := \
chg.o ltp2011_dms.o\
dms_shell.o\

OBJPOT := \
calc_h2o_n_dip.o
#getdms.o\

# Targets:
.PHONY : clean 

%.o : %.f90
	$(FC) -c $(FFLAGS) $< 
%.o : %.f
	$(FC) -c $(FFLAGS) $< 

#all : getdms.x 

#getdms.x :  $(OBJ) 
#	$(FC) -o $@  $^ $(FFLAGS) $(LIBS)

all : $(OBJ)
	$(FC) -c $(FFLAGS) $(LIBS) calc_h2o_n_dip.f90
	python3 -m numpy.f2py -lgomp -c $(OBJ) $(LIBDIR)/libdms2b.a -m h2o_n_dip calc_h2o_n_dip.f90
	cp *.so ./..

clean :
	rm *.o *.mod *.x *.so ../*.so
