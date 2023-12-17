main :
	gfortran -c -fPIC CH4_PES.f90
	gfortran -c -fPIC calc_ch4_pot.f
	python3 -m numpy.f2py -c -I. CH4_PES.o -m ch4_pot calc_ch4_pot.f

clean :
	rm *.o *.so
