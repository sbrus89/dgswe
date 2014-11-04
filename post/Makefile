FC = ifort

dgpost : dgpost.F90
	$(FC) $< -o $@
	cp dgpost ../work/
	
.PHONY : clean
	rm *.o dgpost