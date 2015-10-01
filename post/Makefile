FC = ifort

.PHONY : dgpost clean

dgpost : dgpost.F90
	$(FC) $< -o $@
	cp dgpost ../work/
	
clean :
	rm *.o dgpost
