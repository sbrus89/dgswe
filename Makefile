CC = ifort
CFLAGS1 = -c -O3 -132 #-vec-report3 #-Dadcirc
# CFLAGS1 = -c -C -g -132 -traceback
CFLAGS2 = -O3 -o 
LIB = -llapack
#LIB = -Wl,--start-group   /afs/crc.nd.edu/x86_64_linux/intel/12.1/mkl/lib/intel64/libmkl_intel_lp64.a /afs/crc.nd.edu/x86_64_linux/intel/12.1/mkl/lib/intel64/libmkl_sequential.a  /afs/crc.nd.edu/x86_64_linux/intel/12.1/mkl/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm

objects = globals.o allocation.o basis.o evaluate.o read_input.o read_grid.o read_solution.o kdtree2.o find_nesting.o area_qpts.o error.o


main: $(objects)
	$(CC) $(CFLAGS2) error $(objects) $(LIB)

globals.o: globals.f90
	$(CC) $(CFLAGS1) $<
	
allocation.o: allocation.f90 globals.f90
	$(CC) $(CFLAGS1) $<
	
basis.o: basis.f90 globals.f90
	$(CC) $(CFLAGS1) $<
	
evaluate.o: evaluate.F90 basis.f90 globals.f90
	$(CC) $(CFLAGS1) $<

read_input.o: read_input.f90 globals.f90
	$(CC) $(CFLAGS1) $<

read_grid.o: read_grid.f90 globals.f90
	$(CC) $(CFLAGS1) $<

read_solution.o: read_solution.F90 globals.f90
	$(CC) $(CFLAGS1) $<

kdtree2.o: kdtree2.F
	$(CC) $(CFLAGS1) $<

find_nesting.o: find_nesting.f90 globals.f90 kdtree2.F
	$(CC) $(CFLAGS1) $<

area_qpts.o: area_qpts.f90
	$(CC) $(CFLAGS1) $<

error.o: error.F90 globals.f90 evaluate.F90 basis.f90 kdtree2.F
	$(CC) $(CFLAGS1) $<

.PHONY : clean
clean : 
	rm error *.o *.mod	
