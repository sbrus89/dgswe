CC = ifort
CFLAGS1 = -c -O3 -C -g -132 -traceback
CFLAGS2 = -O3 -o 
LIB = -llapack

# objects = globals.o basis.o evaluate.o read_input.o read_grid.o read_stations.o kdtree2.o error.o
objects = globals.o allocation.o read_input.o read_grid.o kdtree2.o error.o


main: $(objects)
	$(CC) $(CFLAGS2) error $(objects) $(LIB)

globals.o: globals.f90
	$(CC) $(CFLAGS1) $<
	
allocation.o: allocation.f90 globals.f90
	$(CC) $(CFLAGS1) $<
	
basis.o: basis.f90 globals.f90
	$(CC) $(CFLAGS1) $<
	
evaluate.o: evaluate.f90 basis.f90 globals.f90
	$(CC) $(CFLAGS1) $<

read_input.o: read_input.f90 globals.f90
	$(CC) $(CFLAGS1) $<

read_grid.o: read_grid.f90 globals.f90
	$(CC) $(CFLAGS1) $<

read_stations.o: read_stations.f90 globals.f90
	$(CC) $(CFLAGS1) $<

kdtree2.o: kdtree2.F
	$(CC) $(CFLAGS1) $<

error.o: error.f90 globals.f90 evaluate.f90 basis.f90 kdtree2.F
	$(CC) $(CFLAGS1) $<

.PHONY : clean
clean : 
	rm error *.o *.mod	