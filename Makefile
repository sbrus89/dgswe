CC = ifort
CFLAGS1 = -c -O2 -C -g -132 -traceback
CFLAGS2 = -O2 -o 

objects = globals.o read_input.o read_grid.o read_stations.o kdtree2.o stations.o

main: $(objects)
	$(CC) $(CFLAGS2) stations $(objects)

globals.o: globals.f90
	$(CC) $(CFLAGS1) $<

read_input.o: read_input.f90 globals.f90
	$(CC) $(CFLAGS1) $<	

read_grid.o: read_grid.f90 globals.f90
	$(CC) $(CFLAGS1) $<

read_stations.o: read_stations.f90 globals.f90
	$(CC) $(CFLAGS1) $<

kdtree2.o: kdtree2.F
	$(CC) $(CFLAGS1) $<

stations.o: stations.f90 globals.f90 kdtree2.F
	$(CC) $(CFLAGS1) $<

.PHONY : clean
clean : 
	rm stations *.o *.mod	