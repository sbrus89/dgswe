CC = ifort
CFLAGS1 = -c -132 -heap-arrays #-C -traceback -g 
CFLAGS2 = -O3 -o 
LIB = -llapack

objects = kdtree2.o globals.o allocation.o basis.o find_element.o write_results.o connect.o evaluate.o read_input.o read_grid.o rimls.o

main: $(objects)
	$(CC) $(CFLAGS2) rimls $(objects) $(LIB)

kdtree2.o: kdtree2.F
	$(CC) $(CFLAGS1) $<	

globals.o: globals.f90 kdtree2.F
	$(CC) $(CFLAGS1) $<

allocation.o: allocation.f90 globals.f90	
	$(CC) $(CFLAGS1) $<

basis.o: basis.f90 globals.f90
	$(CC) $(CFLAGS1) $<
	
find_element.o: find_element.f90 globals.f90 basis.f90 kdtree2.F
	$(CC) $(CFLAGS1) $<

write_results.o: write_results.f90 globals.f90
	$(CC) $(CFLAGS1) $<

connect.o: connect.f90 globals.f90
	$(CC) $(CFLAGS1) $<

evaluate.o: evaluate.f90 basis.f90 globals.f90 kdtree2.F
	$(CC) $(CFLAGS1) $<		

read_input.o: read_input.f90 globals.f90
	$(CC) $(CFLAGS1) $<	

read_grid.o: read_grid.f90 globals.f90
	$(CC) $(CFLAGS1) $<

rimls.o: rimls.f90 globals.f90 evaluate.f90 allocation.f90 write_results.f90 basis.f90 kdtree2.F
	$(CC) $(CFLAGS1) $<

.PHONY : clean
clean : 
	rm rimls *.o *.mod	
