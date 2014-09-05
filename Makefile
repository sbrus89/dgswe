CC = ifort
CFLAGS1 = -c -O2 -C -g -traceback
CFLAGS2 = -O2
LIB = -llapack

all: main

main: globals read_grid read_input connect write_grid spline 
	$(CC) $(CFLAGS2) globals.o read_grid.o read_input.o connect.o write_grid.o spline.o -o spline

globals: globals.f90
	$(CC) $(CFLAGS1) globals.f90

spline: spline.f90 globals.f90
	$(CC) $(CFLAGS1) spline.f90

read_grid: read_grid.f90 globals.f90
	$(CC) $(CFLAGS1) read_grid.f90
	
read_input: read_input.f90 globals.f90
	$(CC) $(CFLAGS1) read_input.f90
	
connect: connect.f90 globals.f90
	$(CC) $(CFLAGS1) connect.f90
	
write_grid: write_grid.f90 globals.f90
	$(CC) $(CFLAGS1) write_grid.f90

.PHONY : clean
clean : 
	rm spline *.o *.mod	