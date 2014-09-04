CC = ifort
CFLAGS1 = -c -O2 -C -g -traceback
CFLAGS2 = -O2
LIB = -llapack

all: main

main: globals read_grid read_input connect spline 
	$(CC) $(CFLAGS2) globals.o read_grid.o read_input.o connect.o spline.o -o spline

globals: globals.f90
	$(CC) $(CFLAGS1) globals.f90

spline: spline.f90
	$(CC) $(CFLAGS1) spline.f90

read_grid: read_grid.f90
	$(CC) $(CFLAGS1) read_grid.f90
	
read_input: read_input.f90
	$(CC) $(CFLAGS1) read_input.f90
	
connect: connect.f90
	$(CC) $(CFLAGS1) connect.f90

.PHONY : clean
clean : 
	rm spline *.o *.mod	