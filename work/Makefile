VPATH := ../src 

CC = ifort
CFLAGS1 = -c -O3 -132  #-Dadcirc
#  CFLAGS1 = -c -C -g -132 -traceback -Dadcirc
CFLAGS2 = -O3 -o 
LIB = -llapack
#LIB = -Wl,--start-group   /afs/crc.nd.edu/x86_64_linux/intel/12.1/mkl/lib/intel64/libmkl_intel_lp64.a /afs/crc.nd.edu/x86_64_linux/intel/12.1/mkl/lib/intel64/libmkl_sequential.a  /afs/crc.nd.edu/x86_64_linux/intel/12.1/mkl/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm

GIT_SHA := $(shell git rev-parse HEAD)    
GIT_MOD := $(shell git diff-index --quiet HEAD; echo $$?)
GIT_BRANCH := $(shell git rev-parse --abbrev-ref HEAD )

GIT_MOD := $(strip $(subst 1,+,$(GIT_MOD)))
GIT_MOD := $(strip $(subst 0,,$(GIT_MOD)))

objects = globals.o version.o allocation.o basis.o evaluate.o read_input.o read_grid.o read_solution.o kdtree2.o find_nesting.o area_qpts.o error.o


main: version $(objects)
	@echo "Git SHA: $(GIT_SHA)"
	@echo "Modified: $(GIT_MOD) \n"
	$(CC) $(CFLAGS2) error $(objects) $(LIB)

version:
	@sed -i '7 c \      PRINT*, "  Branch: $(GIT_BRANCH)" ' ../src/version.F90
	@sed -i '8 c \      PRINT*, "  SHA: $(GIT_SHA) $(GIT_MOD)" ' ../src/version.F90

globals.o: globals.f90
	$(CC) $(CFLAGS1) $<
	
version.o: version.F90
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

error.o: error.F90 globals.f90 allocation.f90 evaluate.F90 basis.f90 kdtree2.F
	$(CC) $(CFLAGS1) $<

.PHONY : clean
clean : 
	rm error *.o *.mod	
