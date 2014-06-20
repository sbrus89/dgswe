FC = ifort

VPATH := ../src

RK = -Drk33
LIB_PATH = ../metis-4.0.3/
LIB = metis
LAPACK_LIB = lapack

CFLAGS_VEC =  -O3 -xHost -vec-report6 -c
CFLAGS =  -O3 -xHost -c 
CFLAGS_LINK = -O3 -xHost -o 
ALIGN = #-DALIGN32

# ### OpenMP compilation ###
# CFLAGS_VEC =  -O3 -xHost -vec-report6 -openmp -Dopenmp -c
# CFLAGS =  -O3 -xHost -openmp -Dopenmp -c 
# CFLAGS_LINK = -O3 -xHost -openmp -Dopenmp -o 
# ALIGN = -DALIGN32

#  ### Mic compilation ###
# CFLAGS_VEC =  -O3 -mmic -vec-report6 -openmp -Dopenmp -c
# CFLAGS =  -O3 -mmic -openmp -Dopenmp -c 
# CFLAGS_LINK = -O3 -mmic -openmp -Dopenmp  -o 
# ALIGN = -DALIGN64

# ### Profile compilation ###
# CFLAGS_VEC =  -O3 -xHost -vec-report3 -profile-loops=all -c 
# CFLAGS =  -O3 -xHost -profile-loops=all -c 
# CFLAGS_LINK = -O3 -xHost -profile-loops=all -o 

# ### IPO compilation ###
# CFLAGS =  -O3 -xHost -ipo -c 
# CFLAGS_VEC =  -O3 -xHost -ipo -c 
# CFLAGS_LINK = -O3 -xHost -ipo -vec-report3 -o 

### Default compilation ###
# CFLAGS_VEC = -vec-report3 -c 
# CFLAGS = -c
# CFLAGS_LINK = -o

## Debug compilation ###
# CFLAGS_VEC = -c -C -g -traceback
# CFLAGS = -c -C -g -traceback
# CFLAGS_LINK = -o

### Other compiler flags ###
# -no-prec-div
# -opt-report3 -opt-report-file=opt_report 
# -vec-report3

# objects = globals.o basis.o area_qpts.o edge_qpts.o read_input.o read_grid.o read_forcing.o  \
# 	alloc_arrays.o ptr_arrays.o initial.o connect.o element_data.o edge_partition2.o  \
# 	rk.o rhs2.o rhs3.o rhs4.o numerical_flux.o swe_tri.o

objects = globals.o basis.o area_qpts.o edge_qpts.o read_input.o read_grid.o read_forcing.o  \
	alloc_arrays.o ptr_arrays.o initial.o connect.o element_data.o area_transformation.o edge_transformation.o edge_partition2.o  \
	area_integration.o edge_integration.o rhs2.o  numerical_flux.o rk.o swe_tri.o

dgswe : $(objects)
	$(FC) $(CFLAGS_LINK) dgswe $(objects) -L $(LIB_PATH) -l$(LIB) -l$(LAPACK_LIB)

globals.o : globals.F90
	$(FC) $(CFLAGS) $(ALIGN) $<

basis.o : basis.f90 globals.F90
	$(FC) $(CFLAGS) $<

#file_setup.o : file_setup.f90 globals.F90
#	$(FC) $(CFLAGS) -lnetcdf file_setup.f90

area_qpts.o : area_qpts.f90 globals.F90
	$(FC) $(CFLAGS) $<

edge_qpts.o : edge_qpts.f90 globals.F90
	$(FC) $(CFLAGS) $<

read_input.o : read_input.f90 globals.F90
	$(FC) $(CFLAGS) $<

read_grid.o : read_grid.f90 globals.F90
	$(FC) $(CFLAGS) $<

read_forcing.o : read_forcing.f90 globals.F90
	$(FC) $(CFLAGS) $<

alloc_arrays.o : alloc_arrays.f90 globals.F90
	$(FC) $(CFLAGS) $<

ptr_arrays.o : ptr_arrays.f90 globals.F90
	$(FC) $(CFLAGS) $<

initial.o : initial.f90 globals.F90
	$(FC) $(CFLAGS) $<

connect.o : connect.f90 globals.F90
	$(FC) $(CFLAGS) $<

element_data.o : element_data.f90 globals.F90
	$(FC) $(CFLAGS) $<

edge_partition2.o : edge_partition2.f90 globals.F90
	$(FC) $(CFLAGS) $<

area_transformation.o : area_transformation.F90 globals.F90 basis.f90
	$(FC) $(CFLAGS) $<

edge_transformation.o : edge_transformation.F90 globals.F90 basis.f90
	$(FC) $(CFLAGS) $<

area_integration.o : area_integration.F90 globals.F90
	$(FC) $(CFLAGS_VEC) $<

edge_integration.o : edge_integration.F90 globals.F90
	$(FC) $(CFLAGS_VEC) $<

rk.o : rk.F90 globals.F90
	$(FC) $(CFLAGS_VEC) $(RK) $<

rhs2.o : rhs2.f90 globals.F90
	$(FC) $(CFLAGS_VEC) $<

rhs3.o : rhs3.F90 globals.F90
	$(FC) $(CFLAGS_VEC) $<

rhs4.o : rhs4.F90 globals.F90
	$(FC) $(CFLAGS_VEC) $<	

numerical_flux.o : numerical_flux.f90 globals.F90
	$(FC) $(CFLAGS) $<

swe_tri.o : swe_tri.F90 globals.F90 basis.f90
	$(FC) $(CFLAGS) $(RK) $<

.PHONY : clean
clean : 
	rm dgswe *.o *.mod

.PHONY : metis
metis : 
	$(MAKE) -C ../metis-4.0.3/

.PHONY : cleanmetis
cleanmetis : 
	$(MAKE) -C ../metis-4.0.3/ clean

