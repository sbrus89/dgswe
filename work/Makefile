FC = ifort

VPATH := ../src

CFLAGS_VEC =  -O3 -xHost -vec-report6 -fno-alias -c
CFLAGS =  -O3 -xHost -c 
CFLAGS_LINK = -O3 -xHost -o 

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

# ## Debug compilation ###
# CFLAGS_VEC = -c -C -g -traceback
# CFLAGS = -c -C -g -traceback
# CFLAGS_LINK = -o

### Other compiler flags ###
# -no-prec-div
# -opt-report3 -opt-report-file=opt_report 
# -vec-report3

objects = globals.o basis.o area_qpts.o edge_qpts.o read_grid.o read_forcing.o alloc_arrays.o ptr_arrays.o initial.o connect.o element_data.o rk.o rhs2.o numerical_flux.o swe_tri.o

dgswe : $(objects)
	$(FC) $(CFLAGS_LINK) dgswe $(objects)

globals.o : globals.f90
	$(FC) $(CFLAGS) $<

basis.o : basis.f90 globals.f90
	$(FC) $(CFLAGS) $<

#file_setup.o : file_setup.f90 globals.f90
#	$(FC) $(CFLAGS) -lnetcdf file_setup.f90

area_qpts.o : area_qpts.f90 globals.f90
	$(FC) $(CFLAGS) $<

edge_qpts.o : edge_qpts.f90 globals.f90
	$(FC) $(CFLAGS) $<

read_grid.o : read_grid.f90 globals.f90
	$(FC) $(CFLAGS) $<

read_forcing.o : read_forcing.f90 globals.f90
	$(FC) $(CFLAGS) $<

alloc_arrays.o : alloc_arrays.f90 globals.f90
	$(FC) $(CFLAGS) $<

ptr_arrays.o : ptr_arrays.f90 globals.f90
	$(FC) $(CFLAGS) $<

initial.o : initial.f90 globals.f90
	$(FC) $(CFLAGS) $<

connect.o : connect.f90 globals.f90
	$(FC) $(CFLAGS) $<

element_data.o : element_data.f90 globals.f90
	$(FC) $(CFLAGS) $<

rk.o : rk.F90 globals.f90
	$(FC) $(CFLAGS_VEC) -Drk33 $<

rhs2.o : rhs2.f90 globals.f90
	$(FC) $(CFLAGS_VEC) $<

numerical_flux.o : numerical_flux.f90 globals.f90
	$(FC) $(CFLAGS) $<

swe_tri.o : swe_tri.f90 globals.f90
	$(FC) $(CFLAGS_VEC) $<

.PHONY : clean
clean : 
	rm dgswe $(objects) *.mod

