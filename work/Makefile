FC = ifort

VPATH := ../src
ODIR = odir/

RK = -Drk33
LIB_PATH = ../metis-4.0.3/
LIB = metis
LAPACK_LIB = lapack

CFLAGS_VEC =  -O3 -xHost -vec-report6 -I$(ODIR) -c
CFLAGS =  -O3 -xHost -I$(ODIR) -c 
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
# ALIGN  = # -DALIGN32

# ### IPO compilation ###
# CFLAGS =  -O3 -xHost -ipo -c 
# CFLAGS_VEC =  -O3 -xHost -ipo -c 
# CFLAGS_LINK = -O3 -xHost -ipo -vec-report3 -o 

### Default compilation ###
# CFLAGS_VEC = -vec-report3 -c 
# CFLAGS = -c
# CFLAGS_LINK = -o

# Debug compilation ###
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

objects = globals.o allocation.o basis.o messenger2.o area_qpts.o edge_qpts.o read_input.o read_grid.o read_forcing.o modal2nodal.o \
	initial.o connect.o curvilinear.o element_data.o area_transformation.o edge_transformation.o edge_partition2.o  \
	area_integration.o edge_integration.o rhs2.o  numerical_flux.o rk.o swe_tri.o

obj := $(patsubst %.o, $(ODIR)%.o,$(objects))

# $< = first prerequisite/depenency
# $@ = name of the target


dgswe : $(ODIR) $(obj)
	$(FC) $(CFLAGS_LINK) dgswe $(obj) -L $(LIB_PATH) -l$(LIB) -l$(LAPACK_LIB)

.PHONY : $(ODIR)
$(ODIR) :
	mkdir -p $@

$(ODIR)globals.o : globals.F90
	$(FC) $(CFLAGS) $(ALIGN) $< -o $@
	mv globals.mod $(ODIR)

$(ODIR)allocation.o : allocation.f90 globals.F90
	$(FC) $(CFLAGS) $< -o $@
	mv allocation.mod $(ODIR)

$(ODIR)basis.o : basis.f90 globals.F90 allocation.f90
	$(FC) $(CFLAGS) $< -o $@
	mv basis.mod $(ODIR)

$(ODIR)messenger2.o : messenger2.F90 globals.F90
	$(FC) $(CFLAGS) $< -o $@
	mv messenger2.mod $(ODIR)

$(ODIR)area_qpts.o : area_qpts.f90 globals.F90 allocation.f90
	$(FC) $(CFLAGS) $< -o $@

$(ODIR)edge_qpts.o : edge_qpts.f90 globals.F90 allocation.f90
	$(FC) $(CFLAGS) $< -o $@

$(ODIR)read_input.o : read_input.F90 globals.F90
	$(FC) $(CFLAGS) $(RK) $< -o $@

$(ODIR)read_grid.o : read_grid.f90 globals.F90 allocation.f90
	$(FC) $(CFLAGS) $< -o $@

$(ODIR)read_forcing.o : read_forcing.f90 globals.F90 allocation.f90
	$(FC) $(CFLAGS) $< -o $@

$(ODIR)modal2nodal.o : modal2nodal.F90 globals.F90 basis.f90
	$(FC) $(CFLAGS) $< -o $@

$(ODIR)initial.o : initial.f90 globals.F90 allocation.f90
	$(FC) $(CFLAGS) $< -o $@

$(ODIR)connect.o : connect.f90 globals.F90 allocation.f90
	$(FC) $(CFLAGS) $< -o $@

$(ODIR)curvilinear.o : curvilinear.F90 globals.F90 basis.f90
	$(FC) $(CFLAGS) $< -o $@

$(ODIR)element_data.o : element_data.f90 globals.F90 allocation.f90
	$(FC) $(CFLAGS) $< -o $@

$(ODIR)edge_partition2.o : edge_partition2.f90 globals.F90 allocation.f90
	$(FC) $(CFLAGS) $< -o $@

$(ODIR)area_transformation.o : area_transformation.F90 globals.F90 basis.f90
	$(FC) $(CFLAGS) $< -o $@

$(ODIR)edge_transformation.o : edge_transformation.F90 globals.F90 basis.f90
	$(FC) $(CFLAGS) $< -o $@

$(ODIR)area_integration.o : area_integration.F90 globals.F90
	$(FC) $(CFLAGS_VEC) $< -o $@

$(ODIR)edge_integration.o : edge_integration.F90 globals.F90
	$(FC) $(CFLAGS_VEC) $< -o $@

$(ODIR)rk.o : rk.F90 globals.F90
	$(FC) $(CFLAGS_VEC) $(RK) $< -o $@

$(ODIR)rhs2.o : rhs2.f90 globals.F90
	$(FC) $(CFLAGS_VEC) $< -o $@

$(ODIR)rhs3.o : rhs3.F90 globals.F90
	$(FC) $(CFLAGS_VEC) $< -o $@

$(ODIR)rhs4.o : rhs4.F90 globals.F90
	$(FC) $(CFLAGS_VEC) $< -o $@

$(ODIR)numerical_flux.o : numerical_flux.f90 globals.F90
	$(FC) $(CFLAGS) $< -o $@

$(ODIR)swe_tri.o : swe_tri.F90 globals.F90 basis.f90 allocation.f90
	$(FC) $(CFLAGS) $(RK) $< -o $@

.PHONY : clean
clean : 
	rm -r dgswe odir

.PHONY : metis
metis : 
	$(MAKE) -C ../metis-4.0.3/

.PHONY : cleanmetis
cleanmetis : 
	$(MAKE) -C ../metis-4.0.3/ clean

