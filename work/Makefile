SHELL = /bin/sh

RK = -Drk33


VPATH := ../src
ODIR = odir/

NAME := $(shell hostname)
$(warning (INFO) Machine hostname is $(NAME))

########################################################################
#                  Compiler Flags (Build Specific)                     #
########################################################################

ifeq ($(BUILDTYPE),dgswe)
  FC = ifort
  CFLAGS_VEC  =  -O3 -xHost -vec-report6 -I$(ODIR) -c
  CFLAGS      =  -O3 -xHost -I$(ODIR) -c 
  CFLAGS_LINK = -O3 -xHost -o 
  $(warning (INFO) Compling serial program with optimization flags)  
endif

ifeq ($(BUILDTYPE),pdgswe)
  FC = mpif90
  CFLAGS_VEC  =  -O3 -xHost -vec-report6 -DCMPI -I$(ODIR) -c
  CFLAGS      =  -O3 -xHost -DCMPI -I$(ODIR) -c 
  CFLAGS_LINK = -O3 -xHost -o 
  $(warning (INFO) Compling parallel program with optimization flags)    
endif

ifeq ($(BUILDTYPE),dgswe_debug)
  FC = ifort
  CFLAGS_VEC  =  -C -g -traceback -I$(ODIR) -c
  CFLAGS      =  -C -g -traceback -I$(ODIR) -c 
  CFLAGS_LINK = -o 
  $(warning (INFO) Compling serial program with debug flags)    
endif

ifeq ($(BUILDTYPE),pdgswe_debug)
  FC = mpif90
  CFLAGS_VEC  =  -C -g -traceback -DCMPI -I$(ODIR) -c
  CFLAGS      =  -C -g -traceback -DCMPI -I$(ODIR) -c 
  CFLAGS_LINK =  -o 
  $(warning (INFO) Compling parallel program with debug flags)    
endif

# ### OpenMP compilation ###
# CFLAGS_VEC =  -O3 -xHost -vec-report6 -openmp -Dopenmp -c
# CFLAGS =  -O3 -xHost -openmp -Dopenmp -c 
# CFLAGS_LINK = -O3 -xHost -openmp -Dopenmp -o 

#  ### Mic compilation ###
# CFLAGS_VEC =  -O3 -mmic -vec-report6 -openmp -Dopenmp -c
# CFLAGS =  -O3 -mmic -openmp -Dopenmp -c 
# CFLAGS_LINK = -O3 -mmic -openmp -Dopenmp  -o 

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

### Other compiler flags ###
# -no-prec-div
# -opt-report3 -opt-report-file=opt_report 
# -vec-report3

########################################################################
#                  Library Links (Platform Specific)                   #
########################################################################

NAME := $(strip $(subst ., , $(NAME)))
ifeq ($(NAME),chl-tilos)
  LAPACK_LIB = -llapack
  ALIGN = #-DALIGN64
endif

ifeq ($(filter crc nd edu,$(NAME)),crc nd edu)
  LAPACK_LIB = -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a  $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm
  ALIGN = #-DALIGN32
endif

ifeq ($(filter tacc utexas edu,$(NAME)),tacc utexas edu)
  LAPACK_LIB =  -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm
  ALIGN = #-DALIGN64
endif

METIS_LIB = -L../metis-4.0.3 -lmetis

########################################################################
#                               Objects                                #
########################################################################

objects = globals.o messenger2.o allocation.o basis.o  area_qpts.o edge_qpts.o read_input.o read_grid.o read_forcing.o modal2nodal.o \
	initial.o connect.o curvilinear.o element_data.o area_transformation.o edge_transformation.o edge_partition2.o  \
	area_integration.o edge_integration.o rhs2.o  numerical_flux.o rk.o swe_tri.o

obj := $(patsubst %.o, $(ODIR)%.o,$(objects))

########################################################################
#                          Executable Targets                          #
########################################################################

.PHONY : $(ODIR) dgswe pdgswe dgswe_debug pdgswe_debug clean metis cleanmetis prep cleanprep post cleanpost

ifeq ($(MAKELEVEL),0)

dgswe : 
	$(MAKE) BUILDTYPE=dgswe $@

pdgswe : 
	$(MAKE) BUILDTYPE=pdgswe $@

dgswe_debug : 
	$(MAKE) BUILDTYPE=dgswe_debug $@

pdgswe_debug : 
	$(MAKE) BUILDTYPE=pdgswe_debug $@

metis : 
	$(MAKE) -C ../metis-4.0.3/

prep : 
	$(MAKE) -C ../prep/work/ dgprep

post : 
	$(MAKE) -C ../post dgpost	

else

dgswe : $(ODIR) $(obj)
	$(FC) $(CFLAGS_LINK) dgswe $(obj) $(METIS_LIB) $(LAPACK_LIB)

pdgswe : $(ODIR) $(obj)
	$(FC) $(CFLAGS_LINK) pdgswe $(obj) $(METIS_LIB) $(LAPACK_LIB)

dgswe_debug : $(ODIR) $(obj)
	$(FC) $(CFLAGS_LINK) dgswe_debug $(obj) $(METIS_LIB) $(LAPACK_LIB)

pdgswe_debug : $(ODIR) $(obj)
	$(FC) $(CFLAGS_LINK) pdgswe_debug $(obj) $(METIS_LIB) $(LAPACK_LIB)	

endif


########################################################################
#                            Clean Targets                             #
########################################################################

.PHONY : $(ODIR)
$(ODIR) :
	mkdir -p $@
	
clean : 
	rm -r pdgswe* dgswe* odir

cleanmetis : 
	$(MAKE) -C ../metis-4.0.3/ clean

cleanprep : 
	$(MAKE) -C ../prep/work/ clean

cleanpost : 
	$(MAKE) -C ../post clean

########################################################################
#                            File Targets                              #
########################################################################

# $< = first prerequisite/depenency
# $@ = name of the target	

$(ODIR)globals.o : globals.F90
	$(FC) $(CFLAGS) $(ALIGN) $< -o $@
	mv globals.mod $(ODIR)

$(ODIR)messenger2.o : messenger2.F90 globals.F90
	$(FC) $(CFLAGS) $(ALIGN) $< -o $@
	mv messenger2.mod $(ODIR)	

$(ODIR)allocation.o : allocation.f90 globals.F90 messenger2.F90
	$(FC) $(CFLAGS) $< -o $@
	mv allocation.mod $(ODIR)

$(ODIR)basis.o : basis.f90 globals.F90 allocation.f90
	$(FC) $(CFLAGS) $< -o $@
	mv basis.mod $(ODIR)

$(ODIR)area_qpts.o : area_qpts.f90 globals.F90 allocation.f90
	$(FC) $(CFLAGS) $< -o $@

$(ODIR)edge_qpts.o : edge_qpts.f90 globals.F90 allocation.f90
	$(FC) $(CFLAGS) $< -o $@

$(ODIR)read_input.o : read_input.F90 globals.F90 messenger2.F90
	$(FC) $(CFLAGS) $(RK) $< -o $@

$(ODIR)read_grid.o : read_grid.f90 globals.F90 allocation.f90 messenger2.F90
	$(FC) $(CFLAGS) $< -o $@

$(ODIR)read_forcing.o : read_forcing.f90 globals.F90 allocation.f90
	$(FC) $(CFLAGS) $< -o $@

$(ODIR)modal2nodal.o : modal2nodal.F90 globals.F90 basis.f90
	$(FC) $(CFLAGS) $< -o $@

$(ODIR)initial.o : initial.f90 globals.F90 allocation.f90
	$(FC) $(CFLAGS) $< -o $@

$(ODIR)connect.o : connect.F90 globals.F90 allocation.f90
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

$(ODIR)rhs2.o : rhs2.F90 globals.F90 messenger2.F90
	$(FC) $(CFLAGS_VEC) $< -o $@

$(ODIR)rhs3.o : rhs3.F90 globals.F90
	$(FC) $(CFLAGS_VEC) $< -o $@

$(ODIR)rhs4.o : rhs4.F90 globals.F90
	$(FC) $(CFLAGS_VEC) $< -o $@

$(ODIR)numerical_flux.o : numerical_flux.f90 globals.F90
	$(FC) $(CFLAGS) $< -o $@

$(ODIR)swe_tri.o : swe_tri.F90 globals.F90 basis.f90 allocation.f90 messenger2.F90
	$(FC) $(CFLAGS) $(RK) $< -o $@