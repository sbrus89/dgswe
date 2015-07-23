SHELL = /bin/sh

VPATH := ../src/
ODIR = odir/

FULLNAME := $(shell hostname)
# GIT_SHA := $(shell git rev-parse --short HEAD)
GIT_SHA := $(shell git rev-parse HEAD)
GIT_MOD := $(shell git diff-index --quiet HEAD; echo $$?)
GIT_BRANCH := $(shell git rev-parse --abbrev-ref HEAD )

GIT_MOD := $(strip $(subst 1,+,$(GIT_MOD)))
GIT_MOD := $(strip $(subst 0,,$(GIT_MOD)))

########################################################################
#                  Compiler Flags (Build Specific)                     #
########################################################################

ifeq ($(BUILDTYPE),dgswe)
  FC = ifort
  
#   CFLAGS_VEC  =  -O3 -xHost -align array32byte -align rec32byte -vec-report6  -I$(ODIR) -c
#   CFLAGS      =  -O3 -xHost -align array32byte -align rec32byte  -I$(ODIR) -c 
#   MESSAGE = Compling serial program with optimization flags  
#   CFLAGS_LINK = -O3 -xHost -o 
  
  CFLAGS_VEC  =  -O0  -I$(ODIR) -c
  CFLAGS      =  -O0  -I$(ODIR) -c 
  MESSAGE = Compling serial program without optimization flags  
  CFLAGS_LINK = -o   

#   ### Default compilation ###
#   CFLAGS_VEC = -vec-report3 -I$(ODIR) -c 
#   CFLAGS = -I$(ODIR) -c
#   CFLAGS_LINK = -o

#   ### IPO compilation ###
#   CFLAGS =  -O3 -xHost -ipo -I$(ODIR) -c 
#   CFLAGS_VEC =  -O3 -xHost -ipo -I$(ODIR) -c 
#   CFLAGS_LINK = -O3 -xHost -ipo -vec-report3 -I$(ODIR) -o 
  
#   ### Profile compilation ###
#   CFLAGS_VEC =  -O3 -xHost -align array32byte -vec-report3 -profile-loops=all -I$(ODIR) -c 
#   CFLAGS =  -O3 -xHost -align array32byte -profile-loops=all -I$(ODIR) -c 
#   CFLAGS_LINK = -O3 -xHost -profile-loops=all -I$(ODIR) -o

  EXE_NAME = dgswe
  
endif

ifeq ($(BUILDTYPE),dgswe_omp)
  FC = ifort
  CFLAGS_VEC  =  -O3 -xHost -align array32byte -openmp -Dopenmp -vec-report6 -I$(ODIR) -c
  CFLAGS      =  -O3 -xHost -align array32byte -openmp -Dopenmp -I$(ODIR) -c 
  CFLAGS_LINK = -O3 -xHost -openmp -Dopenmp -o 
  
#   ### Mic compilation ###
#   CFLAGS_VEC =  -O3 -mmic -vec-report6 -openmp -Dopenmp -c
#   CFLAGS =  -O3 -mmic -openmp -Dopenmp -c 
#   CFLAGS_LINK = -O3 -mmic -openmp -Dopenmp  -o   
  
  EXE_NAME = dgswe_omp
  MESSAGE = Compling OpenMP program with optimization flags
endif

ifeq ($(BUILDTYPE),pdgswe)
  FC = mpif90
  
#   CFLAGS_VEC  =  -O3 -xHost -vec-report6 -DCMPI -I$(ODIR) -c
#   CFLAGS      =  -O3 -xHost -DCMPI -I$(ODIR) -c 
#   CFLAGS_LINK = -O3 -xHost -o 
#   MESSAGE = Compling parallel program with optimization flags
  
  CFLAGS_VEC  =  -O0 -DCMPI -I$(ODIR) -c
  CFLAGS      =  -O0 -DCMPI -I$(ODIR) -c 
  CFLAGS_LINK =  -o 
  MESSAGE = Compling parallel program without optimization flags  
  
  EXE_NAME = pdgswe
endif

ifeq ($(BUILDTYPE),dgswe_debug)
  FC = ifort
  CFLAGS_VEC  =  -C -g -traceback -I$(ODIR) -c
  CFLAGS      =  -C -g -traceback -I$(ODIR) -c 
  CFLAGS_LINK = -o 
  EXE_NAME = dgswe_debug
  MESSAGE =  Compling serial program with debug flags 
endif

ifeq ($(BUILDTYPE),pdgswe_debug)
  FC = mpif90
  CFLAGS_VEC  =  -C -g -traceback -DCMPI -I$(ODIR) -c
  CFLAGS      =  -C -g -traceback -DCMPI -I$(ODIR) -c 
  CFLAGS_LINK =  -o 
  EXE_NAME = pdgswe_debug
  MESSAGE = Compling parallel program with debug flags
endif

### Other compiler flags ###
# -no-prec-div
# -opt-report3 -opt-report-file=opt_report 
# -vec-report3

########################################################################
#                  Library Links (Platform Specific)                   #
########################################################################

NAME := $(strip $(subst ., , $(FULLNAME)))
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

objects = globals.o messenger2.o read_input.o allocation.o basis.o  area_qpts.o edge_qpts.o read_grid.o read_forcing.o modal2nodal.o \
	initial.o connect.o curvilinear.o element_data.o shape_functions.o vandermonde.o area_transformation.o edge_transformation.o \
	bathymetry_interp.o edge_partition2.o area_integration.o edge_integration.o rhs2.o  numerical_flux.o rk.o write_output.o version.o swe_tri.o

obj := $(patsubst %.o, $(ODIR)%.o,$(objects))

########################################################################
#                          Executable Targets                          #
########################################################################

.PHONY : $(ODIR) dgswe dgswe_omp pdgswe dgswe_debug pdgswe_debug clean metis cleanmetis prep cleanprep post cleanpost

ifeq ($(MAKELEVEL),0)

dgswe : 
	@echo "\n\n"
	$(MAKE) BUILDTYPE=dgswe dgswe
	
dgswe_omp : 
	@echo "\n\n"
	$(MAKE) BUILDTYPE=dgswe_omp dgswe	

pdgswe : 
	@echo "\n\n"
	$(MAKE) BUILDTYPE=pdgswe dgswe

dgswe_debug : 
	@echo "\n\n"
	$(MAKE) BUILDTYPE=dgswe_debug dgswe

pdgswe_debug : 
	@echo "\n\n"
	$(MAKE) BUILDTYPE=pdgswe_debug dgswe

metis : 
	@echo "\n\n"
	$(MAKE) -C ../metis-4.0.3/

prep : 
	@echo "\n\n"
	$(MAKE) -C ../prep/work/ dgprep

post : 
	@echo "\n\n"
	$(MAKE) -C ../post dgpost	

else

dgswe : $(ODIR)  $(obj)
	@echo  " \n\033[35m$(MESSAGE) on $(FULLNAME) \033[0m \n"
	@echo "Git SHA: $(GIT_SHA)"
	@echo "Modified: $(GIT_MOD)"
	$(FC) $(CFLAGS_LINK) $(EXE_NAME) $(obj) $(METIS_LIB) $(LAPACK_LIB)
endif


########################################################################
#                            Clean Targets                             #
########################################################################

$(ODIR) :
	mkdir -p $@
	@sed -i '8 c \      gitBranch = "$(GIT_BRANCH)" ' $(VPATH)version.F90
	@sed -i '9 c \      gitSHA = "$(GIT_SHA) $(GIT_MOD)" ' $(VPATH)version.F90
	
clean : 
	@echo "\n\n"
	rm -r odir
	
clobber :
	@echo "\n\n"
	rm -r odir
	rm dgswe pdgswe dgswe_debug pdgswe_debug

cleanmetis : 
	@echo "\n\n"
	$(MAKE) -C ../metis-4.0.3/ clean

cleanprep : 
	@echo "\n\n"
	$(MAKE) -C ../prep/work/ clean

cleanpost : 
	@echo "\n\n"
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

$(ODIR)basis.o : basis.f90 globals.F90 allocation.f90 messenger2.F90
	$(FC) $(CFLAGS) $< -o $@
	mv basis.mod $(ODIR)
	
$(ODIR)read_input.o : read_input.F90 globals.F90 messenger2.F90
	$(FC) $(CFLAGS) $< -o $@	
	mv read_dginp.mod $(ODIR)		

$(ODIR)area_qpts.o : area_qpts.f90 globals.F90 allocation.f90 messenger2.F90
	$(FC) $(CFLAGS) $< -o $@

$(ODIR)edge_qpts.o : edge_qpts.f90 globals.F90 allocation.f90 messenger2.F90
	$(FC) $(CFLAGS) $< -o $@

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

$(ODIR)edge_partition2.o : edge_partition2.f90 globals.F90 allocation.f90 messenger2.F90
	$(FC) $(CFLAGS) $< -o $@

$(ODIR)vandermonde.o : vandermonde.F90 globals.F90 basis.f90
	$(FC) $(CFLAGS) $< -o $@
	
$(ODIR)shape_functions.o : shape_functions.F90 globals.F90 basis.f90
	$(FC) $(CFLAGS) $< -o $@
	
$(ODIR)area_transformation.o : area_transformation.F90 globals.F90 basis.f90
	$(FC) $(CFLAGS) $< -o $@

$(ODIR)edge_transformation.o : edge_transformation.F90 globals.F90 basis.f90
	$(FC) $(CFLAGS) $< -o $@

$(ODIR)area_integration.o : area_integration.F90 globals.F90
	$(FC) $(CFLAGS_VEC) $< -o $@

$(ODIR)edge_integration.o : edge_integration.F90 globals.F90
	$(FC) $(CFLAGS_VEC) $< -o $@
	
$(ODIR)bathymetry_interp.o : bathymetry_interp.F90 globals.F90
	$(FC) $(CFLAGS_VEC) $< -o $@	

$(ODIR)rk.o : rk.F90 globals.F90
	$(FC) $(CFLAGS_VEC) $< -o $@

$(ODIR)rhs2.o : rhs2.F90 globals.F90 messenger2.F90
	$(FC) $(CFLAGS_VEC) $< -o $@

$(ODIR)rhs3.o : rhs3.F90 globals.F90
	$(FC) $(CFLAGS_VEC) $< -o $@

$(ODIR)rhs4.o : rhs4.F90 globals.F90
	$(FC) $(CFLAGS_VEC) $< -o $@

$(ODIR)numerical_flux.o : numerical_flux.f90 globals.F90
	$(FC) $(CFLAGS) $< -o $@
	
$(ODIR)write_output.o : write_output.F90 globals.F90 messenger2.F90
	$(FC) $(CFLAGS) $< -o $@	
	
$(ODIR)version.o : version.F90 globals.F90 messenger2.F90
	$(FC) $(CFLAGS) $< -o $@		

$(ODIR)swe_tri.o : swe_tri.F90 globals.F90 basis.f90 allocation.f90 messenger2.F90
	$(FC) $(CFLAGS) $< -o $@