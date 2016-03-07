SHELL = /bin/sh

VPATH := ../src/

FULLNAME := $(shell hostname)

########################################################################
#                  Compiler Flags (Build Specific)                     #
########################################################################

ifeq ($(findstring debug,$(OPTS)),debug)
  FFLAGS += -C -g -traceback
  FLAG_TYPE += debug and bounds checking
endif
ifeq ($(findstring strict,$(OPTS)),strict)
  FFLAGS += -O0 -fp-model strict
  FLAG_TYPE += strict floating point 
endif
ifeq ($(findstring default,$(OPTS)),default)
  FFLAGS += -O2
  FLAG_TYPE += default  
endif

FFLAGS ?= -O3 -xHost
FLAG_TYPE ?= optimized

ifeq ($(findstring prof,$(OPTS)),prof)
  FFLAGS += -profile-loops=all
  FLAG_TYPE += profiling  
endif
ifeq ($(findstring align,$(OPTS)),align)
  FFLAGS += -align array32byte -align rec32byte
  FLAG_TYPE += aligned  
endif
ifeq ($(findstring ipo,$(OPTS)),ipo)
  FFLAGS += -ipo
  FLAG_TYPE += ipo  
endif

ifeq ($(findstring mpi,$(OPTS)),mpi)
  FFLAGS += -DCMPI
endif
ifeq ($(findstring omp,$(OPTS)),omp)
  FFLAGS += -openmp -Dopenmp
endif
# ifeq ($(findstring mic,$(OPTS)),mic)
#   FFLAGS += -mmic
# endif

ifeq ($(findstring netcdf,$(OPTS)),netcdf)
  FFLAGS += -DNETCDF
  NC_FLAG = true
endif







FFLAGS += -I$(ODIR)
OPT_REPORT = -vec_report6 
# OPT_REPORT = -qopt-report -qopt-report-phase=vec

ifeq ($(BUILDTYPE),dgswe)
  FC := ifort
  ODIR = odir_dgswe/

  EXE_NAME = dgswe
  MESSAGE = Compling serial program with $(FLAG_TYPE) flags
endif

ifeq ($(findstring mpi,$(OPTS)),mpi)
  FC := mpif90
  ODIR = odir_dgswe_mpi/
  
  EXE_NAME = dgswe_mpi
  MESSAGE = Compling MPI program with $(FLAG_TYPE) flags  
endif

ifeq ($(findstring omp,$(OPTS)),omp)
  FC := ifort
  ODIR = odir_dgswe_omp/
  
  EXE_NAME = dgswe_omp
  MESSAGE = Compling OpenMP program with $(FLAG_TYPE) flags  
endif

ifeq ($(BUILDTYPE),dgprep)
  FC := ifort
  ODIR = odir_dgprep/

  EXE_NAME = dgprep
  MESSAGE = Compling prep program with $(FLAG_TYPE) flags
endif

ifeq ($(BUILDTYPE),dgpost)
  FC := ifort
  ODIR = odir_dgpost/

  EXE_NAME = dgpost
  MESSAGE = Compling post program with $(FLAG_TYPE) flags
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
  NC_LIB = -lnetcdf -lnetcdff 
  NC_INC = -I/usr/local/include
endif

ifeq ($(NAME),sbrus-laptop)
  LAPACK_LIB = -llapack
endif

ifeq ($(filter crc nd edu,$(NAME)),crc nd edu)
    LAPACK_LIB = -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a  $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm
#   LAPACK_LIB = -L./ -llapack -lrefblas -ltmglib
#   NC_LIB = -L/afs/crc.nd.edu/x86_64_linux/netcdf/4.3.2/lib -lnetcdf -lnetcdff
#   NC_INC = -I/afs/crc.nd.edu/x86_64_linux/netcdf/4.3.2/include
endif

ifeq ($(filter tacc utexas edu,$(NAME)),tacc utexas edu)
  LAPACK_LIB =  -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm
#   NC_LIB = -L${TACC_NETCDF_LIB} -lnetcdf -lnetcdff
#   NC_INC = -I${TACC_NETCDF_INC}
endif

ifeq ($(NC_FLAG),true)
  NETCDF_LIB = $(NC_LIB)
  NETCDF_INC = $(NC_INC)
else 
  NETCDF_LIB = 
  NETCDF_INC = 
endif

METIS_LIB = -L../metis-4.0.3 -lmetis




########################################################################
#                               Objects                                #
########################################################################

dgswe_objects = globals.o lapack_interfaces.o messenger2.o allocation.o  \
                read_input.o read_grid.o read_forcing.o bathymetry_interp.o \
                connect.o initial.o interp_forcing.o curvilinear.o \
                edge_partition2.o metis2.o \
                vandermonde.o shape_functions.o transformation.o \
                element_data.o area_transformation.o normals.o edge_transformation.o \
                basis.o eval_basis.o area_qpts.o edge_qpts.o \
                rhs2.o rk.o \
                area_integration.o edge_integration.o  numerical_flux.o linear_solve.o \
                boundary_edge_land.o boundary_edge_flow.o boundary_edge_elev.o  \
                output.o quit.o modal2nodal.o version.o dgswe.o

dgswe_obj := $(patsubst %.o, $(ODIR)%.o,$(dgswe_objects))

dgprep_objects = globals.o messenger2.o read_input.o allocation.o \
                 read_grid.o read_forcing.o  curvilinear.o\
                 connect.o edge_qpts.o metis2.o decomp2.o \
                 basis.o vandermonde.o shape_functions.o transformation.o \
                 normals.o bathymetry_interp.o edge_transformation.o \
                 modal2nodal.o write_files.o version.o quit.o dgprep.o
                 
dgprep_obj := $(patsubst %.o, $(ODIR)%.o,$(dgprep_objects))  

dgpost_objects = globals.o messenger2.o version.o dgpost.o

dgpost_obj := $(patsubst %.o, $(ODIR)%.o,$(dgpost_objects))  

########################################################################
#                          Executable Targets                          #
########################################################################

.PHONY : $(ODIR) version dgswe 
.PHONY : metis dgprep dgpost 
.PHONY : clean cleanmetis cleandgprep cleandgpost

ifeq ($(MAKELEVEL),0)

dgswe : 
	@printf "\n\n"
	$(MAKE) BUILDTYPE=dgswe OPTS=$(OPTS) dgswe
	
metis : 
	@printf "\n\n"
	$(MAKE) -C ../metis-4.0.3/

dgprep : 
	@printf "\n\n"
	$(MAKE) BUILDTYPE=dgprep OPTS=$(OPTS) dgprep

dgpost : 
	@printf "\n\n"
	$(MAKE) BUILDTYPE=dgpost OPTS=$(OPTS) dgpost

else

FC_VERSION := $(shell $(FC) -v 2>&1 | head -2)
DATE := $(shell date)

# GIT_SHA := $(shell git rev-parse --short HEAD)
GIT_SHA := $(shell git rev-parse HEAD)
GIT_BRANCH := $(shell git rev-parse --abbrev-ref HEAD )
GIT_MOD := $(shell git diff-index --quiet HEAD; echo $$?)
GIT_MOD := $(strip $(subst 1,+,$(GIT_MOD)))
GIT_MOD := $(strip $(subst 0,,$(GIT_MOD)))
GIT_MFILES := $(strip $(shell git ls-files ../* -m))

dgswe : version $(ODIR)  $(dgswe_obj)
	@printf  " \n\033[35m$(MESSAGE) on $(FULLNAME) \033[0m \n"
	@printf "Git SHA: $(GIT_SHA) \n"
	@printf "Modified: $(GIT_MOD) \n\n"
	$(FC) $(FFLAGS) -o $(EXE_NAME) $(dgswe_obj) $(METIS_LIB) $(LAPACK_LIB) $(NETCDF_LIB)
	@printf "\n\n"
	
dgprep : version $(ODIR)  $(dgprep_obj)
	@printf " \n\033[35m$(MESSAGE) on $(FULLNAME) \033[0m \n"
	@printf "Git SHA: $(GIT_SHA) \n"
	@printf "Modified: $(GIT_MOD) \n\n"
	$(FC) $(FFLAGS) -o $(EXE_NAME) $(dgprep_obj) $(METIS_LIB) $(LAPACK_LIB)
	@printf "\n\n"
	
dgpost: version $(ODIR)  $(dgpost_obj)
	@printf " \n\033[35m$(MESSAGE) on $(FULLNAME) \033[0m \n"
	@printf "Git SHA: $(GIT_SHA) \n"
	@printf "Modified: $(GIT_MOD) \n\n"
	$(FC) $(FFLAGS) -o $(EXE_NAME) $(dgpost_obj)
	@printf "\n\n"


endif


########################################################################
#                      Preliminary Targets                             #
########################################################################

$(ODIR) :
	mkdir -p $@
	
version :
	@sed -i '10 c \      gitBranch = "$(GIT_BRANCH)" ' $(VPATH)version.F90
	@sed -i '11 c \      gitSHA = "$(GIT_SHA) $(GIT_MOD)" ' $(VPATH)version.F90
	@sed -i '12 c \      compiler_version = "$(FC_VERSION)" ' $(VPATH)version.F90
	@sed -i '13 c \      compiler_flags = "$(FFLAGS)" ' $(VPATH)version.F90
	@sed -i '14 c \      modified_files = "$(GIT_MFILES)" ' $(VPATH)version.F90
	@sed -i '15 c \      compile_date = "$(DATE)" ' $(VPATH)version.F90
	@sed -i '16 c \      host = "$(NAME)" ' $(VPATH)version.F90	
	
########################################################################
#                            Clean Targets                             #
########################################################################	
	
	
clean : 
	@printf "\n\n"
	rm odir*/*.o odir*/*.mod
	
clobber :
	@printf "\n\n"
	rm -rf odir*
	rm dgswe dgswe_mpi dgswe_omp

cleanmetis : 
	@printf "\n\n"
	$(MAKE) -C ../metis-4.0.3/ clean

cleandgprep : 
	@printf "\n\n"
	rm -r odir_dgprep/*.o odir_dgprep/*.mod

cleandgpost : 
	@printf "\n\n"
	$(MAKE) -C ../post clean
	
help : 
	@printf "\n\n"
	@printf "This Makefile uses commands of the form: \n"
	@printf "    make target OPTS=opt1,opt2,... \n"
	@printf "\n"
	@printf "The executable targets are: \n"
	@printf "    dgswe  : makes the main program executable \n"
	@printf "    dgprep : makes the domain decomposition executable \n"
	@printf "    dgpost : makes the solution globalization executable \n"
	@printf "    metis  : makes the METIS library used by dgprep \n"
	@printf "\n"
	@printf "The clean targets are: \n"
	@printf "    clean       : removes object files for dgswe \n"
	@printf "    clobber     : removes all object file directories and executables for dgswe \n"
	@printf "    cleandgprep : removes object files for dgprep \n"
	@printf "    cleandgpost : cleans the dgpost complication directory \n"
	@printf "    cleanmetis  : cleans the METIS compliation directory \n"
	@printf "\n"
	@printf "Valid options for the OPTS variable are: \n"
	@printf " - The following options affect the type of executable produced: \n"
	@printf "     mpi : compiles the MPI parallel version of dgswe \n"
	@printf "     omp : compiles the OpenMP paralel version of dgswe \n"
	@printf "     ** not specifying one of these options makes the serial version \n"
	@printf "\n"
	@printf " - The following options specify the un-optimized compiler flags used: \n"
	@printf "     debug   : compiles with bounds checking and debugging flags \n"
	@printf "     strict  : compiles with strict floating point flags \n"
	@printf "     default : compiles with default ifort flags \n"
	@printf "     ** not specifying one of these options defaults to: -O3 -xHost \n"
	@printf "\n"
	@printf " - The following options specify the optimized compiler flags used: \n"
	@printf "     align : compiles with array and derived type alignment flags \n" 
	@printf "     prof  : compiles with profiling flags \n"
	@printf "     ipo   : complies with inter procedural optimization flags \n"



########################################################################
#                            File Targets                              #
########################################################################

# $< = first prerequisite/depenency
# $@ = name of the target	

$(ODIR)globals.o : globals.F90
	$(FC) $(FFLAGS) -c $< -o $@
	mv globals.mod $(ODIR)

$(ODIR)lapack_interfaces.o : lapack_interfaces.F90
	$(FC) $(FFLAGS) -c $< -o $@
	mv lapack_interfaces.mod $(ODIR)

$(ODIR)messenger2.o : messenger2.F90                           $(ODIR)globals.o $(ODIR)quit.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv messenger2.mod $(ODIR)

$(ODIR)allocation.o : allocation.F90                           $(ODIR)globals.o $(ODIR)messenger2.o $(ODIR)read_input.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv allocation.mod $(ODIR)

$(ODIR)basis.o : basis.F90                                     $(ODIR)globals.o $(ODIR)lapack_interfaces.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv basis.mod $(ODIR)
	
$(ODIR)read_input.o : read_input.F90                           $(ODIR)globals.o $(ODIR)messenger2.o $(ODIR)quit.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv read_dginp.mod $(ODIR)

$(ODIR)vandermonde.o : vandermonde.F90                         $(ODIR)globals.o $(ODIR)basis.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv vandermonde.mod $(ODIR)
	
$(ODIR)shape_functions.o : shape_functions.F90                 $(ODIR)globals.o $(ODIR)messenger2.o $(ODIR)basis.o $(ODIR)lapack_interfaces.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv shape_functions.mod $(ODIR)

$(ODIR)transformation.o : transformation.F90                   $(ODIR)globals.o $(ODIR)read_input.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv transformation.mod $(ODIR)
	
$(ODIR)output.o : output.F90                                   $(ODIR)globals.o $(ODIR)messenger2.o $(ODIR)read_input.o
	$(FC) $(FFLAGS) $(NETCDF_INC) -c $< -o $@
	mv output.mod $(ODIR)
	
$(ODIR)quit.o : quit.F90                                       $(ODIR)globals.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv quit.mod $(ODIR)
	

$(ODIR)area_qpts.o : area_qpts.F90                             $(ODIR)globals.o $(ODIR)messenger2.o $(ODIR)read_input.o $(ODIR)allocation.o
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)edge_qpts.o : edge_qpts.F90                             $(ODIR)globals.o $(ODIR)messenger2.o $(ODIR)read_input.o $(ODIR)allocation.o
	$(FC) $(FFLAGS) -c $< -o $@
	
$(ODIR)eval_basis.o : eval_basis.F90                           $(ODIR)globals.o $(ODIR)messenger2.o $(ODIR)read_input.o $(ODIR)allocation.o $(ODIR)basis.o
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)read_grid.o : read_grid.F90                             $(ODIR)globals.o $(ODIR)messenger2.o $(ODIR)read_input.o $(ODIR)allocation.o $(ODIR)quit.o
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)read_forcing.o : read_forcing.F90                       $(ODIR)globals.o $(ODIR)messenger2.o $(ODIR)read_input.o $(ODIR)allocation.o $(ODIR)quit.o
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)modal2nodal.o : modal2nodal.F90                         $(ODIR)globals.o  $(ODIR)read_input.o $(ODIR)basis.o
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)interp_forcing.o : interp_forcing.F90                   $(ODIR)globals.o $(ODIR)allocation.o
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)initial.o : initial.F90                                 $(ODIR)globals.o $(ODIR)read_input.o $(ODIR)allocation.o $(ODIR)basis.o
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)connect.o : connect.F90                                 $(ODIR)globals.o $(ODIR)messenger2.o $(ODIR)read_input.o $(ODIR)allocation.o
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)curvilinear.o : curvilinear.F90                         $(ODIR)globals.o $(ODIR)read_input.o $(ODIR)basis.o $(ODIR)transformation.o
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)element_data.o : element_data.F90                       $(ODIR)globals.o $(ODIR)read_input.o $(ODIR)allocation.o $(ODIR)basis.o $(ODIR)vandermonde.o $(ODIR)shape_functions.o
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)edge_partition2.o : edge_partition2.F90                 $(ODIR)globals.o $(ODIR)messenger2.o $(ODIR)read_input.o $(ODIR)allocation.o
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)metis2.o : metis2.F90                                   $(ODIR)globals.o $(ODIR)messenger2.o 
	$(FC) $(FFLAGS) -c $< -o $@	
	
$(ODIR)area_transformation.o : area_transformation.F90         $(ODIR)globals.o $(ODIR)transformation.o $(ODIR)lapack_interfaces.o
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)edge_transformation.o : edge_transformation.F90         $(ODIR)globals.o $(ODIR)basis.o
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)area_integration.o : area_integration.F90               $(ODIR)globals.o $(ODIR)read_input.o
	$(FC) $(FFLAGS) $(OPT_REPORT) -c $< -o $@
	
$(ODIR)normals.o : normals.F90                                 $(ODIR)globals.o $(ODIR)transformation.o
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)edge_integration.o : edge_integration.F90               $(ODIR)globals.o $(ODIR)messenger2.o
	$(FC) $(FFLAGS) $(OPT_REPORT) -c $< -o $@
	
$(ODIR)bathymetry_interp.o : bathymetry_interp.F90             $(ODIR)globals.o $(ODIR)lapack_interfaces.o $(ODIR)read_input.o $(ODIR)transformation.o
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)rk.o : rk.F90                                           $(ODIR)globals.o $(ODIR)messenger2.o $(ODIR)read_input.o $(ODIR)output.o $(ODIR)quit.o
	$(FC) $(FFLAGS) $(OPT_REPORT) -c $< -o $@

$(ODIR)rhs2.o : rhs2.F90                                       $(ODIR)globals.o $(ODIR)messenger2.o $(ODIR)read_input.o
	$(FC) $(FFLAGS) $(OPT_REPORT) -c $< -o $@

$(ODIR)rhs3.o : rhs3.F90                                       $(ODIR)globals.o $(ODIR)messenger2.o $(ODIR)read_input.o
	$(FC) $(FFLAGS) $(OPT_REPORT) -c $< -o $@

$(ODIR)numerical_flux.o : numerical_flux.F90                   $(ODIR)globals.o
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)boundary_edge_land.o : boundary_edge_land.F90           $(ODIR)globals.o
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)boundary_edge_flow.o : boundary_edge_flow.F90           $(ODIR)globals.o
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)boundary_edge_elev.o : boundary_edge_elev.F90           $(ODIR)globals.o
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)linear_solve.o : linear_solve.F90                       $(ODIR)globals.o
	$(FC) $(FFLAGS) -c $< -o $@	
	
$(ODIR)version.o : version.F90                                 $(ODIR)globals.o $(ODIR)messenger2.o 
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)dgswe.o : dgswe.F90                                     $(ODIR)globals.o $(ODIR)messenger2.o $(ODIR)read_input.o $(ODIR)allocation.o $(ODIR)basis.o $(ODIR)output.o 
	$(FC) $(FFLAGS) -c $< -o $@

# prep file targets

$(ODIR)decomp2.o : ../prep/src/decomp2.F90                     $(ODIR)globals.o $(ODIR)read_input.o $(ODIR)messenger2.o 
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)write_files.o : ../prep/src/write_files.F90             $(ODIR)globals.o $(ODIR)read_input.o $(ODIR)messenger2.o 
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)dgprep.o : ../prep/src/dgprep.F90                        $(ODIR)globals.o $(ODIR)allocation.o $(ODIR)read_input.o $(ODIR)messenger2.o $(ODIR)vandermonde.o $(ODIR)shape_functions.o
	$(FC) $(FFLAGS) -c $< -o $@

# post file targets

$(ODIR)dgpost.o : ../post/dgpost.F90                        
	$(FC) $(FFLAGS) -c $< -o $@
