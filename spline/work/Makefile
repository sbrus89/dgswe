SHELL = /bin/sh

VPATH := ../src/
ODIR  = odir/

FC = ifort
FFLAGS = -O2 -C -g -traceback

FULLNAME := $(shell hostname)
NAME := $(strip $(subst ., , $(FULLNAME)))
FC_VERSION := $(shell $(FC) -v 2>&1)
DATE := $(shell date)
GIT_SHA := $(shell git rev-parse HEAD)
GIT_BRANCH := $(shell git rev-parse --abbrev-ref HEAD )
GIT_MOD := $(shell git diff-index --quiet HEAD; echo $$?)
GIT_MOD := $(strip $(subst 1,+,$(GIT_MOD)))
GIT_MOD := $(strip $(subst 0,,$(GIT_MOD)))
GIT_MFILES := $(strip $(shell git ls-files ../* -m))

FFLAGS += -I$(ODIR)

########################################################################
#                  Library Links (Platform Specific)                   #
########################################################################

NAME := $(strip $(subst ., , $(FULLNAME)))
ifeq ($(NAME),chl-tilos)
  LAPACK_LIB = -llapack
endif

ifeq ($(NAME),sbrus-laptop)
  LAPACK_LIB = -llapack
endif

ifeq ($(filter crc nd edu,$(NAME)),crc nd edu)
 LAPACK_LIB = -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a  $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm
endif

ifeq ($(filter tacc utexas edu,$(NAME)),tacc utexas edu)
  LAPACK_LIB =  -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm
endif

########################################################################
#                               Objects                                #
########################################################################

spline_objects = globals.o kdtree2.o lapack_interfaces.o basis.o allocation.o edge_connectivity_mod.o\
                 find_element_mod.o calc_spline.o check.o find_edge.o grid_file_mod.o\
                 read_grid.o read_input.o connect.o write_output.o version.o spline.o \
                 quit.o spherical_mod.o shape_functions_mod.o transformation.o vandermonde.o
                 
spline_obj := $(patsubst %.o, $(ODIR)%.o,$(spline_objects)) 


########################################################################
#                          Executable Targets                          #
########################################################################

.PHONY : spline

spline: version $(ODIR) $(spline_obj)
	@echo "\n\n"
	@echo "Git SHA: $(GIT_SHA)"
	@echo "Modified: $(GIT_MOD) \n"
	$(FC) $(FFLAGS) -o spline $(spline_obj) $(LAPACK_LIB)

########################################################################
#                      Premilinary Targets                             #
########################################################################	
.PHONY : $(ODIR) version

$(ODIR) :
	mkdir -p $@	
	
version :
	@sed -i '10 c \      gitBranch = "$(GIT_BRANCH)" ' ../../src/version.F90
	@sed -i '11 c \      gitSHA = "$(GIT_SHA) $(GIT_MOD)" ' ../../src/version.F90
	@sed -i '12 c \      compiler_version = "$(FC_VERSION)" ' ../../src/version.F90
	@sed -i '13 c \      compiler_flags = "$(FFLAGS)" ' ../../src/version.F90
	@sed -i '14 c \      modified_files = "$(GIT_MFILES)" ' ../../src/version.F90
	@sed -i '15 c \      compile_date = "$(DATE)" ' ../../src/version.F90
	@sed -i '16 c \      host = "$(NAME)" ' ../../src/version.F90


########################################################################
#                            Clean Targets                             #
########################################################################
.PHONY : clean

clean : 
	rm -r spline odir	

########################################################################
#                            File Targets                              #
########################################################################

# $< = first prerequisite/depenency
# $@ = name of the target	

$(ODIR)kdtree2.o : ../../src/kdtree2.F
	$(FC) $(FFLAGS) -132 -c $< -o $@
	mv kdtree2_precision_module.mod $(ODIR)
	mv kdtree2_module.mod $(ODIR)
	mv kdtree2_priority_queue_module.mod $(ODIR)
	
$(ODIR)lapack_interfaces.o : ../../src/lapack_interfaces.F90
	$(FC) $(FFLAGS) -c $< -o $@
	mv lapack_interfaces.mod $(ODIR)		

$(ODIR)globals.o : globals.f90                                $(ODIR)kdtree2.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv globals.mod $(ODIR)
	
$(ODIR)grid_file_mod.o : ../../src/grid_file_mod.F90          $(ODIR)globals.o $(ODIR)quit.o 
	$(FC) $(FFLAGS) -c $< -o $@
	mv grid_file_mod.mod $(ODIR)	
	
$(ODIR)vandermonde.o : ../../src/vandermonde.F90              $(ODIR)globals.o $(ODIR)basis.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv vandermonde.mod $(ODIR)		
	
$(ODIR)shape_functions_mod.o : ../../src/shape_functions_mod.F90 $(ODIR)globals.o $(ODIR)basis.o $(ODIR)lapack_interfaces.o $(ODIR)vandermonde.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv shape_functions_mod.mod $(ODIR)	

$(ODIR)transformation.o : ../../src/transformation.F90        $(ODIR)globals.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv transformation.mod $(ODIR)		

$(ODIR)allocation.o : allocation.f90                          $(ODIR)globals.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv allocation.mod $(ODIR)
	
$(ODIR)basis.o : ../../src/basis.F90                          $(ODIR)globals.o  $(ODIR)lapack_interfaces.o 
	$(FC) $(FFLAGS) -c $< -o $@
	mv basis.mod $(ODIR)
	
$(ODIR)calc_spline.o : calc_spline.f90                        $(ODIR)globals.o  $(ODIR)kdtree2.o $(ODIR)lapack_interfaces.o 
	$(FC) $(FFLAGS) -c $< -o $@
	mv calc_spline.mod $(ODIR)
	
$(ODIR)check.o : check.f90                                    $(ODIR)globals.o $(ODIR)calc_spline.o  $(ODIR)find_element_mod.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv check.mod $(ODIR)
	
$(ODIR)find_element_mod.o : ../../src/find_element_mod.F90    $(ODIR)globals.o  $(ODIR)kdtree2.o  $(ODIR)shape_functions_mod.o $(ODIR)transformation.o $(ODIR)vandermonde.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv find_element.mod $(ODIR)
	
$(ODIR)spherical_mod.o : ../../src/spherical_mod.F90          $(ODIR)globals.o 
	$(FC) $(FFLAGS) -c $< -o $@
	mv spherical_mod.mod $(ODIR)
	
$(ODIR)quit.o : ../../src/quit.F90                            $(ODIR)globals.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv quit.mod $(ODIR)	

$(ODIR)edge_connectivity_mod.o : ../../src/edge_connectivity_mod.F90 $(ODIR)globals.o $(ODIR)quit.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv edge_connectivity_mod.mod $(ODIR)	

$(ODIR)spline.o : spline.f90                                  $(ODIR)globals.o  $(ODIR)allocation.o  $(ODIR)calc_spline.o  $(ODIR)check.o $(ODIR)find_element_mod.o
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)read_grid.o : read_grid.f90                            $(ODIR)globals.o  $(ODIR)grid_file_mod.o $(ODIR)spherical_mod.o 
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)read_input.o : read_input.f90                          $(ODIR)globals.o
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)connect.o : connect.f90                                $(ODIR)globals.o $(ODIR)edge_connectivity_mod.o
	$(FC) $(FFLAGS) -c $< -o $@
	
$(ODIR)find_edge.o : find_edge.f90                            $(ODIR)globals.o
	$(FC) $(FFLAGS) -c $< -o $@	

$(ODIR)write_output.o : write_output.f90                      $(ODIR)globals.o
	$(FC) $(FFLAGS) -c $< -o $@
	
$(ODIR)version.o : ../../src/version.F90                      $(ODIR)globals.o
	$(FC) $(FFLAGS) -c $< -o $@
