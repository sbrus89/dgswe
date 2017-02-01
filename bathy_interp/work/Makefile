SHELL = /bin/sh

VPATH := ../src/
ODIR = odir/

FC = ifort
FFLAGS = -O3 -132 -I$(ODIR) #-Dadcirc
# FFLAGS = -O0 -fp-model strict -132 -I$(ODIR) #-Dadcirc
# FFLAGS = -132 -C -g -traceback -I$(ODIR) #-Dadcirc

FULLNAME := $(shell hostname)
NAME := $(strip $(subst ., , $(FULLNAME)))
FC_VERSION := $(shell $(FC) -v 2>&1)
DATE := $(shell date)
GIT_CHK := $(shell git checkout -- ../../src/version.F90)
GIT_SHA := $(shell git rev-parse HEAD)
GIT_BRANCH := $(shell git rev-parse --abbrev-ref HEAD )
GIT_MOD := $(shell git diff-index --quiet HEAD; echo $$?)
GIT_MOD := $(strip $(subst 1,+,$(GIT_MOD)))
GIT_MOD := $(strip $(subst 0,,$(GIT_MOD)))
GIT_MFILES := $(strip $(shell git ls-files ../* -m))

########################################################################
#                  Library Links (Platform Specific)                   #
########################################################################

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

objects = globals.o kdtree2.o allocation.o basis.o  \
          connect.o curvilinear.o evaluate.o find_element_mod.o sort_mod.o \
          read_input.o read_grid.o write_results.o system.o \
          bathy_interp.o bathymetry_interp_mod.o lapack_interfaces.o interp_check.o \
          grid_file_mod.o shape_functions_mod.o curvilinear_nodes_mod.o  \
          area_qpts_mod.o edge_qpts_mod.o \
          transformation.o vandermonde.o edge_connectivity_mod.o quit.o version.o

obj := $(patsubst %.o, $(ODIR)%.o,$(objects))

########################################################################
#                          Executable Target                           #
########################################################################
.PHONY :  version $(ODIR) bathy 

bathy:  version $(ODIR) $(obj)
	@echo "\n\n"
	@echo "Git SHA: $(GIT_SHA)"
	@echo "Modified: $(GIT_MOD) \n"
	$(FC) $(FFLAGS) -o bathy $(obj) $(LAPACK_LIB)
	
########################################################################
#                        Premilinary Targets                           #
########################################################################	

$(ODIR) :
	mkdir -p $@
	
version:
	@sed -i '25 c \      gitBranch = "$(GIT_BRANCH)" ' ../../src/version.F90
	@sed -i '26 c \      gitSHA = "$(GIT_SHA) $(GIT_MOD)" ' ../../src/version.F90
	@sed -i '27 c \      compiler_version = "$(FC_VERSION)" ' ../../src/version.F90
	@sed -i '28 c \      compiler_flags = "$(FFLAGS)" ' ../../src/version.F90
	@sed -i '29 c \      modified_files = "$(GIT_MFILES)" ' ../../src/version.F90
	@sed -i '30 c \      compile_date = "$(DATE)" ' ../../src/version.F90
	@sed -i '31 c \      host = "$(NAME)" ' ../../src/version.F90

	
########################################################################
#                            File Targets                              #
########################################################################	

$(ODIR)kdtree2.o: ../../src/kdtree2.F90
	$(FC) $(FFLAGS) -c $< -o $@
	mv kdtree2_*.mod  $(ODIR)
	
$(ODIR)version.o: ../../src/version.F90
	$(FC) $(FFLAGS) -c $< -o $@
	mv version.mod $(ODIR)

$(ODIR)system.o: ../../src/system.F90
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)lapack_interfaces.o : ../../src/lapack_interfaces.F90
	$(FC) $(FFLAGS) -c $< -o $@
	mv lapack_interfaces.mod $(ODIR)

$(ODIR)globals.o: globals.f90                                   $(ODIR)kdtree2.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv globals.mod $(ODIR)
	
$(ODIR)grid_file_mod.o : ../../src/grid_file_mod.F90            $(ODIR)globals.o $(ODIR)quit.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv grid_file_mod.mod $(ODIR)
	
$(ODIR)edge_connectivity_mod.o : ../../src/edge_connectivity_mod.F90 $(ODIR)quit.o $(ODIR)globals.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv edge_connectivity_mod.mod $(ODIR)	
	
$(ODIR)curvilinear_nodes_mod.o : ../../src/curvilinear_nodes_mod.F90  $(ODIR)globals.o $(ODIR)shape_functions_mod.o $(ODIR)transformation.o $(ODIR)basis.o $(ODIR)lapack_interfaces.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv curvilinear_nodes_mod.mod $(ODIR)
	
$(ODIR)shape_functions_mod.o : ../../src/shape_functions_mod.F90 $(ODIR)globals.o $(ODIR)basis.o $(ODIR)lapack_interfaces.o $(ODIR)vandermonde.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv shape_functions_mod.mod $(ODIR)

$(ODIR)transformation.o : ../../src/transformation.F90          $(ODIR)globals.o $(ODIR)basis.o $(ODIR)lapack_interfaces.o $(ODIR)vandermonde.o $(ODIR)shape_functions_mod.o $(ODIR)area_qpts_mod.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv transformation.mod $(ODIR)
	
$(ODIR)area_qpts_mod.o : ../../src/area_qpts_mod.F90            $(ODIR)globals.o $(ODIR)edge_qpts_mod.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv area_qpts_mod.mod $(ODIR)

$(ODIR)edge_qpts_mod.o : ../../src/edge_qpts_mod.F90            $(ODIR)globals.o 
	$(FC) $(FFLAGS) -c $< -o $@
	mv edge_qpts_mod.mod $(ODIR)		
	
$(ODIR)vandermonde.o : ../../src/vandermonde.F90                $(ODIR)globals.o $(ODIR)basis.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv vandermonde.mod $(ODIR)
	
$(ODIR)bathymetry_interp_mod.o : ../../src/bathymetry_interp_mod.F90 $(ODIR)globals.o $(ODIR)basis.o $(ODIR)transformation.o $(ODIR)shape_functions_mod.o $(ODIR)lapack_interfaces.o $(ODIR)vandermonde.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv bathymetry_interp_mod.mod $(ODIR)
	
$(ODIR)quit.o : ../../src/quit.F90                              $(ODIR)globals.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv quit.mod $(ODIR)	

$(ODIR)allocation.o: allocation.f90                             $(ODIR)globals.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv allocation.mod $(ODIR)

$(ODIR)basis.o: ../../src/basis.F90                             $(ODIR)globals.o $(ODIR)lapack_interfaces.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv basis.mod $(ODIR)

$(ODIR)evaluate.o: evaluate.F90                                 $(ODIR)globals.o $(ODIR)find_element_mod.o $(ODIR)shape_functions_mod.o $(ODIR)transformation.o $(ODIR)bathymetry_interp_mod.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv evaluate.mod $(ODIR)

$(ODIR)read_grid.o: read_grid.f90                               $(ODIR)globals.o  $(ODIR)allocation.o $(ODIR)grid_file_mod.o
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)find_element_mod.o: ../../src/find_element_mod.F90       $(ODIR)globals.o  $(ODIR)kdtree2.o  $(ODIR)shape_functions_mod.o $(ODIR)transformation.o $(ODIR)vandermonde.o  $(ODIR)sort_mod.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv find_element.mod $(ODIR)
	
$(ODIR)sort_mod.o : ../../src/sort_mod.F90                      $(ODIR)globals.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv sort_mod.mod $(ODIR)	

$(ODIR)read_input.o: read_input.f90                             $(ODIR)globals.o
	$(FC) $(FFLAGS) -c $< -o $@
	
$(ODIR)write_results.o: write_results.f90                       $(ODIR)globals.o  $(ODIR)version.o
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)connect.o: connect.f90                                   $(ODIR)globals.o  $(ODIR)edge_connectivity_mod.o
	$(FC) $(FFLAGS) -c $< -o $@
	
$(ODIR)curvilinear.o: curvilinear.f90                           $(ODIR)globals.o  $(ODIR)curvilinear_nodes_mod.o
	$(FC) $(FFLAGS) -c $< -o $@
	
$(ODIR)interp_check.o: interp_check.f90                         $(ODIR)globals.o 
	$(FC) $(FFLAGS) -c $< -o $@
	
$(ODIR)bathy_interp.o: bathy_interp.F90                         $(ODIR)globals.o  $(ODIR)allocation.o  $(ODIR)basis.o  $(ODIR)grid_file_mod.o  $(ODIR)evaluate.o  $(ODIR)find_element_mod.o $(ODIR)bathymetry_interp_mod.o $(ODIR)version.o
	$(FC) $(FFLAGS) -c $< -o $@
	


########################################################################
#                            Clean Targets                             #
########################################################################

.PHONY : clean
clean : 
	rm -r bathy $(ODIR)
