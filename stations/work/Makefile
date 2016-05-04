FC = ifort
FFLAGS = -O3 -132 -I$(ODIR) -C -g -traceback
LIB = -llapack

VPATH := ../src/
ODIR = odir/

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

objects = globals.o quit.o lapack_interfaces.o basis.o read_input.o find_stations.o create_stations.o\
          messenger2.o grid_file_mod.o spherical_mod.o read_grid.o find_element_mod.o\
          edge_connectivity_mod.o read_stations.o kdtree2.o version.o stations.o \
          bathymetry_interp_mod.o transformation.o shape_functions_mod.o vandermonde.o

obj := $(patsubst %.o, $(ODIR)%.o,$(objects))

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

stations: version $(ODIR) $(obj)
	$(FC) $(CFLAGS2) -o stations $(obj) $(LIB)
	
$(ODIR)lapack_interfaces.o : ../../src/lapack_interfaces.F90
	$(FC) $(FFLAGS) -c $< -o $@
	mv lapack_interfaces.mod $(ODIR)
	
$(ODIR)version.o : ../../src/version.F90                                $(ODIR)globals.o
	$(FC) $(FFLAGS) -c $< -o $@	

$(ODIR)globals.o: ../../src/globals.F90
	$(FC) $(FFLAGS) -c  $< -o $@
	mv globals.mod $(ODIR) 
	
$(ODIR)basis.o: ../../src/basis.F90                                      $(ODIR)globals.o $(ODIR)lapack_interfaces.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv basis.mod $(ODIR)
	
$(ODIR)quit.o : ../../src/quit.F90                                       $(ODIR)globals.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv quit.mod $(ODIR)	
	
$(ODIR)messenger2.o : ../../src/messenger2.F90                           $(ODIR)globals.o $(ODIR)quit.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv messenger2.mod $(ODIR)	
	
$(ODIR)read_input.o : ../../src/read_input.F90                           $(ODIR)globals.o $(ODIR)messenger2.o $(ODIR)quit.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv read_dginp.mod $(ODIR)	

$(ODIR)read_grid.o : ../../src/read_grid.F90                             $(ODIR)globals.o $(ODIR)messenger2.o $(ODIR)read_input.o $(ODIR)grid_file_mod.o $(ODIR)spherical_mod.o
	$(FC) $(FFLAGS) -c $< -o $@
	
$(ODIR)spherical_mod.o : ../../src/spherical_mod.F90                     $(ODIR)globals.o 
	$(FC) $(FFLAGS) -c $< -o $@
	mv spherical_mod.mod $(ODIR)	
	
$(ODIR)grid_file_mod.o : ../../src/grid_file_mod.F90                     $(ODIR)globals.o $(ODIR)quit.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv grid_file_mod.mod $(ODIR)		
	
$(ODIR)edge_connectivity_mod.o : ../../src/edge_connectivity_mod.F90     $(ODIR)quit.o $(ODIR)globals.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv edge_connectivity_mod.mod $(ODIR)	
	
$(ODIR)find_element_mod.o: ../../src/find_element_mod.F90                $(ODIR)kdtree2.o  $(ODIR)globals.o $(ODIR)shape_functions_mod.o $(ODIR)transformation.o $(ODIR)vandermonde.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv find_element.mod $(ODIR) 	
	
$(ODIR)vandermonde.o : ../../src/vandermonde.F90                         $(ODIR)globals.o $(ODIR)basis.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv vandermonde.mod $(ODIR)	
	
$(ODIR)shape_functions_mod.o : ../../src/shape_functions_mod.F90         $(ODIR)globals.o $(ODIR)basis.o $(ODIR)lapack_interfaces.o $(ODIR)vandermonde.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv shape_functions_mod.mod $(ODIR)	
	
$(ODIR)bathymetry_interp_mod.o : ../../src/bathymetry_interp_mod.F90     $(ODIR)globals.o $(ODIR)basis.o $(ODIR)transformation.o $(ODIR)shape_functions_mod.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv bathymetry_interp_mod.mod $(ODIR)			

$(ODIR)transformation.o : ../../src/transformation.F90                   $(ODIR)globals.o $(ODIR)basis.o $(ODIR)lapack_interfaces.o $(ODIR)vandermonde.o $(ODIR)shape_functions_mod.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv transformation.mod $(ODIR)		

$(ODIR)read_stations.o: ../../src/read_stations.F90                      $(ODIR)globals.o $(ODIR)read_input.o
	$(FC) $(FFLAGS) -c $< -o $@
	
$(ODIR)find_stations.o: ../../src/find_stations.F90                      $(ODIR)globals.o  $(ODIR)basis.o  $(ODIR)read_input.o $(ODIR)bathymetry_interp_mod.o $(ODIR)find_element_mod.o $(ODIR)shape_functions_mod.o
	$(FC) $(FFLAGS) -c $< -o $@	
	
$(ODIR)create_stations.o: create_stations.f90                            $(ODIR)globals.o $(ODIR)read_input.o
	$(FC) $(FFLAGS) -c $< -o $@	

$(ODIR)kdtree2.o: ../../src/kdtree2.F
	$(FC) $(FFLAGS) -c $< -o $@
	mv kdtree2_precision_module.mod $(ODIR)
	mv kdtree2_module.mod $(ODIR)
	mv kdtree2_priority_queue_module.mod $(ODIR)

$(ODIR)stations.o: stations.f90                                          $(ODIR)globals.o   $(ODIR)read_input.o $(ODIR)edge_connectivity_mod.o $(ODIR)messenger2.o $(ODIR)create_stations.o
	$(FC) $(FFLAGS) -c $< -o $@

.PHONY : clean
clean : 
	rm stations odir/*.o odir/*.mod