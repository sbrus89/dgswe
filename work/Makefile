SHELL = /bin/sh

VPATH := ../src 
ODIR = odir/

GIT_SHA := $(shell git rev-parse HEAD)    
GIT_MOD := $(shell git diff-index --quiet HEAD; echo $$?)
GIT_BRANCH := $(shell git rev-parse --abbrev-ref HEAD )

GIT_MOD := $(strip $(subst 1,+,$(GIT_MOD)))
GIT_MOD := $(strip $(subst 0,,$(GIT_MOD)))

########################################################################
#                  Compiler Flags (Build Specific)                     #
########################################################################

FC = ifort
FFLAGS1 = -c -O3 -132 -I$(ODIR) #-Dadcirc
#FFLAGS1 = -c -132 -C -g -traceback -I$(ODIR) -Dadcirc
FFLAGS2 = -O3 -o 

LIB = -llapack
#LIB = -Wl,--start-group   /afs/crc.nd.edu/x86_64_linux/intel/12.1/mkl/lib/intel64/libmkl_intel_lp64.a /afs/crc.nd.edu/x86_64_linux/intel/12.1/mkl/lib/intel64/libmkl_sequential.a  /afs/crc.nd.edu/x86_64_linux/intel/12.1/mkl/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm

########################################################################
#                               Objects                                #
########################################################################

objects = globals.o version.o allocation.o basis.o evaluate.o read_input.o read_grid.o read_solution.o kdtree2.o find_nesting.o area_qpts.o error.o

obj := $(patsubst %.o, $(ODIR)%.o,$(objects))

########################################################################
#                          Executable Target                           #
########################################################################
.PHONY : error $(ODIR) version

error: version $(ODIR) $(obj)
	@echo "\n\n"
	@echo "Git SHA: $(GIT_SHA)"
	@echo "Modified: $(GIT_MOD) \n"
	$(FC) $(FFLAGS2) error $(obj) $(LIB)

$(ODIR) :
	mkdir -p $@	

version:
	@sed -i '7 c \      PRINT*, "  Branch: $(GIT_BRANCH)" ' ../src/version.F90
	@sed -i '8 c \      PRINT*, "  SHA: $(GIT_SHA) $(GIT_MOD)" ' ../src/version.F90
	
########################################################################
#                            File Targets                              #
########################################################################	

$(ODIR)globals.o: globals.f90
	$(FC) $(FFLAGS1) $< -o $@
	mv globals.mod $(ODIR)
	
$(ODIR)allocation.o: allocation.f90 globals.f90
	$(FC) $(FFLAGS1) $< -o $@
	mv allocation.mod $(ODIR)
	
$(ODIR)basis.o: basis.f90 globals.f90
	$(FC) $(FFLAGS1) $< -o $@
	mv basis.mod $(ODIR)
	
$(ODIR)evaluate.o: evaluate.F90 basis.f90 globals.f90
	$(FC) $(FFLAGS1) $< -o $@
	mv evaluate.mod $(ODIR)

$(ODIR)read_input.o: read_input.f90 globals.f90
	$(FC) $(FFLAGS1) $< -o $@

$(ODIR)read_grid.o: read_grid.f90 globals.f90
	$(FC) $(FFLAGS1) $< -o $@
	mv read_grid.mod $(ODIR)

$(ODIR)read_solution.o: read_solution.F90 globals.f90
	$(FC) $(FFLAGS1) $< -o $@

$(ODIR)kdtree2.o: kdtree2.F
	$(FC) $(FFLAGS1) $< -o $@
	mv kdtree2_*.mod $(ODIR)

$(ODIR)find_nesting.o: find_nesting.f90 globals.f90 kdtree2.F
	$(FC) $(FFLAGS1) $< -o $@

$(ODIR)area_qpts.o: area_qpts.f90
	$(FC) $(FFLAGS1) $< -o $@

$(ODIR)version.o: version.F90
	$(FC) $(FFLAGS1) $< -o $@	

$(ODIR)error.o: error.F90 globals.f90 allocation.f90 evaluate.F90 basis.f90 kdtree2.F version.F90
	$(FC) $(FFLAGS1) $< -o $@
	
########################################################################
#                            Clean Targets                             #
########################################################################

.PHONY : clean
clean : 
	rm -r error $(ODIR)
