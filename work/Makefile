FC = ifort
FFLAGS = -132 -heap-arrays -I$(ODIR) -traceback -g -C

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

objects = kdtree2.o globals.o allocation.o basis.o find_element.o write_results.o connect.o evaluate.o read_input.o read_grid.o curvilinear.o version.o rimls.o
obj := $(patsubst %.o, $(ODIR)%.o,$(objects))

########################################################################
#                          Executable Targets                          #
########################################################################
.PHONY : rimls 

rimls: version $(ODIR) $(obj)
	$(FC) $(FFLAGS) -o rimls $(obj) $(LAPACK_LIB)
	
########################################################################
#                      Premilinary Targets                             #
########################################################################	
.PHONY : $(ODIR) version

$(ODIR) :
	mkdir -p $@	
	
version :
	@sed -i '10 c \      gitBranch = "$(GIT_BRANCH)" ' $(VPATH)version.f90
	@sed -i '11 c \      gitSHA = "$(GIT_SHA) $(GIT_MOD)" ' $(VPATH)version.f90
	@sed -i '12 c \      compiler_version = "$(FC_VERSION)" ' $(VPATH)version.f90
	@sed -i '13 c \      compiler_flags = "$(FFLAGS)" ' $(VPATH)version.f90
	@sed -i '14 c \      modified_files = "$(GIT_MFILES)" ' $(VPATH)version.f90
	@sed -i '15 c \      compile_date = "$(DATE)" ' $(VPATH)version.f90
	@sed -i '16 c \      host = "$(NAME)" ' $(VPATH)version.f90			
	
########################################################################
#                            Clean Targets                             #
########################################################################
.PHONY : clean

clean : 
	rm -r rimls odir	
	
########################################################################
#                            File Targets                              #
########################################################################	

$(ODIR)kdtree2.o: kdtree2.F
	$(FC) $(FFLAGS) -c $< -o $@
	mv kdtree2_precision_module.mod $(ODIR)
	mv kdtree2_module.mod $(ODIR)
	mv kdtree2_priority_queue_module.mod $(ODIR)

$(ODIR)globals.o: globals.f90                                   $(ODIR)kdtree2.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv globals.mod $(ODIR)	

$(ODIR)allocation.o: allocation.f90                             $(ODIR)globals.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv allocation.mod $(ODIR)

$(ODIR)basis.o: basis.f90                                       $(ODIR)globals.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv basis.mod $(ODIR)	
	
$(ODIR)find_element.o: find_element.f90                         $(ODIR)kdtree2.o  $(ODIR)globals.o  $(ODIR)basis.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv find_element.mod $(ODIR) 

$(ODIR)write_results.o: write_results.f90                       $(ODIR)globals.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv write_results.mod $(ODIR)

$(ODIR)connect.o: connect.f90                                   $(ODIR)globals.o
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)evaluate.o: evaluate.f90                                 $(ODIR)kdtree2.o  $(ODIR)globals.o  $(ODIR)basis.o  $(ODIR)find_element.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv evaluate.mod $(ODIR)
	
$(ODIR)read_input.o: read_input.f90                             $(ODIR)globals.o
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)read_grid.o: read_grid.f90                               $(ODIR)globals.o  $(ODIR)allocation.o
	$(FC) $(FFLAGS) -c $< -o $@
	
$(ODIR)curvilinear.o : curvilinear.f90                           $(ODIR)globals.o
	$(FC) $(FFLAGS) -c $< -o $@	
	
$(ODIR)version.o : version.f90                                  $(ODIR)globals.o
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)rimls.o: rimls.f90                                       $(ODIR)globals.o  $(ODIR)allocation.o  $(ODIR)basis.o  $(ODIR)write_results.o
	$(FC) $(FFLAGS) -c $< -o $@
	

