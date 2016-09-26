SHELL = /bin/sh

FC = ifort
FFLAGS = -I$(ODIR) -traceback -g -C

ODIR = odir/
FULLNAME := $(shell hostname)
FC_VERSION := $(shell $(FC) -v 2>&1 | head -2)
DATE := $(shell date)

GIT_SHA := $(shell git rev-parse HEAD)
GIT_BRANCH := $(shell git rev-parse --abbrev-ref HEAD )
GIT_MOD := $(shell git diff-index --quiet HEAD; echo $$?)
GIT_MOD := $(strip $(subst 1,+,$(GIT_MOD)))
GIT_MOD := $(strip $(subst 0,,$(GIT_MOD)))
GIT_MFILES := $(strip $(shell git ls-files ../* -m))

########################################################################
#                               Objects                                #
########################################################################

convert_objects =  globals.o  quit.o \
           spherical_mod.o \
           grid_file_mod.o \
           convert.o 
          
convert_obj := $(patsubst %.o, $(ODIR)%.o,$(convert_objects))

adjust_objects =  globals.o  quit.o \
           grid_file_mod.o \
           adjust.o 
          
adjust_obj := $(patsubst %.o, $(ODIR)%.o,$(adjust_objects))

convert_bou_objects =  globals.o  quit.o \
           grid_file_mod.o \
           convert_open_bou.o 
          
convert_bou_obj := $(patsubst %.o, $(ODIR)%.o,$(convert_bou_objects))

refine_objects = globals.o quit.o \
                 grid_file_mod.o edge_connectivity_mod.o \
                 version.o system.o refine.o
                 
refine_obj := $(patsubst %.o, $(ODIR)%.o,$(refine_objects))   

copy_objects = globals.o quit.o kdtree2.o\
               grid_file_mod.o copy_nodes.o
               
copy_obj := $(patsubst %.o, $(ODIR)%.o,$(copy_objects))                

########################################################################
#                          Executable Targets                          #
########################################################################
.PHONY : convert adjust

convert:  $(ODIR) $(convert_obj)
	@echo "\n\n"
	$(FC) $(FFLAGS) -o convert $(convert_obj) 
	
adjust:  $(ODIR) $(adjust_obj)
	@echo "\n\n"
	$(FC) $(FFLAGS) -o adjust $(adjust_obj)
	
convert_bou: $(ODIR) $(convert_bou_obj)
	@echo "\n\n"
	$(FC) $(FFLAGS) -o convert_bou $(convert_bou_obj)
	
refine: $(ODIR) $(refine_obj)	
	@echo "\n\n"
	$(FC) $(FFLAGS) -o refine $(refine_obj)
	
copy: $(ODIR) $(copy_obj)	
	@echo "\n\n"
	$(FC) $(FFLAGS) -o copy $(copy_obj)	
        
########################################################################
#                      Premilinary Targets                             #
########################################################################        
.PHONY : $(ODIR) 

$(ODIR) :
	mkdir -p $@     
        
version :
	@sed -i '5 c \      CHARACTER(:), ALLOCATABLE, PARAMETER :: gitBranch = "$(GIT_BRANCH)" ' ../src/version.F90
	@sed -i '6 c \      CHARACTER(:), ALLOCATABLE, PARAMETER :: gitSHA = "$(GIT_SHA) $(GIT_MOD)" ' ../src/version.F90
	@sed -i '7 c \      CHARACTER(:), ALLOCATABLE, PARAMETER :: compiler_version = "$(FC_VERSION)" ' ../src/version.F90
	@sed -i '8 c \      CHARACTER(:), ALLOCATABLE, PARAMETER :: compiler_flags = "$(FFLAGS)" ' ../src/version.F90
	@sed -i '9 c \      CHARACTER(:), ALLOCATABLE, PARAMETER :: modified_files = "$(GIT_MFILES)" ' ../src/version.F90
	@sed -i '10 c \      CHARACTER(:), ALLOCATABLE, PARAMETER :: compile_date = "$(DATE)" ' ../src/version.F90
	@sed -i '11 c \      CHARACTER(:), ALLOCATABLE, PARAMETER :: host = "$(NAME)" ' ../src/version.F90	        
########################################################################
#                            Clean Targets                             #
########################################################################
.PHONY : clean

clean : 
	rm -r convert adjust odir refine copy      
        
########################################################################
#                            File Targets                              #
########################################################################        

$(ODIR)globals.o: ../src/globals.F90                                   
	$(FC) $(FFLAGS) -c $< -o $@
	mv globals.mod $(ODIR)         	  	
	
$(ODIR)version.o : ../src/version.F90                                
	$(FC) $(FFLAGS) -c $< -o $@
	mv version.mod $(ODIR)
	
$(ODIR)system.o : ../src/system.F90                                
	$(FC) $(FFLAGS) -c $< -o $@	
	
$(ODIR)kdtree2.o: ../src/kdtree2.F
	$(FC) $(FFLAGS) -132 -c $< -o $@
	mv kdtree2_precision_module.mod $(ODIR)
	mv kdtree2_module.mod $(ODIR)
	mv kdtree2_priority_queue_module.mod $(ODIR)	

$(ODIR)grid_file_mod.o : ../src/grid_file_mod.F90            $(ODIR)globals.o $(ODIR)quit.o $(ODIR)spherical_mod.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv grid_file_mod.mod $(ODIR)        
        
$(ODIR)quit.o : ../src/quit.F90                              $(ODIR)globals.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv quit.mod $(ODIR)     
        
$(ODIR)spherical_mod.o : ../src/spherical_mod.F90            $(ODIR)globals.o 
	$(FC) $(FFLAGS) -c $< -o $@
	mv spherical_mod.mod $(ODIR)            
	
$(ODIR)edge_connectivity_mod.o : ../src/edge_connectivity_mod.F90 $(ODIR)globals.o $(ODIR)quit.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv edge_connectivity_mod.mod $(ODIR)

$(ODIR)convert.o: convert2cartesian.F90                      $(ODIR)globals.o $(ODIR)grid_file_mod.o $(ODIR)spherical_mod.o
	$(FC) $(FFLAGS) -c $< -o $@
        
$(ODIR)adjust.o: adjust_channel_boundaries.f90               $(ODIR)globals.o $(ODIR)grid_file_mod.o 
	$(FC) $(FFLAGS) -c $< -o $@       
	
$(ODIR)convert_open_bou.o: convert_open_boundaries.F90       $(ODIR)globals.o $(ODIR)grid_file_mod.o
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)refine.o: refine_grid.F90                             $(ODIR)globals.o $(ODIR)grid_file_mod.o $(ODIR)edge_connectivity_mod.o $(ODIR)version.o
	$(FC) $(FFLAGS) -c $< -o $@
	
$(ODIR)copy_nodes.o: copy_nodes.F90                          $(ODIR)globals.o $(ODIR)grid_file_mod.o $(ODIR)kdtree2.o
	$(FC) $(FFLAGS) -c $< -o $@