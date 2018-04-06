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

########################################################################
#                               Objects                                #
########################################################################

convert_objects =  globals.o  quit.o \
           spherical_mod.o \
           grid_file_mod.o \
           convert.o 
          
convert_obj := $(patsubst %.o, $(ODIR)%.o,$(convert_objects))

########################################################################

adjust_objects =  globals.o  quit.o \
           grid_file_mod.o \
           adjust.o 
          
adjust_obj := $(patsubst %.o, $(ODIR)%.o,$(adjust_objects))

########################################################################

equad_objects =  globals.o  quit.o \
           grid_file_mod.o basis.o\
           elongated_quad_channel.o 
          
equad_obj := $(patsubst %.o, $(ODIR)%.o,$(equad_objects))

########################################################################

convert_bou_objects =  globals.o  quit.o \
           grid_file_mod.o \
           convert_open_bou.o 
          
convert_bou_obj := $(patsubst %.o, $(ODIR)%.o,$(convert_bou_objects))

########################################################################

refine_objects = globals.o quit.o \
                 grid_file_mod.o edge_connectivity_mod.o \
                 version.o system.o refine.o
                 
refine_obj := $(patsubst %.o, $(ODIR)%.o,$(refine_objects))  

########################################################################

copy_objects = globals.o quit.o kdtree2.o\
               grid_file_mod.o copy_nodes.o
               
copy_obj := $(patsubst %.o, $(ODIR)%.o,$(copy_objects))  

########################################################################

bou2sta_objects =  globals.o  quit.o \
           grid_file_mod.o \
           convert_bou2sta.o 
                   
bou2sta_obj := $(patsubst %.o, $(ODIR)%.o,$(bou2sta_objects))

########################################################################

consolidate_objects =  globals.o  quit.o \
           grid_file_mod.o edge_connectivity_mod.o \
           triangulation.o tripack.o consolidate.o \
           kdtree2.o
                   
consolidate_obj := $(patsubst %.o, $(ODIR)%.o,$(consolidate_objects))

########################################################################

cfl_objects =  cfl_checker.o globals.o  quit.o  grid_file_mod.o \
               read_input.o  spherical_mod.o  system.o allocation.o messenger2.o version.o
                   
cfl_obj := $(patsubst %.o, $(ODIR)%.o,$(cfl_objects))

########################################################################

remove_bou_objects = globals.o grid_file_mod.o quit.o spherical_mod.o remove_bou.o

remove_bou_obj := $(patsubst %.o, $(ODIR)%.o,$(remove_bou_objects))

########################################################################

write_bou_file_objects = globals.o grid_file_mod.o quit.o spherical_mod.o \
                         allocation.o messenger2.o read_input.o \
                         system.o version.o write_bou_file.o

write_bou_file_obj := $(patsubst %.o, $(ODIR)%.o,$(write_bou_file_objects))   

########################################################################

find_const_objects = globals.o grid_file_mod.o quit.o kdtree2.o find_const.o

find_const_obj := $(patsubst %.o, $(ODIR)%.o,$(find_const_objects))   

########################################################################

create_fort13_objects = globals.o nodal_attributes_mod.o grid_file_mod.o quit.o kdtree2.o create_fort13.o

create_fort13_obj := $(patsubst %.o, $(ODIR)%.o,$(create_fort13_objects))   

########################################################################

essential2natural_objects = globals.o grid_file_mod.o quit.o essential2natural.o

essential2natural_obj := $(patsubst %.o, $(ODIR)%.o,$(essential2natural_objects))  

########################################################################

stations_objects = globals.o kdtree2.o basis.o grid_file_mod.o quit.o sort_mod.o \
                   edge_connectivity_mod.o find_element_mod.o shape_functions_mod.o \
                   transformation.o vandermonde.o area_qpts_mod.o edge_qpts_mod.o stations.o  

stations_obj := $(patsubst %.o, $(ODIR)%.o,$(stations_objects))   


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
	
equad:  $(ODIR) $(equad_obj)
	@echo "\n\n"
	$(FC) $(FFLAGS) -o equad $(equad_obj)	$(LAPACK_LIB)
	
convert_bou: $(ODIR) $(convert_bou_obj)
	@echo "\n\n"
	$(FC) $(FFLAGS) -o convert_bou $(convert_bou_obj)
	
refine: $(ODIR) $(refine_obj)	
	@echo "\n\n"
	$(FC) $(FFLAGS) -o refine $(refine_obj)
	
copy: $(ODIR) $(copy_obj)	
	@echo "\n\n"
	$(FC) $(FFLAGS) -o copy $(copy_obj)	
	
bou2sta: $(ODIR) $(bou2sta_obj)	
	@echo "\n\n"
	$(FC) $(FFLAGS) -o bou2sta $(bou2sta_obj)
	
consolidate: $(ODIR) $(consolidate_obj)
	@echo "\n\n"
	$(FC) $(FFLAGS) -o consolidate $(consolidate_obj)
	
cfl_checker: $(ODIR) $(cfl_obj)
	@echo "\n\n"
	$(FC) $(FFLAGS) -o cfl_checker $(cfl_obj)
	
remove_bou: $(ODIR) $(remove_bou_obj)
	@echo "\n\n"
	$(FC) $(FFLAGS) -o remove_bou $(remove_bou_obj)  
	
write_bou: $(ODIR) $(write_bou_file_obj)
	@echo "\n\n"
	$(FC) $(FFLAGS) -o write_bou $(write_bou_file_obj)  	
	
find_const: $(ODIR) $(find_const_obj)
	@echo "\n\n"
	$(FC) $(FFLAGS) -o find_const $(find_const_obj)  		
	
create_fort13: $(ODIR) $(create_fort13_obj)
	@echo "\n\n"
	$(FC) $(FFLAGS) -o create_fort13 $(create_fort13_obj) 
	
essential2natural: $(ODIR) $(essential2natural_obj)
	@echo "\n\n"
	$(FC) $(FFLAGS) -o essential2natural $(essential2natural_obj) 
	
stations: $(ODIR) $(stations_obj)
	@echo "\n\n"
	$(FC) $(FFLAGS) -o stations $(stations_obj) $(LAPACK_LIB)
	
########################################################################
#                      Premilinary Targets                             #
########################################################################        
.PHONY : $(ODIR) 

$(ODIR) :
	mkdir -p $@     
        
version :
	@sed -i '25 c \      gitBranch = "$(GIT_BRANCH)" ' ../../src/version.F90
	@sed -i '26 c \      gitSHA = "$(GIT_SHA) $(GIT_MOD)" ' ../../src/version.F90
	@sed -i '27 c \      compiler_version = "$(FC_VERSION)" ' ../../src/version.F90
	@sed -i '28 c \      compiler_flags = "$(FFLAGS)" ' ../../src/version.F90
	@sed -i '29 c \      modified_files = "$(GIT_MFILES)" ' ../../src/version.F90
	@sed -i '30 c \      compile_date = "$(DATE)" ' ../../src/version.F90
	@sed -i '31 c \      host = "$(NAME)" ' ../../src/version.F90
     
########################################################################
#                            Clean Targets                             #
########################################################################
.PHONY : clean

clean : 
	rm -r convert adjust equad odir refine copy consolidate cfl_checker remove_bou write_bou find_const create_fort13
        
########################################################################
#                            File Targets                              #
########################################################################        

$(ODIR)globals.o: ../src/globals.F90                                   
	$(FC) $(FFLAGS) -c $< -o $@
	mv globals.mod $(ODIR)    
	
$(ODIR)version.o : ../src/version.F90                                
	$(FC) $(FFLAGS) -c $< -o $@
	mv version.mod $(ODIR)
	
$(ODIR)messenger2.o : ../src/messenger2.F90                           $(ODIR)globals.o $(ODIR)quit.o $(ODIR)version.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv messenger2.mod $(ODIR)	
	
$(ODIR)allocation.o : ../src/allocation.F90                           $(ODIR)globals.o $(ODIR)messenger2.o $(ODIR)read_input.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv allocation.mod $(ODIR)	
	
$(ODIR)system.o : ../src/system.F90                                
	$(FC) $(FFLAGS) -c $< -o $@	
	
$(ODIR)kdtree2.o: ../src/kdtree2.F90
	$(FC) $(FFLAGS) -132 -c $< -o $@
	mv kdtree2_precision_module.mod $(ODIR)
	mv kdtree2_module.mod $(ODIR)
	mv kdtree2_priority_queue_module.mod $(ODIR)	
	
$(ODIR)sort_mod.o : ../src/sort_mod.F90                      $(ODIR)globals.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv sort_mod.mod $(ODIR)	
	
$(ODIR)area_qpts_mod.o : ../src/area_qpts_mod.F90            $(ODIR)globals.o $(ODIR)edge_qpts_mod.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv area_qpts_mod.mod $(ODIR)
	
$(ODIR)edge_qpts_mod.o : ../src/edge_qpts_mod.F90            $(ODIR)globals.o 
	$(FC) $(FFLAGS) -c $< -o $@
	mv edge_qpts_mod.mod $(ODIR)		
	
$(ODIR)lapack_interfaces.o : ../src/lapack_interfaces.F90
	$(FC) $(FFLAGS) -c $< -o $@
	mv lapack_interfaces.mod $(ODIR)
	
$(ODIR)vandermonde.o : ../src/vandermonde.F90                $(ODIR)globals.o $(ODIR)basis.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv vandermonde.mod $(ODIR)	
	
$(ODIR)basis.o : ../src/basis.F90                            $(ODIR)globals.o $(ODIR)lapack_interfaces.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv basis.mod $(ODIR)	

$(ODIR)grid_file_mod.o : ../src/grid_file_mod.F90            $(ODIR)globals.o $(ODIR)quit.o $(ODIR)spherical_mod.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv grid_file_mod.mod $(ODIR)        
        
$(ODIR)nodal_attributes_mod.o : ../src/nodal_attributes.F90  $(ODIR)globals.o $(ODIR)kdtree2.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv nodal_attributes_mod.mod $(ODIR)	        
        
$(ODIR)quit.o : ../src/quit.F90                              $(ODIR)globals.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv quit.mod $(ODIR)     
        
$(ODIR)spherical_mod.o : ../src/spherical_mod.F90            $(ODIR)globals.o 
	$(FC) $(FFLAGS) -c $< -o $@
	mv spherical_mod.mod $(ODIR)            
	
$(ODIR)edge_connectivity_mod.o : ../src/edge_connectivity_mod.F90 $(ODIR)globals.o $(ODIR)quit.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv edge_connectivity_mod.mod $(ODIR)
	
$(ODIR)read_input.o : ../src/read_input.F90                  $(ODIR)globals.o $(ODIR)quit.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv read_dginp.mod $(ODIR)	
	
$(ODIR)triangulation.o : ../src/triangulation.F90            $(ODIR)globals.o         
	$(FC) $(FFLAGS) -c $< -o $@
	mv triangulation.mod $(ODIR)
	
$(ODIR)shape_functions_mod.o : ../src/shape_functions_mod.F90 $(ODIR)globals.o $(ODIR)basis.o $(ODIR)lapack_interfaces.o $(ODIR)vandermonde.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv shape_functions_mod.mod $(ODIR)	
	
$(ODIR)find_element_mod.o: ../src/find_element_mod.F90       $(ODIR)kdtree2.o  $(ODIR)globals.o $(ODIR)shape_functions_mod.o $(ODIR)transformation.o $(ODIR)vandermonde.o $(ODIR)sort_mod.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv find_element.mod $(ODIR) 	
	
$(ODIR)transformation.o : ../src/transformation.F90          $(ODIR)basis.o $(ODIR)vandermonde.o $(ODIR)lapack_interfaces.o $(ODIR)shape_functions_mod.o $(ODIR)area_qpts_mod.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv transformation.mod $(ODIR)	
	
$(ODIR)tripack.o : ../tripack/Src/Fortran77/Sp/src.f90                                
	$(FC) $(FFLAGS) -c $< -o $@	

$(ODIR)convert.o: convert2cartesian.F90                      $(ODIR)globals.o $(ODIR)grid_file_mod.o $(ODIR)spherical_mod.o
	$(FC) $(FFLAGS) -c $< -o $@
        
$(ODIR)adjust.o: adjust_channel_boundaries.f90               $(ODIR)globals.o $(ODIR)grid_file_mod.o 
	$(FC) $(FFLAGS) -c $< -o $@       
	
$(ODIR)elongated_quad_channel.o: elongated_quad_channel.F90  $(ODIR)globals.o $(ODIR)grid_file_mod.o $(ODIR)basis.o
	$(FC) $(FFLAGS) -c $< -o $@       	
	
$(ODIR)convert_open_bou.o: convert_open_boundaries.F90       $(ODIR)globals.o $(ODIR)grid_file_mod.o
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)refine.o: refine_grid.F90                             $(ODIR)globals.o $(ODIR)grid_file_mod.o $(ODIR)edge_connectivity_mod.o $(ODIR)version.o
	$(FC) $(FFLAGS) -c $< -o $@
	
$(ODIR)copy_nodes.o: copy_nodes.F90                          $(ODIR)globals.o $(ODIR)grid_file_mod.o $(ODIR)kdtree2.o
	$(FC) $(FFLAGS) -c $< -o $@
	
$(ODIR)convert_bou2sta.o: convert_bou2sta.F90                $(ODIR)globals.o $(ODIR)grid_file_mod.o
	$(FC) $(FFLAGS) -c $< -o $@	
	
$(ODIR)consolidate.o: consolidate.F90                        $(ODIR)globals.o $(ODIR)grid_file_mod.o $(ODIR)edge_connectivity_mod.o $(ODIR)triangulation.o $(ODIR)kdtree2.o
	$(FC) $(FFLAGS) -c $< -o $@
	
$(ODIR)cfl_checker.o: cfl_checker.f90                        $(ODIR)globals.o $(ODIR)grid_file_mod.o $(ODIR)read_input.o $(ODIR)allocation.o
	$(FC) $(FFLAGS) -c $< -o $@	
	
$(ODIR)remove_bou.o: remove_bou_type.F90                     $(ODIR)globals.o $(ODIR)grid_file_mod.o 
	$(FC) $(FFLAGS) -c $< -o $@
	
$(ODIR)write_bou_file.o: write_bou_file.F90                  $(ODIR)globals.o $(ODIR)grid_file_mod.o $(ODIR)allocation.o
	$(FC) $(FFLAGS) -c $< -o $@	
	
$(ODIR)find_const.o: find_constraints.F90                    $(ODIR)globals.o $(ODIR)grid_file_mod.o 
	$(FC) $(FFLAGS) -c $< -o $@
	
$(ODIR)create_fort13.o: create_fort13.F90                    $(ODIR)nodal_attributes_mod.o $(ODIR)grid_file_mod.o $(ODIR)kdtree2.o
	$(FC) $(FFLAGS) -c $< -o $@	
	
$(ODIR)essential2natural.o: essential2natural.F90            $(ODIR)globals.o $(ODIR)grid_file_mod.o 
	$(FC) $(FFLAGS) -c $< -o $@
	
$(ODIR)stations.o: stations.F90                              $(ODIR)globals.o $(ODIR)grid_file_mod.o $(ODIR)edge_connectivity_mod.o $(ODIR)shape_functions_mod.o
	$(FC) $(FFLAGS) -c $< -o $@	