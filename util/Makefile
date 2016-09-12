FC = ifort
FFLAGS = -I$(ODIR) #-traceback -g -C

ODIR = odir/


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
        
########################################################################
#                      Premilinary Targets                             #
########################################################################        
.PHONY : $(ODIR) 

$(ODIR) :
	mkdir -p $@     
        
        
########################################################################
#                            Clean Targets                             #
########################################################################
.PHONY : clean

clean : 
	rm -r convert adjust odir        
        
########################################################################
#                            File Targets                              #
########################################################################        

$(ODIR)globals.o: ../src/globals.F90                                   
	$(FC) $(FFLAGS) -c $< -o $@
	mv globals.mod $(ODIR)          

$(ODIR)grid_file_mod.o : ../src/grid_file_mod.F90            $(ODIR)globals.o $(ODIR)quit.o $(ODIR)spherical_mod.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv grid_file_mod.mod $(ODIR)        
        
$(ODIR)quit.o : ../src/quit.F90                              $(ODIR)globals.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv quit.mod $(ODIR)     
        
$(ODIR)spherical_mod.o : ../src/spherical_mod.F90            $(ODIR)globals.o 
	$(FC) $(FFLAGS) -c $< -o $@
	mv spherical_mod.mod $(ODIR)            

$(ODIR)convert.o: convert2cartesian.F90                      $(ODIR)globals.o $(ODIR)grid_file_mod.o $(ODIR)spherical_mod.o
	$(FC) $(FFLAGS) -c $< -o $@
        
$(ODIR)adjust.o: adjust_channel_boundaries.f90               $(ODIR)globals.o $(ODIR)grid_file_mod.o 
	$(FC) $(FFLAGS) -c $< -o $@       
	
$(ODIR)convert_open_bou.o: convert_open_boundaries.F90            $(ODIR)globals.o $(ODIR)grid_file_mod.o
	$(FC) $(FFLAGS) -c $< -o $@	

