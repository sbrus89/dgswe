SHELL = /bin/sh

VPATH := ../src
ODIR = odir/

########################################################################
#                          Compiler Flags                              #
########################################################################

FC = ifort
CFLAGS = -c -g -traceback -I$(ODIR)

########################################################################
#                           Library Links                              #
########################################################################

LIB_PATH = ../../metis-4.0.3/
LIB = metis

########################################################################
#                               Objects                                #
########################################################################

objects = globals.o read_input_vars.o read_input.o allocation.o read_grid.o read_forcing.o  \
          edge_connect.o metis2.o decomp2.o write_files.o dgprep.o
obj = $(patsubst %.o, $(ODIR)%.o,$(objects))

########################################################################
#                          Executable Targets                          #
########################################################################

.PHONY. : $(ODIR) dgprep clean metis cleanmetis

dgprep: $(ODIR) $(obj)
	$(FC) -o dgprep $(obj)  -L$(LIB_PATH) -l$(LIB)
	cp dgprep ../../work
metis : 
	$(MAKE) -C ../../metis-4.0.3/

########################################################################
#                            Clean Targets                             #
########################################################################

$(ODIR) : 
	mkdir -p $@

cleanmetis : 
	$(MAKE) -C ../../metis-4.0.3/ clean

clean : 
	rm dgprep odir/*

########################################################################
#                            File Targets                              #
########################################################################

# $< = first prerequisite/depenency
# $@ = name of the target		

$(ODIR)globals.o : globals.F90
	$(FC) $(CFLAGS) $< -o $@
	mv globals.mod $(ODIR)
	
$(ODIR)allocation.o : allocation.f90                                   $(ODIR)globals.o $(ODIR)read_input.o
	$(FC) $(CFLAGS) $< -o $@
	mv allocation.mod $(ODIR)

$(ODIR)read_input_vars.o : read_input_vars.F90 
	$(FC) $(CFLAGS) $< -o $@
	mv messenger2.mod $(ODIR)
	
# This is shared with dgswe 
$(ODIR)read_input.o : ../../src/read_input.F90                         $(ODIR)globals.o $(ODIR)read_input_vars.o         
	$(FC) $(CFLAGS) $< -o $@
	mv read_dginp.mod $(ODIR)

# This is shared with dgswe 
$(ODIR)read_input.o : ../../src/read_input.F90                         $(ODIR)globals.o $(ODIR)read_input_vars.o         
	$(FC) $(CFLAGS) $< -o $@
	mv read_dginp.mod $(ODIR)

$(ODIR)read_grid.o : read_grid.F90                                     $(ODIR)globals.o $(ODIR)allocation.o
	$(FC) $(CFLAGS) $< -o $@

$(ODIR)read_forcing.o : read_forcing.F90                               $(ODIR)globals.o $(ODIR)allocation.o $(ODIR)read_input.o
	$(FC) $(CFLAGS) $< -o $@

$(ODIR)edge_connect.o : edge_connect.F90                               $(ODIR)globals.o $(ODIR)allocation.o
	$(FC) $(CFLAGS) $< -o $@

$(ODIR)metis2.o : metis2.F90                                           $(ODIR)globals.o
	$(FC) $(CFLAGS) $< -o $@

$(ODIR)decomp2.o : decomp2.F90                                         $(ODIR)globals.o
	$(FC) $(CFLAGS) $< -o $@

$(ODIR)write_files.o : write_files.F90                                 $(ODIR)globals.o $(ODIR)read_input.o
	$(FC) $(CFLAGS) $< -o $@

$(ODIR)dgprep.o : dgprep.F90                                           $(ODIR)globals.o $(ODIR)allocation.o $(ODIR)read_input.o
	$(FC) $(CFLAGS) $< -o $@




