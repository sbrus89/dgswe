FC = ifort

VPATH := ../src

LIB_PATH = ../../metis-4.0.3/
LIB = metis
ODIR = odir/
CFLAGS = -c -C -g -traceback -I$(ODIR)

objects = globals.o allocation.o read_input.o read_grid.o read_forcing.o edge_connect.o metis2.o decomp2.o write_files.o dgprep.o
obj = $(patsubst %.o, $(ODIR)%.o,$(objects))

dgprep: $(ODIR) $(obj)
	$(FC) -o dgprep $(obj)  -L$(LIB_PATH) -l$(LIB)
	cp dgprep ../../work

.PHONY : $(ODIR)
$(ODIR) : 
	mkdir -p $@

$(ODIR)globals.o : globals.F90
	$(FC) $(CFLAGS) $< -o $@
	mv globals.mod $(ODIR)
	
$(ODIR)allocation.o : allocation.f90
	$(FC) $(CFLAGS) $< -o $@
	mv allocation.mod $(ODIR)	

$(ODIR)read_input.o : read_input.F90 globals.F90
	$(FC) $(CFLAGS) $< -o $@

$(ODIR)read_grid.o : read_grid.F90 globals.F90
	$(FC) $(CFLAGS) $< -o $@

$(ODIR)read_forcing.o : read_forcing.F90 globals.F90
	$(FC) $(CFLAGS) $< -o $@

$(ODIR)edge_connect.o : edge_connect.F90 globals.F90
	$(FC) $(CFLAGS) $< -o $@

$(ODIR)metis2.o : metis2.F90 globals.F90
	$(FC) $(CFLAGS) $< -o $@

$(ODIR)decomp2.o : decomp2.F90 globals.F90
	$(FC) $(CFLAGS) $< -o $@

$(ODIR)write_files.o : write_files.F90 globals.F90
	$(FC) $(CFLAGS) $< -o $@

$(ODIR)dgprep.o : dgprep.F90 globals.F90
	$(FC) $(CFLAGS) $< -o $@

.PHONY : metis
metis : 
	$(MAKE) -C ../metis-4.0.3/

.PHONY : cleanmetis
cleanmetis : 
	$(MAKE) -C ../metis-4.0.3/ clean

.PHONY : clean
clean : 
	rm dgprep odir/*
