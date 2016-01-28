SHELL = /bin/sh

VPATH := ../src/
ODIR  = odir/

FC = ifort
FFLAGS = -O2 -C -g -traceback
LIB = -llapack


spline_objects = globals.o read_grid.o read_input.o connect.o write_grid.o spline.o
spline_obj := $(patsubst %.o, $(ODIR)%.o,$(spline_objects)) 
FFLAGS += -I$(ODIR)




.PHONY : $(ODIR) clean spline

spline: $(ODIR) $(spline_obj)
	$(FC) $(FFLAGS) -o spline $(spline_obj)

$(ODIR) :
	mkdir -p $@





# $< = first prerequisite/depenency
# $@ = name of the target	


$(ODIR)globals.o: globals.f90
	$(FC) $(FFLAGS) -c $< -o $@
	mv globals.mod $(ODIR)

$(ODIR)spline.o: spline.f90                                  $(ODIR)globals.o
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)read_grid.o: read_grid.f90                            $(ODIR)globals.o
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)read_input.o: read_input.f90                          $(ODIR)globals.o
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)connect.o: connect.f90                                $(ODIR)globals.o
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)write_grid.o: write_grid.f90                          $(ODIR)globals.o
	$(FC) $(FFLAGS) -c $< -o $@





clean : 
	rm -rf spline odir