SHELL = /bin/sh

VPATH := ../src/
ODIR  = odir/

FC = ifort
FFLAGS = -O2 -C -g -traceback
LIB = -llapack


spline_objects = globals.o kdtree2.o allocation.o calc_spline.o check.o \
                 read_grid.o read_input.o connect.o write_output.o spline.o
spline_obj := $(patsubst %.o, $(ODIR)%.o,$(spline_objects)) 
FFLAGS += -I$(ODIR)




.PHONY : $(ODIR) clean spline

spline: $(ODIR) $(spline_obj)
	$(FC) $(FFLAGS) -o spline $(spline_obj)

$(ODIR) :
	mkdir -p $@





# $< = first prerequisite/depenency
# $@ = name of the target	

$(ODIR)kdtree2.o: kdtree2.F
	$(FC) $(FFLAGS) -132 -c $< -o $@
	mv kdtree2_precision_module.mod $(ODIR)
	mv kdtree2_module.mod $(ODIR)
	mv kdtree2_priority_queue_module.mod $(ODIR)

$(ODIR)globals.o: globals.f90                                $(ODIR)kdtree2.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv globals.mod $(ODIR)

$(ODIR)allocation.o : allocation.f90                         $(ODIR)globals.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv allocation.mod $(ODIR)
	
$(ODIR)calc_spline.o: calc_spline.f90                        $(ODIR)globals.o  $(ODIR)kdtree2.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv calc_spline.mod $(ODIR)
	
$(ODIR)check.o: check.f90                                    $(ODIR)globals.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv check.mod $(ODIR)

$(ODIR)spline.o: spline.f90                                  $(ODIR)globals.o  $(ODIR)allocation.o  $(ODIR)calc_spline.o  $(ODIR)check.o
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)read_grid.o: read_grid.f90                            $(ODIR)globals.o  $(ODIR)allocation.o 
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)read_input.o: read_input.f90                          $(ODIR)globals.o
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)connect.o: connect.f90                                $(ODIR)globals.o
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)write_output.o: write_output.f90                      $(ODIR)globals.o
	$(FC) $(FFLAGS) -c $< -o $@





clean : 
	rm -rf spline odir