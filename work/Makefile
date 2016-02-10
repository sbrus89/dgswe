SHELL = /bin/sh

VPATH := ../src/
ODIR  = odir/

FC = ifort
FFLAGS = -O2 -C -g -traceback
LIB = -llapack


spline_objects = globals.o kdtree2.o lapack_interfaces.o basis.o allocation.o  \
                 evaluate.o calc_spline.o check.o find_element.o \
                 read_grid.o read_input.o connect.o write_output.o spline.o
spline_obj := $(patsubst %.o, $(ODIR)%.o,$(spline_objects)) 
FFLAGS += -I$(ODIR)




.PHONY : $(ODIR) clean spline

spline: $(ODIR) $(spline_obj)
	$(FC) $(FFLAGS) -o spline $(spline_obj) $(LIB)

$(ODIR) :
	mkdir -p $@





# $< = first prerequisite/depenency
# $@ = name of the target	

$(ODIR)kdtree2.o : kdtree2.F
	$(FC) $(FFLAGS) -132 -c $< -o $@
	mv kdtree2_precision_module.mod $(ODIR)
	mv kdtree2_module.mod $(ODIR)
	mv kdtree2_priority_queue_module.mod $(ODIR)
	
$(ODIR)lapack_interfaces.o : lapack_interfaces.f90
	$(FC) $(FFLAGS) -c $< -o $@
	mv lapack_interfaces.mod $(ODIR)		

$(ODIR)globals.o : globals.f90                                $(ODIR)kdtree2.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv globals.mod $(ODIR)

$(ODIR)allocation.o : allocation.f90                          $(ODIR)globals.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv allocation.mod $(ODIR)
	
$(ODIR)basis.o : basis.f90                                    $(ODIR)globals.o  $(ODIR)lapack_interfaces.o 
	$(FC) $(FFLAGS) -c $< -o $@
	mv basis.mod $(ODIR)
	
$(ODIR)calc_spline.o : calc_spline.f90                        $(ODIR)globals.o  $(ODIR)kdtree2.o $(ODIR)lapack_interfaces.o 
	$(FC) $(FFLAGS) -c $< -o $@
	mv calc_spline.mod $(ODIR)
	
$(ODIR)check.o : check.f90                                    $(ODIR)globals.o $(ODIR)calc_spline.o  $(ODIR)find_element.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv check.mod $(ODIR)
	
$(ODIR)evaluate.o : evaluate.f90                              $(ODIR)globals.o $(ODIR)basis.o  $(ODIR)lapack_interfaces.o 
	$(FC) $(FFLAGS) -c $< -o $@
	mv evaluate.mod $(ODIR)
	
$(ODIR)find_element.o : find_element.f90                      $(ODIR)globals.o  $(ODIR)kdtree2.o  $(ODIR)basis.o $(ODIR)lapack_interfaces.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv find_element.mod $(ODIR)

$(ODIR)spline.o : spline.f90                                  $(ODIR)globals.o  $(ODIR)allocation.o  $(ODIR)calc_spline.o  $(ODIR)check.o $(ODIR)find_element.o
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)read_grid.o : read_grid.f90                            $(ODIR)globals.o  $(ODIR)allocation.o 
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)read_input.o : read_input.f90                          $(ODIR)globals.o
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)connect.o : connect.f90                                $(ODIR)globals.o
	$(FC) $(FFLAGS) -c $< -o $@

$(ODIR)write_output.o : write_output.f90                      $(ODIR)globals.o
	$(FC) $(FFLAGS) -c $< -o $@





clean : 
	rm -rf spline odir