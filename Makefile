CC = ifort
CFLAGS1 = -c -132 -heap-arrays -I$(ODIR) #-C -traceback -g 
CFLAGS2 = -O3 -o -I$(ODIR) 
LIB = -llapack
ODIR = odir/

objects = kdtree2.o globals.o allocation.o basis.o find_element.o write_results.o connect.o evaluate.o read_input.o read_grid.o rimls.o
obj := $(patsubst %.o, $(ODIR)%.o,$(objects))

main: $(ODIR) $(obj)
	$(CC) $(CFLAGS2) rimls $(obj) $(LIB)

$(ODIR)kdtree2.o: kdtree2.F
	$(CC) $(CFLAGS1) $< -o $@
	mv kdtree2_precision_module.mod $(ODIR)
	mv kdtree2_module.mod $(ODIR)
	mv kdtree2_priority_queue_module.mod $(ODIR)

$(ODIR)globals.o: globals.f90 kdtree2.F
	$(CC) $(CFLAGS1) $< -o $@
	mv globals.mod $(ODIR)	

$(ODIR)allocation.o: allocation.f90 globals.f90
	$(CC) $(CFLAGS1) $< -o $@
	mv allocation.mod $(ODIR)

$(ODIR)basis.o: basis.f90 globals.f90
	$(CC) $(CFLAGS1) $< -o $@
	mv basis.mod $(ODIR)	
	
$(ODIR)find_element.o: find_element.f90 globals.f90 basis.f90 kdtree2.F
	$(CC) $(CFLAGS1) $< -o $@
	mv find_element.mod $(ODIR) 

$(ODIR)write_results.o: write_results.f90 globals.f90
	$(CC) $(CFLAGS1) $< -o $@
	mv write_results.mod $(ODIR)

$(ODIR)connect.o: connect.f90 globals.f90
	$(CC) $(CFLAGS1) $< -o $@

$(ODIR)evaluate.o: evaluate.f90 basis.f90 globals.f90 kdtree2.F
	$(CC) $(CFLAGS1) $< -o $@
	mv evaluate.mod $(ODIR)
	
$(ODIR)read_input.o: read_input.f90 globals.f90
	$(CC) $(CFLAGS1) $< -o $@

$(ODIR)read_grid.o: read_grid.f90 globals.f90
	$(CC) $(CFLAGS1) $< -o $@

$(ODIR)rimls.o: rimls.f90 globals.f90 evaluate.f90 allocation.f90 write_results.f90 basis.f90 kdtree2.F
	$(CC) $(CFLAGS1) $< -o $@

.PHONY : clean
clean : 
	rm -r rimls odir	
	
.PHONY : $(ODIR)
$(ODIR) :
	mkdir -p $@	
