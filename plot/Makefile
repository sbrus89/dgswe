FC = ifort
FFLAGS = -I$(ODIR) #-traceback -g -C

ODIR = odir/


########################################################################
#                               Objects                                #
########################################################################


plot_objects =  plot_globals.o globals.o quit.o basis.o version.o read_input.o \
           read_plot_input.o grid_file_mod.o read_write_output.o system.o tripack.o \
           edge_connectivity_mod.o curvilinear_nodes_mod.o bathymetry_interp_mod.o \
           vandermonde.o transformation.o shape_functions_mod.o triangulation.o \
           plot_mod.o evaluate_mod.o axes_mod.o labels_mod.o plot.o 
          
plot_obj := $(patsubst %.o, $(ODIR)%.o,$(plot_objects))


########################################################################
#                          Executable Targets                          #
########################################################################
.PHONY :   plot
	
plot:  $(ODIR) $(plot_obj)
	@echo "\n\n"
	$(FC) $(FFLAGS) -o plot $(plot_obj) -llapack
        
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
	rm -r  plot odir        
        
########################################################################
#                            File Targets                              #
########################################################################        

$(ODIR)plot_globals.o: plot_globals.F90                        $(ODIR)globals.o $(ODIR)read_input.o           
	$(FC) $(FFLAGS) -c $< -o $@
	mv plot_globals.mod $(ODIR)    

$(ODIR)globals.o: ../src/globals.F90                                   
	$(FC) $(FFLAGS) -c $< -o $@
	mv globals.mod $(ODIR)          
	
$(ODIR)lapack_interfaces.o : ../src/lapack_interfaces.F90
	$(FC) $(FFLAGS) -c $< -o $@
	mv lapack_interfaces.mod $(ODIR)	

$(ODIR)grid_file_mod.o : ../src/grid_file_mod.F90            $(ODIR)globals.o $(ODIR)quit.o $(ODIR)spherical_mod.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv grid_file_mod.mod $(ODIR)        
        
$(ODIR)quit.o : ../src/quit.F90                              $(ODIR)globals.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv quit.mod $(ODIR)     
        
$(ODIR)spherical_mod.o : ../src/spherical_mod.F90            $(ODIR)globals.o 
	$(FC) $(FFLAGS) -c $< -o $@
	mv spherical_mod.mod $(ODIR)    
	
$(ODIR)basis.o : ../src/basis.F90                            $(ODIR)globals.o $(ODIR)lapack_interfaces.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv basis.mod $(ODIR)	
	
$(ODIR)read_write_output.o : ../src/read_write_output.F90    $(ODIR)globals.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv read_write_output.mod $(ODIR)	
	
$(ODIR)read_input.o : ../src/read_input.F90                  $(ODIR)globals.o $(ODIR)quit.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv read_dginp.mod $(ODIR)	
	
$(ODIR)version.o : ../src/version.F90                                
	$(FC) $(FFLAGS) -c $< -o $@
	mv version.mod $(ODIR)
	
$(ODIR)triangulation.o : ../src/triangulation.F90                                
	$(FC) $(FFLAGS) -c $< -o $@
	mv triangulation.mod $(ODIR)
	
$(ODIR)tripack.o : ../tripack/Src/Fortran77/Sp/src.f90                                
	$(FC) $(FFLAGS) -c $< -o $@
	
$(ODIR)edge_connectivity_mod.o : ../src/edge_connectivity_mod.F90 $(ODIR)globals.o $(ODIR)quit.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv edge_connectivity_mod.mod $(ODIR)
	
$(ODIR)curvilinear_nodes_mod.o : ../src/curvilinear_nodes_mod.F90 $(ODIR)globals.o $(ODIR)transformation.o $(ODIR)shape_functions_mod.o $(ODIR)basis.o $(ODIR)lapack_interfaces.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv curvilinear_nodes_mod.mod $(ODIR)
	
$(ODIR)bathymetry_interp_mod.o : ../src/bathymetry_interp_mod.F90 $(ODIR)globals.o $(ODIR)basis.o $(ODIR)transformation.o $(ODIR)shape_functions_mod.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv bathymetry_interp_mod.mod $(ODIR)
	
$(ODIR)shape_functions_mod.o : ../src/shape_functions_mod.F90     $(ODIR)globals.o $(ODIR)basis.o $(ODIR)lapack_interfaces.o $(ODIR)vandermonde.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv shape_functions_mod.mod $(ODIR)	
	
$(ODIR)transformation.o : ../src/transformation.F90               $(ODIR)basis.o $(ODIR)vandermonde.o $(ODIR)lapack_interfaces.o $(ODIR)shape_functions_mod.o
	$(FC) $(FFLAGS) -c $< -o $@
	mv transformation.mod $(ODIR)		
	
$(ODIR)vandermonde.o : ../src/vandermonde.F90                     $(ODIR)globals.o $(ODIR)basis.o 
	$(FC) $(FFLAGS) -c $< -o $@
	mv vandermonde.mod $(ODIR)		
	
$(ODIR)system.o : ../src/system.F90                                
	$(FC) $(FFLAGS) -c $< -o $@	
	
$(ODIR)read_plot_input.o : read_plot_input.F90               $(ODIR)plot_globals.o $(ODIR)plot_mod.o              
	$(FC) $(FFLAGS) -c $< -o $@		

$(ODIR)convert.o: convert2cartesian.F90                      $(ODIR)globals.o $(ODIR)grid_file_mod.o $(ODIR)spherical_mod.o
	$(FC) $(FFLAGS) -c $< -o $@
        
$(ODIR)adjust.o: adjust_channel_boundaries.f90               $(ODIR)globals.o $(ODIR)grid_file_mod.o 
	$(FC) $(FFLAGS) -c $< -o $@        

$(ODIR)plot_mod.o: plot_mod.F90                              $(ODIR)globals.o $(ODIR)plot_globals.o $(ODIR)evaluate_mod.o $(ODIR)labels_mod.o $(ODIR)axes_mod.o  $(ODIR)read_input.o $(ODIR)basis.o  $(ODIR)transformation.o
	$(FC) $(FFLAGS) -c $< -o $@ 	
	mv plot_mod.mod $(ODIR)
	
$(ODIR)evaluate_mod.o: evaluate_mod.F90                      $(ODIR)globals.o $(ODIR)triangulation.o  $(ODIR)basis.o 
	$(FC) $(FFLAGS) -c $< -o $@ 	
	mv evaluate_mod.mod $(ODIR)	
	
$(ODIR)labels_mod.o: labels_mod.F90                          $(ODIR)globals.o $(ODIR)plot_globals.o
	$(FC) $(FFLAGS) -c $< -o $@ 	
	mv labels_mod.mod $(ODIR)
	
$(ODIR)axes_mod.o: axes_mod.F90                              $(ODIR)globals.o $(ODIR)plot_globals.o
	$(FC) $(FFLAGS) -c $< -o $@ 	
	mv axes_mod.mod $(ODIR)	
	
$(ODIR)plot.o: plot_dg.F90                                   $(ODIR)plot_globals.o $(ODIR)globals.o $(ODIR)grid_file_mod.o $(ODIR)plot_mod.o $(ODIR)evaluate_mod.o $(ODIR)labels_mod.o $(ODIR)basis.o $(ODIR)read_write_output.o $(ODIR)read_input.o
	$(FC) $(FFLAGS) -c $< -o $@       
