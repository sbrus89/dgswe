      PROGRAM plot_dg
      
      USE plot_globals
      USE globals, ONLY: ne,nn,nverts,ndof,mndof,np,mnp,nnds,mnnds,nel_type,ect,xy,depth, &
                         nope,neta,obseg,obnds,nvel,nbou,fbseg,fbnds,bndxy,grid_name, &
                         el_type,elxy,elhb, &
                         nepn,epn,mnepn,ned,ged2el,ged2el,ged2led,ged2nn,ed_type,recv_edge, &
                         nied,iedn,nobed,obedn,nfbed,fbedn,nnfbed,nfbedn,nfbednn, &
                         psic,psiv, &
                         deg2rad, &
                         nsta,xysta
      USE basis, ONLY: element_nodes,element_basis
      USE read_dginp, ONLY: read_input,out_direc,p,ctp,hbp,tf, &
                            grid_file,curve_file,cb_file_exists,bathy_file,hb_file_exists, &
                            sphi0,slam0, &
                            sta_opt,stations_file
      USE plot_mod, ONLY: read_colormap,setup_cbounds,plot_ref_el, &
                          scale_factors,zoom_box,make_plot,make_movie, &
                          plot_filled_contours
      USE evaluate_mod, ONLY: evaluate_basis,evaluate_plotting_nodes
      USE labels_mod, ONLY: latex_axes_labels,run_latex,close_tex, &
                            latex_element_labels,latex_node_labels  
      USE axes_mod, ONLY: write_all_axes   
      USE edge_connectivity_mod
      USE curvilinear_nodes_mod
      USE transformation
      USE shape_functions_mod
      USE version
      USE initialize
      USE google_map
      
      IMPLICIT NONE
      
      INTEGER :: ord
      INTEGER :: start_snap,end_snap
      INTEGER :: sol_snap
      INTEGER :: ncall
      CHARACTER(3) :: nout
      CHARACTER(25) :: fname
      REAL(rp) :: H
      LOGICAL :: file_exists

      
      space = 0            

      
      CALL version_information(6)
      
      CALL read_plot_input()
      CALL setup_plot_types()      
  
      
      IF (mesh%plot_sol_option == 0 .and. zeta%plot_sol_option == 0 .and. &
          bathy%plot_sol_option == 0 .and. vel%plot_sol_option == 0) THEN
        
        PRINT("(A)"), "No plot options have been specified"
        STOP
      ENDIF    
      
 

      CALL read_dgsweinp(sol1,input_path,substitute_path,replace_path,sub_path)      
      CALL sizes(sol1)           
      CALL read_grid(sol1)
      CALL connectivity(sol1)
      CALL curvilinear(sol1)
      CALL find_output_type(sol1)
      
      IF (sol_diff_option == 1) THEN
        CALL read_dgsweinp(sol2,input_path2,substitute_path,replace_path,sub_path)      
        CALL sizes(sol2)           
        CALL read_grid(sol2)
        CALL connectivity(sol2)
        CALL curvilinear(sol2)    
        CALL find_output_type(sol2)        
      ENDIF
      
      spherical_flag = 0
      IF (abs(sol1%slam0) > 0d0 .and. abs(sol1%sphi0) > 0d0) THEN
        spherical_flag = 1
      ENDIF
      
      PRINT("(A)"), "Calculating additional ploting nodes..."
      nord = (p_high-p_low+1)/p_skip      
      mnpp = (p_high+1)**2      
      sol1%nord = 1
      sol1%mnpp = mnpp
      CALL evaluate_plotting_nodes(sol1%nel_type,p_high,p_skip,sol1%nord,sol1%mnpp,sol1%np,sol1%mnnds, &
                                   sol1%pplt,sol1%npplt,sol1%r,sol1%s,sol1%psic,sol1%nptri,sol1%rect)                                
      IF (sol_diff_option == 1) THEN
        PRINT("(A)"), "Calculating additional ploting nodes..."
        sol2%nord = 1
        sol2%mnpp = mnpp
        CALL evaluate_plotting_nodes(sol2%nel_type,p_high,p_skip,sol2%nord,sol2%mnpp,sol2%np,sol2%mnnds, &
                                     sol2%pplt,sol2%npplt,sol2%r,sol2%s,sol2%psic,sol2%nptri,sol2%rect)    
      ENDIF


      
      PRINT("(A)"), "Calculating additional ploting point coordinates..."      
      ALLOCATE(xyplt(mnpp,sol1%ne,2))
      DO el = 1,sol1%ne      
        et = sol1%el_type(el)                          
        nnd = sol1%nnds(et)
        npts = sol1%npplt(et)
        DO pt = 1,npts              
          CALL element_transformation(nnd,sol1%elxy(:,el,1),sol1%elxy(:,el,2),sol1%psic(:,pt,et),xpt,ypt)           
          xyplt(pt,el,1) = xpt
          xyplt(pt,el,2) = ypt
        ENDDO
      ENDDO       
             
                 
             
      
      PRINT("(A)"), "Evaluating reference element coordinate information..."
      CALL evaluate_basis(sol1%output_type,sol1%p,sol1%nord,mnpp,sol1%mndof,sol1%nel_type,sol1%npplt,sol1%r,sol1%s,ndof_sol,phi_sol)   
      
      
      CALL read_colormap(cmap_file)
      
      
      
      PRINT("(A)"), "Finding zoom box..."
      CALL zoom_box(sol1%ne,sol1%npplt,sol1%nnds,sol1%el_type,sol1%psic,sol1%elxy,xbox_min,xbox_max,ybox_min,ybox_max, &
                                                                        xmin,xmax,ymin,ymax,sol1%el_in)
      IF (sol_diff_option == 1) THEN
        PRINT("(A)"), "Finding zoom box..."      
        CALL zoom_box(sol2%ne,sol2%npplt,sol2%nnds,sol2%el_type,sol2%psic,sol2%elxy,xbox_min,xbox_max,ybox_min,ybox_max, &      
                                                                          xmin,xmax,ymin,ymax,sol2%el_in)      
      ENDIF

                                                     
      PRINT("(A)"), "Scaling coordinates..."
      CALL scale_factors(figure_width,figure_height,xmin,xmax,ymin,ymax,ax,bx,ay,by)
                   
      IF (plot_google_map == 1 .and. spherical_flag == 1) THEN
        PRINT("(A)"), "Downloading and satellite image from Google Maps..."
        CALL get_map(xmin,xmax,ymin,ymax,sol1%slam0,sol1%sphi0,lamc,phic,map,map_height,map_width,map_res)    
      ENDIF
      
      
      
      t_start = 0d0
      t_end = tf*86400d0 
      t_snap = 0d0         



      CALL read_solutions(zeta,vel,bathy,sol1)   
      IF (sol_diff_option == 1) THEN
        CALL read_solutions(zeta,vel,bathy,sol2)         
      ENDIF
      
      CALL setup_cbounds(mesh,sol1,1,1)     
      CALL make_plot(1,0d0,mesh,sol1,sol2)  
      

      IF (zeta%plot_sol_option == 0 .and. vel%plot_sol_option == 0 .and. bathy%plot_sol_option == 0) THEN
        STOP
      ENDIF     
      
      IF (adapt_option == 1) THEN
        OPEN(UNIT = 998, FILE="error.out", STATUS="REPLACE")
        WRITE(998,"(A)") "name     snap     error_total     nptri_total     pplt_max     ne_total"
      ENDIF      
      
      CALL setup_cbounds(bathy,sol1,1,1)      
      CALL make_plot(1,t_snap,bathy,sol1,sol2)    
      
      CALL setup_cbounds(cfl,sol1,1,1)      
      CALL make_plot(1,t_snap,cfl,sol1,sol2)          
      
      
      IF (zeta%plot_sol_option == 0 .and. vel%plot_sol_option == 0) THEN
        STOP
      ENDIF
      
      
             
         
!       IF ((snap_start > sol1%nsnap_Z) .or. (snap_start > sol1%nsnap_Qx) .or. (snap_start > sol1%nsnap_Qy)) THEN
!         snap_start = MIN(sol1%nsnap_Z,sol1%nsnap_Qx,sol1%nsnap_Qy)
!       ENDIF                  
!       
!       IF ((snap_end > sol1%nsnap_Z) .or. (snap_end > sol1%nsnap_Qx) .or. (snap_end > sol1%nsnap_Qy)) THEN
!         snap_end = MIN(sol1%nsnap_Z,sol1%nsnap_Qx,sol1%nsnap_Qy)
!       ENDIF
      
      
      

      CALL setup_cbounds(zeta,sol1,snap_start,snap_end)
      CALL setup_cbounds(vel,sol1,snap_start,snap_end)

      
      PRINT("(A)"), " "
      
      
      DO snap = snap_start,snap_end
      
        ! Account for initial condition in dgswe output              
        IF (sol1%output_type == "dgswe") THEN        
          sol_snap = snap + 1        
        ELSE         
          sol_snap = snap        
        ENDIF        
        
        t_snap = sol1%t(sol_snap)      
        PRINT("(A)"), "---------------------------------------------"
        PRINT("(A,I5,A,I5,A,F20.5)"), "Time snap: ", snap,"/",snap_end," t = ", t_snap        
     
        IF (zeta%plot_max_option == 2) THEN
          CALL plot_filled_contours(-1,snap,zeta,sol1,sol2)
        ELSE
          CALL make_plot(snap,t_snap,zeta,sol1,sol2)
        ENDIF
              
        PRINT("(A)"), " "        

        IF (vel%plot_max_option == 2) THEN
          CALL plot_filled_contours(-1,snap,vel,sol1,sol2)
        ELSE
          CALL make_plot(snap,t_snap,vel,sol1,sol2)
        ENDIF
        
        
      ENDDO
      
      
      
      ! Plot maximums
      PRINT("(A)"), "---------------------------------------------"   
      snap = snap_end+1
      IF (zeta%plot_max_option > 0) THEN
        CALL make_plot(snap,t_snap,zeta,sol1,sol2)
      ENDIF
      PRINT("(A)"), " "        
      IF (vel%plot_max_option > 0) THEN
        CALL make_plot(snap,t_snap,vel,sol1,sol2)
      ENDIF 
      
      
      
      
      IF (adapt_option == 1) THEN
        CLOSE(998)
      ENDIF
      
      CALL make_movie(zeta,frmt)
      CALL make_movie(vel,frmt)
  
      PRINT("(A)"), " "      
      END PROGRAM plot_dg
      
