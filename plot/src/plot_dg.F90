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
      USE evaluate_mod, ONLY: evaluate_basis
      USE labels_mod, ONLY: latex_axes_labels,run_latex,close_tex, &
                            latex_element_labels,latex_node_labels  
      USE axes_mod, ONLY: write_all_axes                            
                          
      USE triangulation, ONLY: reference_element_delaunay                            
      USE edge_connectivity_mod
      USE curvilinear_nodes_mod
      USE transformation
      USE shape_functions_mod
      USE version
      USE initialize
      
      IMPLICIT NONE
      
      INTEGER :: ord
      INTEGER :: start_snap,end_snap
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
      
      IF (sol_diff_option == 1) THEN
        CALL read_dgsweinp(sol2,input_path2,substitute_path,replace_path,sub_path)      
        CALL sizes(sol2)           
        CALL read_grid(sol2)
        CALL connectivity(sol2)
        CALL curvilinear(sol2)    
      ENDIF
      

      
      PRINT("(A)"), "Calculating additional ploting nodes..."
      nord = (p_high-p_low+1)/p_skip
      
      mnpp = (p_high+1)**2      
      ALLOCATE(r(mnpp,sol1%nel_type*nord),s(mnpp,sol1%nel_type*nord))      
      ALLOCATE(psic(sol1%mnnds,mnpp,sol1%nel_type*nord))
      ALLOCATE(rect(3,3*mnpp,sol1%nel_type*nord))      
      ALLOCATE(nptri(sol1%nel_type*nord),npplt(sol1%nel_type*nord),pplt(sol1%nel_type*nord))
      ncall = 0 
      DO et = 1,sol1%nel_type 
              
        DO ord = 1,nord
          i = (et-1)*nord+ord
          pplt(i) = (ord-1)*p_skip+p_low
          CALL element_nodes(et,space,pplt(i),npplt(i),r(:,i),s(:,i))                  
          CALL shape_functions_area_eval(et,sol1%np(et),nnd,npplt(i),r(:,i),s(:,i),psic(:,:,i))  
          CALL reference_element_delaunay(et,pplt(i),npplt(i),r(:,i),s(:,i),nptri(i),rect(:,:,i))        
          
          PRINT("(4(A,I4))"), "  number of additional nodes/sub-triangles: ", npplt(i),"/",nptri(i) 
          
!           ncall = ncall + 1
!           WRITE(nout,"(I3.3)") ncall 
!           fname = "ref_el_"//nout//".ps"
!           CALL plot_ref_el(fname,figure_width,et,np(et),nptri(i),rect(:,:,i),r(:,i),s(:,i))
          
!           DO el = 1,nptri(i)
!             PRINT "(4(I5))", el,(rect(nd,el,i), nd=1,3)
!           ENDDO
!           PRINT*, "" 
        ENDDO         
        
      ENDDO                                    
      
      PRINT("(A)"), "Calculating additional ploting point coordinates..."      
      ALLOCATE(xyplt(mnpp,sol1%ne,2))
      DO el = 1,sol1%ne      
        et = sol1%el_type(el)                          
        nnd = sol1%nnds(et)
        i = (et-1)*nord+nord
        npts = npplt(i)
        DO pt = 1,npts              
          CALL element_transformation(nnd,sol1%elxy(:,el,1),sol1%elxy(:,el,2),psic(:,pt,i),xpt,ypt)           
          xyplt(pt,el,1) = xpt
          xyplt(pt,el,2) = ypt
        ENDDO
      ENDDO       
             
                 
             
      
      PRINT("(A)"), "Evaluating reference element coordinate information..."
      CALL evaluate_basis(sol1%p,nord,mnpp,sol1%mndof,sol1%nel_type,npplt,r,s,ndof_sol,phi_sol)   
      
      
      CALL read_colormap(cmap_file)
      
      
      
      PRINT("(A)"), "Finding zoom box..."
      CALL zoom_box(sol1%ne,nord,npplt,sol1%el_type,xyplt,xbox_min,xbox_max,ybox_min,ybox_max, &
                                                     xmin,xmax,ymin,ymax,sol1%el_in)

                                                     
      PRINT("(A)"), "Scaling coordinates..."
      CALL scale_factors(figure_width,figure_height,xmin,xmax,ymin,ymax,ax,bx,ay,by)
          
      
#ifdef adcirc 

#elif dgswem

#else 
      snap_start = snap_start + 1
      snap_end = snap_end + 1      
#endif      
      
      nsnap_Z = snap_end
      nsnap_Qx = snap_end
      nsnap_Qy = snap_end
      nsnap_hb = 1
      
      t_start = 0d0
      t_end = tf*86400d0 
      t_snap = 0d0         



      CALL read_solutions(zeta,vel,bathy,sol1)   
      IF (sol_diff_option == 1) THEN
        CALL read_solutions(zeta,vel,bathy,sol2)         
      ENDIF
      
      CALL setup_cbounds(npplt,mesh,sol1,1,1)     
      CALL make_plot(1,0d0,mesh,sol1)  
      

      IF (zeta%plot_sol_option == 0 .and. vel%plot_sol_option == 0 .and. bathy%plot_sol_option == 0) THEN
        STOP
      ENDIF     
      
      IF (adapt_option == 1) THEN
        OPEN(UNIT = 998, FILE="error.out", STATUS="REPLACE")
        WRITE(998,"(A)") "name     snap     error_total     nptri_total     pplt_max     ne_total"
      ENDIF      
      
      CALL setup_cbounds(npplt,bathy,sol1,1,1)      
      CALL make_plot(1,t_snap,bathy,sol1,sol2)    
      
      CALL setup_cbounds(npplt,cfl,sol1,1,1)      
      CALL make_plot(1,t_snap,cfl,sol1)          
      
      
      IF (zeta%plot_sol_option == 0 .and. vel%plot_sol_option == 0) THEN
        STOP
      ENDIF
      
      
             
         
      IF ((snap_start > nsnap_Z) .or. (snap_start > nsnap_Qx) .or. (snap_start > nsnap_Qy)) THEN
        snap_start = MIN(nsnap_Z,nsnap_Qx,nsnap_Qy)
      ENDIF                  
      
      IF ((snap_end > nsnap_Z) .or. (snap_end > nsnap_Qx) .or. (snap_end > nsnap_Qy)) THEN
        snap_end = MIN(nsnap_Z,nsnap_Qx,nsnap_Qy)
      ENDIF
      
      
      


      CALL setup_cbounds(npplt,zeta,sol1,snap_start,snap_end)
      CALL setup_cbounds(npplt,vel,sol1,snap_start,snap_end)

      
      PRINT("(A)"), " "
      
      
      DO snap = snap_start,snap_end
      
        
        t_snap = sol1%t(snap)      
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
      
