      PROGRAM plot_dg
      
      USE plot_globals
      USE globals, ONLY: ne,nn,nverts,ndof,mndof,np,mnp,nnds,mnnds,nel_type,ect,xy,depth, &
                         nope,neta,obseg,obnds,nvel,nbou,fbseg,fbnds,bndxy,grid_name, &
                         el_type,elxy,elhb, &
                         nepn,epn,mnepn,ned,ged2el,ged2el,ged2led,ged2nn,ed_type,recv_edge, &
                         nied,iedn,nobed,obedn,nfbed,fbedn,nnfbed,nfbedn,nfbednn, &
                         psic,psiv
      USE grid_file_mod
      USE basis, ONLY: element_nodes,element_basis
      USE read_write_output, ONLY: read_solution_full
      USE read_dginp, ONLY: read_input,out_direc,p,ctp,hbp,tf, &
                            grid_file,curve_file,cb_file_exists
      USE plot_mod, ONLY: read_colormap,setup_cbounds, &
                          scale_coordinates,zoom_box,make_plot,make_movie                                                                           
      USE evaluate_mod, ONLY: evaluate_depth_solution,evaluate_velocity_solution, &
                              evaluate_basis,find_solution_minmax
      USE labels_mod, ONLY: latex_axes_labels,run_latex,close_tex, &
                            latex_element_labels,latex_node_labels  
      USE axes_mod, ONLY: write_all_axes                            
                          
      USE triangulation, ONLY: reference_element_delaunay                            
      USE edge_connectivity_mod
      USE curvilinear_nodes_mod
      USE transformation
      USE shape_functions_mod
      USE version
      
      IMPLICIT NONE
      
      INTEGER :: start_snap,end_snap
      
      space = 1  
      
      CALL version_information(6)
      
      CALL read_plot_input()
      
      IF (mesh%plot_sol_option == 0 .and. zeta%plot_sol_option == 0 .and. &
          bathy%plot_sol_option == 0 .and. vel%plot_sol_option == 0) THEN
        
        PRINT("(A)"), "No plot options have been specified"
        STOP
      ENDIF    
      
      CALL read_input(0,input_path)
      IF (substitute_path == 1) THEN
        CALL substitute_partial_path(grid_file,replace_path,sub_path)
        CALL substitute_partial_path(curve_file,replace_path,sub_path)    
      ENDIF
      
      CALL sizes()
      CALL setup_plot_types()

      PRINT("(A)"), grid_file
      CALL read_header(0,grid_file,grid_name,ne,nn)        
      CALL read_coords(nn,xy,depth)
      CALL read_connectivity(ne,ect,el_type) 
      CALL init_element_coordinates(ne,ctp,el_type,nverts,xy,ect,elxy)                  
      CALL read_open_boundaries(nope,neta,obseg,obnds)            
      CALL read_flow_boundaries(nbou,nvel,fbseg,fbnds)                       
      CALL read_curve_file(0,curve_file,ctp,nbou,xy,bndxy,cb_file_exists)      
      CALL print_grid_info(grid_file,grid_name,ne,nn)    
      
      CALL elements_per_node(ne,nn,nverts,el_type,ect,nepn,mnepn,epn)       
      CALL find_edge_pairs(ne,nverts,el_type,ect,nepn,epn,ned,ged2el,ged2nn,ged2led)      
      CALL find_interior_edges(ned,ged2el,nied,iedn,ed_type,recv_edge)      
      CALL find_open_edges(nope,obseg,obnds,ged2nn,nobed,obedn,ed_type,recv_edge)            
      CALL find_flow_edges(nbou,fbseg,fbnds,ged2nn,nnfbed,nfbedn,nfbednn,nfbed,fbedn,recv_edge,ed_type)     
      nred = 0
      CALL print_connect_info(mnepn,ned,nied,nobed,nfbed,nnfbed,nred)

      PRINT("(A)"), "Calculating curved boundary information..."
      CALL shape_functions_linear_at_ctp(nel_type,np,psiv)                   
      CALL eval_coordinates_curved(ctp,nnds,nverts,el_type,xy,ect,fbseg,fbnds, &
                                   nnfbed,nfbedn,nfbednn,ged2el,ged2led, &
                                   psiv,bndxy,elxy)     
      
!       CALL find_element_init(nel_type,nverts,np,nnds,nn,xy,nepn,epn)      
      
      PRINT("(A)"), "Calculating additional ploting point coordinates..."
      ALLOCATE(r(mnpp,nel_type),s(mnpp,nel_type))      
      ALLOCATE(psic(mnnds,mnpp,nel_type))
      ALLOCATE(rect(3,3*mnpp,nel_type))      
      DO et = 1,nel_type     
        CALL element_nodes(et,space,pplt(et),npts,r(:,et),s(:,et))                  
        CALL shape_functions_area_eval(et,np(et),nnd,npts,r(:,et),s(:,et),psic(:,:,et))  
        CALL reference_element_delaunay(npts,r(:,et),s(:,et),nptri(et),rect(:,:,et))        
        
!         DO i = 1,3
!           PRINT "(*(I5))", (rect(i,j,et), j = 1,nptri(et))
!         ENDDO    
!         PRINT*, ""

        PRINT("(4(A,I4))"), "  number of additional nodes/sub-triangles: ", npts,"/",nptri(et)
        
      ENDDO                                    
           
      ALLOCATE(xyplt(mnpp,ne,2))
      DO el = 1,ne      
        et = el_type(el)                          
        nnd = nnds(et)
        npts = npplt(et)
        DO pt = 1,npts              
          CALL element_transformation(nnd,elxy(:,el,1),elxy(:,el,2),psic(:,pt,et),xpt,ypt)           
          xyplt(pt,el,1) = xpt
          xyplt(pt,el,2) = ypt
        ENDDO
      ENDDO       
             
      
      PRINT("(A)"), "Evaluating reference element coordinate information..."
      CALL evaluate_basis(p,mnpp,mndof,nel_type,npplt,r,s,zeta)
      CALL evaluate_basis(hbp,mnpp,mndof,nel_type,npplt,r,s,bathy)
      CALL evaluate_basis(p,mnpp,mndof,nel_type,npplt,r,s,vel)
      
      
      
      CALL read_colormap(cmap_file)
      
      
      
      PRINT("(A)"), "Finding zoom box..."
      CALL zoom_box(ne,el_type,npplt,xyplt,xbox_min,xbox_max,ybox_min,ybox_max, &
                                                     xmin,xmax,ymin,ymax,el_in)
      PRINT("(A)"), "Scaling coordinates..."
      CALL scale_coordinates(ne,nn,el_type,nverts,nnds,npplt,figure_width,xmin,xmax,ymin,ymax,xyplt,xy,elxy)

      
    


      

      snap_start = snap_start + 1
      snap_end = snap_end + 1      
      
      nsnap_Z = snap_end
      nsnap_Qx = snap_end
      nsnap_Qy = snap_end
      nsnap_hb = 1
      
      t_start = 0d0
      t_end = tf*86400d0 
      t_snap = 0d0
      
      snap_char = "0000"      
      
          
        
      IF (zeta%plot_sol_option == 1 .or. vel%plot_sol_option == 1) THEN
        PRINT("(A)"), "Reading zeta solution..."          
        CALL read_solution_full(out_direc,"Z.sol","N",t,Z,nsnap_Z) 
      ENDIF
      IF (vel%plot_sol_option == 1) THEN
        PRINT("(A)"), "Reading Qx and Qy solutions..."       
        CALL read_solution_full(out_direc,"Qx.sol","N",t,Qx,nsnap_Qx)        
        CALL read_solution_full(out_direc,"Qy.sol","N",t,Qy,nsnap_Qy) 
      ELSE 
        ALLOCATE(Qx(1,1,1),Qy(1,1,1))
      ENDIF  
      IF (bathy%plot_sol_option == 1 .or. vel%plot_sol_option == 1) THEN
        PRINT("(A)"), "Reading bathymetry solution..."          
        CALL read_solution_full(out_direc,"hb.sol","N",t,hb,nsnap_hb)  
      ENDIF  
      

      
      

      
      CALL make_plot(-1,t_snap,bathy,hb(:,:,1))                    
      
      CALL make_plot(-1,0d0,mesh)
      
      IF (zeta%plot_sol_option == 0 .and. vel%plot_sol_option == 0) THEN
        STOP
      ENDIF
      
      
             
         
      IF ((snap_start > nsnap_Z) .or. (snap_start > nsnap_Qx) .or. (snap_start > nsnap_Qy)) THEN
        snap_start = MIN(nsnap_Z,nsnap_Qx,nsnap_Qy)
      ENDIF                  
      
      IF ((snap_end > nsnap_Z) .or. (snap_end > nsnap_Qx) .or. (snap_end > nsnap_Qy)) THEN
        snap_end = MIN(nsnap_Z,nsnap_Qx,nsnap_Qy)
      ENDIF
      
      
      

      zeta%sol_max = -1d10
      zeta%sol_min = 1d10
      vel%sol_max = -1d10
      vel%sol_min = 1d10        
      
      IF (zeta%cscale_option == "auto-all" .or. vel%cscale_option == "auto-all") THEN
        DO snap = snap_start,snap_end
        
          CALL find_solution_minmax(ne,el_type,el_in,npplt,zeta,Z(:,:,snap))
        
          CALL find_solution_minmax(ne,el_type,el_in,npplt,vel,zeta%sol_val,bathy%sol_val,Qx(:,:,snap),Qy(:,:,snap))
        
        ENDDO
      ENDIF

      CALL setup_cbounds(zeta,snap_start,snap_end)
      CALL setup_cbounds(vel,snap_start,snap_end)

      
      PRINT("(A)"), " "
      
      

      
      DO snap = snap_start,snap_end
      
        
        t_snap = t(snap)      
        PRINT("(A)"), "---------------------------------------------"
        PRINT("(A,I5,A,I5,A,F20.5)"), "Time snap: ", snap-1,"/",snap_end," t = ", t_snap        
     
      
        CALL make_plot(snap,t_snap,zeta,Z(:,:,snap))
              
        PRINT("(A)"), " "        

        CALL make_plot(snap,t_snap,vel,zeta%sol_val,bathy%sol_val,Qx(:,:,snap),Qy(:,:,snap))
        
      ENDDO
      
      
      CALL make_movie(zeta,frmt)
      CALL make_movie(vel,frmt)
  

      END PROGRAM plot_dg
      
