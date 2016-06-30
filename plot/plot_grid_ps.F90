      PROGRAM plot_grid_ps
      
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
      USE read_dginp, ONLY: read_input,out_direc,p,ctp, &
                            grid_file,curve_file,cb_file_exists
      USE plot_mod, ONLY: read_colormap,write_psheader,close_ps, &
                          scale_coordinates,zoom_box, &
                          evaluate_depth_solution,evaluate_velocity_solution, &
                          plot_contours,plot_mesh
      USE triangulation, ONLY: reference_element_delaunay
      USE edge_connectivity_mod
      USE curvilinear_nodes_mod
      USE transformation
      USE shape_functions_mod
      
      IMPLICIT NONE
      

      
      CALL read_plot_input()
      
      IF (plot_mesh_option == 0 .and. plot_zeta_option == 0 .and. &
          plot_bathy_option == 0 .and. plot_vel_option == 0) THEN
        
        PRINT("(A)"), "No plot options have been specified"
        STOP
      ENDIF    
      
      CALL read_input(0,input_path)
      
      CALL sizes()

      PRINT*, grid_file
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

      PRINT*, "Calculating curved boundary information..."
      CALL shape_functions_linear_at_ctp(nel_type,np,psiv)                   
      CALL eval_coordinates_curved(ctp,nnds,nverts,el_type,xy,ect,fbseg,fbnds, &
                                   nnfbed,nfbedn,nfbednn,ged2el,ged2led, &
                                   psiv,bndxy,elxy)     
      
      PRINT*, "Calculating additional ploting point coordinates..."
      space = 1  
      ALLOCATE(r(mnpp),s(mnpp))      
      ALLOCATE(psic(mnnds,mnpp,nel_type))
      DO et = 1,nel_type     
        CALL element_nodes(et,space,pplt(et),npts,r,s)                  
        CALL shape_functions_area_eval(et,np(et),nnd,npts,r,s,psic(:,:,et))     
      ENDDO                                    
           
      ALLOCATE(xyplt(mnpp,ne,2))
      DO el = 1,ne      
        et = el_type(el)                          
        nnd = nnds(et)
        npts = nplt(et)
        DO pt = 1,npts              
          CALL element_transformation(nnd,elxy(:,el,1),elxy(:,el,2),psic(:,pt,et),xpt,ypt)           
          xyplt(pt,el,1) = xpt
          xyplt(pt,el,2) = ypt
        ENDDO
      ENDDO       
             
      
      PRINT*, "Evaluating reference element coordinate information..."
      ALLOCATE(phi(mndof,mnpp,nel_type))
      ALLOCATE(rect(3,3*mnpp,nel_type))
      DO et = 1,nel_type
        CALL element_nodes(et,1,pplt(et),n,r,s)
        CALL element_basis(et,p,ndf,n,r,s,phi(:,:,et))
        CALL reference_element_delaunay(n,r,s,ntri(et),rect(:,:,et))
        
!         DO i = 1,3
!           PRINT "(*(I5))", (rect(i,j,et), j = 1,ntri(et))
!         ENDDO    
!         PRINT*, ""
      ENDDO      


      
      PRINT*, "Finding zoom box..."
      CALL zoom_box(ne,el_type,nplt,xyplt,xbox_min,xbox_max,ybox_min,ybox_max, &
                                                     xmin,xmax,ymin,ymax,el_in)

      PRINT*, "Scaling coordinates..."
      CALL scale_coordinates(ne,nn,el_type,nverts,nplt,figure_width,xmin,xmax,ymin,ymax,xyplt,xy)

      
      
      CALL read_colormap(cmap_file)

      
      

      snap_start = snap_start + 1
      snap_end = snap_end + 1      
      
      nsnap_Z = snap_end
      nsnap_Qx = snap_end
      nsnap_Qy = snap_end
        
      IF (plot_zeta_option == 1 .or. plot_vel_option == 1) THEN
        PRINT*, "Reading zeta solution..."          
        CALL read_solution_full(out_direc,"Z.sol","N",t,Z,nsnap_Z) 
      ENDIF
      IF (plot_vel_option == 1) THEN
        PRINT*, "Reading Qx and Qy solutions..."       
        CALL read_solution_full(out_direc,"Qx.sol","N",t,Qx,nsnap_Qx)        
        CALL read_solution_full(out_direc,"Qy.sol","N",t,Qy,nsnap_Qy)  
      ENDIF
      IF (plot_bathy_option == 1 .or. plot_vel_option == 1) THEN
        PRINT*, "Reading bathymetry solution..."          
        CALL read_solution_full(out_direc,"hb.sol","N",t,hb,nsnap_hb)      
      ENDIF   
      
      
      
      IF (plot_bathy_option == 1 .or. plot_vel_option == 1) THEN
        PRINT*, "Evaluating bathymetry at additional plotting points..."
        CALL evaluate_depth_solution(ne,el_type,el_in,nplt,ndof,phi,1,hb,hb_val)
      ENDIF      
      
      
      
      IF (plot_bathy_option == 1) THEN
        PRINT*, "Writing bathymetry PostScript file..."
        CALL write_psheader("bathy.ps",hb_unit)
        CALL plot_contours(hb_unit,nplt,ntri,rect,ne,el_type,el_in,xyplt,hb_val)     
        IF (plot_mesh_option == 1) THEN
          CALL plot_mesh(hb_unit,ne,nverts,el_type,el_in,xy,ect)   
        ENDIF
        CALL close_ps(hb_unit)
      ENDIF
      
      IF (plot_mesh_option == 1 .and. plot_zeta_option == 0 .and. &
          plot_bathy_option == 0 .and. plot_vel_option == 0 ) THEN          
        PRINT*, "Writing mesh PostScript file..."   
        CALL write_psheader("mesh.ps",mesh_unit)            
        CALL plot_mesh(mesh_unit,ne,nverts,el_type,el_in,xy,ect)      
        CALL close_ps(mesh_unit)                      
      ENDIF
      
      IF (plot_zeta_option == 0 .and. plot_vel_option == 0) THEN
        STOP
      ENDIF
      
      
             
         
      IF ((snap_start > nsnap_Z) .or. (snap_start > nsnap_Qx) .or. (snap_start > nsnap_Qy)) THEN
        snap_start = MIN(nsnap_Z,nsnap_Qx,nsnap_Qy)
      ENDIF                  
      
      IF ((snap_end > nsnap_Z) .or. (snap_end > nsnap_Qx) .or. (snap_end > nsnap_Qy)) THEN
        snap_end = MIN(nsnap_Z,nsnap_Qx,nsnap_Qy)
      ENDIF
      
      
      
      snap_char = "0000"
      
      DO snap = snap_start,snap_end
      
        PRINT*, "Time snap: ", snap-1        
        WRITE(snap_char,"(I4.4)") snap-1

        IF (plot_zeta_option == 1 .or. plot_vel_option == 1) THEN
          PRINT*, "  Evaluating zeta at additional plotting points..."
          CALL evaluate_depth_solution(ne,el_type,el_in,nplt,ndof,phi,snap,Z,Z_val)
        ENDIF
        
        IF (plot_zeta_option == 1) THEN        
          PRINT*, "  Writing zeta PostScript file..."        
          CALL write_psheader("zeta_"//snap_char//".ps",Z_unit)            
          CALL plot_contours(Z_unit,nplt,ntri,rect,ne,el_type,el_in,xyplt,Z_val)      
          IF (plot_mesh_option == 1) THEN
            CALL plot_mesh(Z_unit,ne,nverts,el_type,el_in,xy,ect)
          ENDIF      
          CALL close_ps(Z_unit)        
        ENDIF        
      
        IF (plot_vel_option == 1) THEN
          PRINT*, "  Evaluating velocity at additional plotting points..."
          CALL evaluate_velocity_solution(ne,el_type,el_in,nplt,ndof,phi,snap,Qx,Qy,Z_val,hb_val,vel_val)
        ENDIF       
        
        IF (plot_vel_option == 1) THEN       
          PRINT*, "  Writing velocity PostScript file..."        
          CALL write_psheader("vel_"//snap_char//".ps",vel_unit)            
          CALL plot_contours(vel_unit,nplt,ntri,rect,ne,el_type,el_in,xyplt,vel_val)      
          IF (plot_mesh_option == 1) THEN
            CALL plot_mesh(vel_unit,ne,nverts,el_type,el_in,xy,ect)
          ENDIF      
          CALL close_ps(vel_unit)        
        ENDIF
            
      ENDDO

      END PROGRAM plot_grid_ps
      
