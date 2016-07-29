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
      USE plot_mod, ONLY: read_colormap,write_psheader,close_ps,convert_ps, &
                          scale_coordinates,zoom_box,make_plot, &                          
                          plot_contours,plot_mesh                          
      USE evaluate_mod, ONLY: evaluate_depth_solution,evaluate_velocity_solution, &
                              evaluate_basis,find_solution_minmax
      USE labels_mod, ONLY: latex_axes_labels,run_latex,close_tex, &
                            latex_element_labels,latex_node_labels  
      USE axes_mod, ONLY: write_all_axes                            
                          
      USE edge_connectivity_mod
      USE curvilinear_nodes_mod
      USE transformation
      USE shape_functions_mod
      
      IMPLICIT NONE
      
      INTEGER :: start_snap,end_snap
      
      space = 0  
      
      CALL read_plot_input()
      
      IF (mesh%plot_sol_option == 0 .and. zeta%plot_sol_option == 0 .and. &
          bathy%plot_sol_option == 0 .and. vel%plot_sol_option == 0) THEN
        
        PRINT("(A)"), "No plot options have been specified"
        STOP
      ENDIF    
      
      CALL read_input(0,input_path)
      
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
             
      
      PRINT("(A)"), "Evaluating reference element coordinate information..."
!       ALLOCATE(phi(mndof,mnpp,nel_type))
!       ALLOCATE(phi_hb(mndof,mnpp,nel_type))      
!       ALLOCATE(rect(3,3*mnpp,nel_type))
!       DO et = 1,nel_type
!         CALL element_nodes(et,space,pplt(et),n,r,s)
!         CALL element_basis(et,p,ndf,n,r,s,phi(:,:,et))
!         CALL element_basis(et,hbp,ndf,n,r,s,phi_hb(:,:,et))        
!         CALL reference_element_delaunay(n,r,s,ntri(et),rect(:,:,et))
!         
! !         DO i = 1,3
! !           PRINT "(*(I5))", (rect(i,j,et), j = 1,ntri(et))
! !         ENDDO    
! !         PRINT*, ""
!         PRINT("(A,I4,A,I4)"), "  number of additional nodes/sub-triangles: ", n,"/",ntri(et)
!       ENDDO      


      CALL evaluate_basis(p,space,mnpp,mndof,nel_type,pplt,zeta)
      CALL evaluate_basis(hbp,space,mnpp,mndof,nel_type,pplt,bathy)
      CALL evaluate_basis(p,space,mnpp,mndof,nel_type,pplt,vel)
      
      
      
      CALL read_colormap(cmap_file)
      
      
      
      PRINT("(A)"), "Finding zoom box..."
      CALL zoom_box(ne,el_type,nplt,xyplt,xbox_min,xbox_max,ybox_min,ybox_max, &
                                                     xmin,xmax,ymin,ymax,el_in)
                                                     
!       PRINT*, xmin,xmax,ymin,ymax                                                     

      PRINT("(A)"), "Scaling coordinates..."
      CALL scale_coordinates(ne,nn,el_type,nverts,nplt,figure_width,xmin,xmax,ymin,ymax,xyplt,xy)

      
    


      

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
      ENDIF  
      IF (bathy%plot_sol_option == 1 .or. vel%plot_sol_option == 1) THEN
        PRINT("(A)"), "Reading bathymetry solution..."          
        CALL read_solution_full(out_direc,"hb.sol","N",t,hb,nsnap_hb)  
      ENDIF  
      

      
      
!       IF (bathy%plot_sol_option == 1 .or. vel%plot_sol_option == 1) THEN
!         PRINT("(A)"), "Evaluating bathymetry at additional plotting points..."
!         CALL evaluate_depth_solution(ne,el_type,el_in,nplt,ndof_hb,phi_hb,1,hb,bathy%sol_val,bathy%sol_min,bathy%sol_max)
!       ENDIF      
      
      
      
!       IF (bathy%plot_sol_option == 1) THEN
!         PRINT("(A)"), "Writing bathymetry PostScript file..."
!         CALL latex_axes_labels(bathy%sol_min,bathy%sol_max,"bathymetry") 
!         CALL run_latex()
!         CALL write_psheader("bathy.ps",bathy%ps_unit)              
!         CALL plot_contours(bathy%ps_unit,nplt,ntri,rect,ne,el_type,el_in,xyplt,bathy%sol_val,bathy%sol_min,bathy%sol_max)     
!         IF (bathy%plot_mesh_option == 1) THEN
!           CALL plot_mesh(bathy%ps_unit,ne,nverts,el_type,el_in,xy,ect)   
!         ENDIF
!         CALL write_all_axes(bathy%ps_unit,"bathy")             
!         CALL close_ps("bathy",bathy%ps_unit)
!         CALL convert_ps("bathy",frmt,density,rm_ps)
!       ENDIF
      
      CALL make_plot(-1,t_snap,bathy,hb(:,:,1))              

      
      IF (mesh%plot_mesh_option == 1) THEN          
        PRINT("(A)"), "Writing mesh PostScript file..."   
        CALL latex_axes_labels()  
        CALL latex_element_labels(ne,mesh%el_label_option,el_type,el_in,nverts,xy,ect,nnfbed,nfbedn,ged2el)
        CALL latex_node_labels(nn,mesh%nd_label_option,xy,nbou,fbseg,fbnds)        
        CALL run_latex()        
        CALL write_psheader("mesh.ps",mesh%ps_unit)            
        CALL plot_mesh(mesh%ps_unit,ne,nverts,el_type,el_in,xy,ect)   
        CALL write_all_axes(mesh%ps_unit)          
        CALL close_ps("mesh",mesh%ps_unit)         
        CALL convert_ps("mesh",frmt,density,0)
      ENDIF
      
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
        
!       DO snap = snap_start,snap_end
!       
!         IF ((zeta%plot_sol_option == 1 .AND. zeta%cscale_option == "auto-all") .OR.  &
!             (vel%plot_sol_option == 1 .AND. vel%cscale_option == "auto-all")) THEN
!               
!           CALL evaluate_depth_solution(ne,el_type,el_in,nplt,ndof,phi,snap,Z,zeta%sol_val,Zsnap_min,Zsnap_max)
!         
!           IF (Zsnap_min < zeta%sol_min) THEN
!             zeta%sol_min = Zsnap_min
!           ENDIF
!           
!           IF (Zsnap_max > zeta%sol_max) THEN
!             zeta%sol_max = Zsnap_max
!           ENDIF
!         ENDIF
!           
! 
!         IF (vel%plot_sol_option == 1 .AND. vel%cscale_option == "auto-all") THEN
!           CALL evaluate_velocity_solution(ne,el_type,el_in,nplt,ndof,phi,snap,Qx,Qy,zeta%sol_val,bathy%sol_val,vel%sol_val,velsnap_min,velsnap_max)
!         
!           IF (velsnap_min < vel%sol_min) THEN
!             vel%sol_min = velsnap_min          
!           ENDIF
!          
!           IF (velsnap_max > vel%sol_max) THEN
!             vel%sol_max = velsnap_max
!           ENDIF
!         ENDIF
!       ENDDO
!       
      
      IF (zeta%cscale_option == "auto-all" .or. vel%cscale_option == "auto-all") THEN
        DO snap = snap_start,snap_end
        
          CALL find_solution_minmax(ne,el_type,el_in,zeta,Z(:,:,snap))
        
          CALL find_solution_minmax(ne,el_type,el_in,vel,zeta%sol_val,bathy%sol_val,Qx(:,:,snap),Qy(:,:,snap))
        
        ENDDO
      ENDIF
      
      
      
      IF (zeta%cscale_option == "file") THEN
        READ(zeta%cscale_unit,*) start_snap,end_snap
        zeta%num_cscale_vals = end_snap - start_snap + 1
        ALLOCATE(zeta%cscale_vals(zeta%num_cscale_vals,3))
        DO i = 1,zeta%num_cscale_vals
          READ(zeta%cscale_unit,*) zeta%cscale_vals(i,1),zeta%cscale_vals(i,2),zeta%cscale_vals(i,3)
        ENDDO
        CLOSE(zeta%cscale_unit)
      ENDIF
      
      IF (vel%cscale_option == "file") THEN
        READ(vel%cscale_unit,*) start_snap,end_snap
        vel%num_cscale_vals = end_snap - start_snap + 1
        ALLOCATE(vel%cscale_vals(vel%num_cscale_vals,3))
        DO i = 1,vel%num_cscale_vals
          READ(vel%cscale_unit,*) vel%cscale_vals(i,1),vel%cscale_vals(i,2),vel%cscale_vals(i,3)
        ENDDO
        CLOSE(vel%cscale_unit)
      ENDIF


      
      PRINT("(A)"), " "
      
      
      IF (zeta%plot_sol_option == 1) THEN   
        OPEN(unit=zeta%cscale_unit,file="zeta.cscale.out")  
        WRITE(zeta%cscale_unit,"(2I5)") snap_start-1,snap_end-1
      ENDIF
      IF (vel%plot_sol_option == 1) THEN     
        OPEN(unit=vel%cscale_unit,file="vel.cscale.out") 
        WRITE(vel%cscale_unit,"(2I5)") snap_start-1,snap_end-1
      ENDIF
      
      

      
      DO snap = snap_start,snap_end
      
        
        t_snap = t(snap)      
        PRINT("(A)"), "---------------------------------------------"
        PRINT("(A,I5,A,I5,A,F20.5)"), "Time snap: ", snap-1,"/",snap_end," t = ", t_snap        
!         WRITE(snap_char,"(I4.4)") snap-1


!         IF (zeta%plot_sol_option == 1 .or. vel%plot_sol_option == 1) THEN
!           PRINT("(A)"), "  Evaluating zeta at additional plotting points..."
!           CALL evaluate_depth_solution(ne,el_type,el_in,nplt,ndof,phi,snap,Z,zeta%sol_val,Zsnap_min,Zsnap_max)
!         ENDIF
!         
!         IF (zeta%plot_sol_option == 1) THEN        
!           PRINT("(A)"), "  Writing zeta PostScript file..."   
!           IF (zeta%cscale_option == "auto-snap") THEN
!             zeta%sol_min = Zsnap_min
!             zeta%sol_max = Zsnap_max            
!           ELSE IF (zeta%cscale_option == "auto-all") THEN
!             ! use existing values
!           ELSE IF (zeta%cscale_option == "file") THEN
! !             DO i = 1,num_cscale_zeta_vals
! !               IF (ABS(cscale_zeta_vals(i,1)-t_snap) < 1d-8) THEN
! !                 Z_min = cscale_zeta_vals(i,2)
! !                 Z_max = cscale_zeta_vals(i,3)              
! !                 EXIT
! !               ENDIF
! !             ENDDO            
! 
!             zeta%sol_min = zeta%cscale_vals(snap-1,2)
!             zeta%sol_max = zeta%cscale_vals(snap-1,3)              
!               
!           ELSE
!             zeta%sol_min = zeta%cscale_min
!             zeta%sol_max = zeta%cscale_max
!           ENDIF
!           
!           PRINT("(2(A,F10.5))"), "  zeta_min = ", zeta%sol_min, "  zeta_max = ", zeta%sol_max
!           WRITE(zeta%cscale_unit,"(3(e24.17,1x))") t_snap,zeta%sol_min,zeta%sol_max
!           
!           CALL latex_axes_labels(zeta%sol_min,zeta%sol_max,"surface elevation",t_snap,t_start,t_end)  
!           CALL run_latex()          
!           CALL write_psheader("zeta_"//snap_char//".ps",zeta%ps_unit)            
!           CALL plot_contours(zeta%ps_unit,nplt,ntri,rect,ne,el_type,el_in,xyplt,zeta%sol_val,zeta%sol_min,zeta%sol_max)      
!           IF (zeta%plot_mesh_option == 1) THEN
!             CALL plot_mesh(zeta%ps_unit,ne,nverts,el_type,el_in,xy,ect)
!           ENDIF      
!           CALL write_all_axes(zeta%ps_unit,"zeta",t_snap,t_start,t_end)               
!           CALL close_ps("zeta_"//snap_char,zeta%ps_unit)
!           CALL convert_ps("zeta_"//snap_char,frmt,density,rm_ps)
!         ENDIF        
      
        CALL make_plot(snap,t_snap,zeta,Z(:,:,snap))
        
      
        PRINT("(A)"), " "
        
      
      
      
!         IF (vel%plot_sol_option == 1) THEN
!           PRINT("(A)"), "  Evaluating velocity at additional plotting points..."
!           CALL evaluate_velocity_solution(ne,el_type,el_in,nplt,ndof,phi,snap,Qx,Qy,zeta%sol_val,bathy%sol_val,vel%sol_val,velsnap_min,velsnap_max)
!         ENDIF       
!         
!         IF (vel%plot_sol_option == 1) THEN    
!         
!           PRINT("(A)"), "  Writing velocity PostScript file..."   
!           IF (vel%cscale_option == "auto-snap") THEN
!             vel%sol_min = velsnap_min
!             vel%sol_max = velsnap_max            
!           ELSE IF (vel%cscale_option == "auto-all") THEN
!             ! use existing values
!           ELSE IF (vel%cscale_option == "file") THEN
! !             DO i = 1,num_cscale_vel_vals
! !               IF (ABS(cscale_vel_vals(i,1)-t_snap) < 1d-8) THEN
! !                 vel_min = cscale_vel_vals(i,2)
! !                 vel_max = cscale_vel_vals(i,3)              
! !                 EXIT
! !               ENDIF
! !             ENDDO        
! 
!            vel%sol_min = vel%cscale_vals(snap-1,2)
!            vel%sol_max = vel%cscale_vals(snap-1,3)
!             
!           ELSE
!             vel%sol_min = vel%cscale_min
!             vel%sol_max = vel%cscale_max
!           ENDIF
!           
!           PRINT("(2(A,F10.5))"), "  vel_min = ", vel%sol_min, "  vel_max = ", vel%sol_max
!           WRITE(vel%cscale_unit,"(3(e24.17,1x))") t_snap,vel%sol_min,vel%sol_max      
!           
!           CALL latex_axes_labels(vel%sol_min,vel%sol_max,"velocity",t_snap,t_start,t_end) 
!           CALL run_latex()          
!           CALL write_psheader("vel_"//snap_char//".ps",vel%ps_unit)        
!           CALL plot_contours(vel%ps_unit,nplt,ntri,rect,ne,el_type,el_in,xyplt,vel%sol_val,vel%sol_min,vel%sol_max)      
!           IF (vel%plot_mesh_option == 1) THEN
!             CALL plot_mesh(vel%ps_unit,ne,nverts,el_type,el_in,xy,ect)
!           ENDIF      
!           CALL write_all_axes(vel%ps_unit,"vel",t_snap,t_start,t_end)          
!           CALL close_ps("vel_"//snap_char,vel%ps_unit)   
!           CALL convert_ps("vel_"//snap_char,frmt,density,rm_ps)
!         ENDIF

        CALL make_plot(snap,t_snap,vel,zeta%sol_val,bathy%sol_val,Qx(:,:,snap),Qy(:,:,snap))
        
      ENDDO
      
      
!       WRITE(start_num,"(I3)") snap_start
!       WRITE(nframes,"(I3)") snap_end-snap_start + 1
      
      IF (zeta%plot_sol_option == 1 .and. make_movie == 1) THEN
        CALL SYSTEM("ffmpeg -i zeta_%04d."//frmt//" -y zeta.mp4")
      ENDIF          
      
      IF (vel%plot_sol_option == 1 .and. make_movie == 1) THEN
       CALL SYSTEM("ffmpeg -i vel_%04d."//frmt//" -y vel.mp4")
      ENDIF
      
  

      END PROGRAM plot_dg
      
