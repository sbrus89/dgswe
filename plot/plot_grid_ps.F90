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
      USE read_dginp, ONLY: read_input,out_direc,p,ctp,hbp,tf, &
                            grid_file,curve_file,cb_file_exists
      USE plot_mod, ONLY: read_colormap,write_psheader,close_ps, &
                          scale_coordinates,zoom_box, &
                          evaluate_depth_solution,evaluate_velocity_solution, &
                          plot_contours,plot_mesh, &
                          latex_all_labels,write_all_axes,close_tex,convert_ps
      USE triangulation, ONLY: reference_element_delaunay
      USE edge_connectivity_mod
      USE curvilinear_nodes_mod
      USE transformation
      USE shape_functions_mod
      
      IMPLICIT NONE
      
      INTEGER :: start_snap,end_snap
      
      space = 0  
      
      CALL read_plot_input()
      
      IF (plot_mesh_option == 0 .and. plot_zeta_option == 0 .and. &
          plot_bathy_option == 0 .and. plot_vel_option == 0) THEN
        
        PRINT("(A)"), "No plot options have been specified"
        STOP
      ENDIF    
      
      CALL read_input(0,input_path)
      
      CALL sizes()

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
      ALLOCATE(phi(mndof,mnpp,nel_type))
      ALLOCATE(phi_hb(mndof,mnpp,nel_type))      
      ALLOCATE(rect(3,3*mnpp,nel_type))
      DO et = 1,nel_type
        CALL element_nodes(et,space,pplt(et),n,r,s)
        CALL element_basis(et,p,ndf,n,r,s,phi(:,:,et))
        CALL element_basis(et,hbp,ndf,n,r,s,phi_hb(:,:,et))        
        CALL reference_element_delaunay(n,r,s,ntri(et),rect(:,:,et))
        
!         DO i = 1,3
!           PRINT "(*(I5))", (rect(i,j,et), j = 1,ntri(et))
!         ENDDO    
!         PRINT*, ""
        PRINT("(A,I4,A,I4)"), "  number of additional nodes/sub-triangles: ", n,"/",ntri(et)
      ENDDO      


      
      PRINT("(A)"), "Finding zoom box..."
      CALL zoom_box(ne,el_type,nplt,xyplt,xbox_min,xbox_max,ybox_min,ybox_max, &
                                                     xmin,xmax,ymin,ymax,el_in)
                                                     
!       PRINT*, xmin,xmax,ymin,ymax                                                     

      PRINT("(A)"), "Scaling coordinates..."
      CALL scale_coordinates(ne,nn,el_type,nverts,nplt,figure_width,xmin,xmax,ymin,ymax,xyplt,xy)

      
      
      CALL read_colormap(cmap_file)

      
      

      snap_start = snap_start + 1
      snap_end = snap_end + 1      
      
      nsnap_Z = snap_end
      nsnap_Qx = snap_end
      nsnap_Qy = snap_end
      nsnap_hb = 1
      
      t_start = 0d0
      t_end = tf*86400d0      
      
          
        
      IF (plot_zeta_option == 1 .or. plot_vel_option == 1) THEN
        PRINT("(A)"), "Reading zeta solution..."          
        CALL read_solution_full(out_direc,"Z.sol","N",t,Z,nsnap_Z) 
      ENDIF
      IF (plot_vel_option == 1) THEN
        PRINT("(A)"), "Reading Qx and Qy solutions..."       
        CALL read_solution_full(out_direc,"Qx.sol","N",t,Qx,nsnap_Qx)        
        CALL read_solution_full(out_direc,"Qy.sol","N",t,Qy,nsnap_Qy)  
      ENDIF  
      IF (plot_bathy_option == 1 .or. plot_vel_option == 1) THEN
        PRINT("(A)"), "Reading bathymetry solution..."          
        CALL read_solution_full(out_direc,"hb.sol","N",t,hb,nsnap_hb)  
      ENDIF  
      

      
      
      IF (plot_bathy_option == 1 .or. plot_vel_option == 1) THEN
        PRINT("(A)"), "Evaluating bathymetry at additional plotting points..."
        CALL evaluate_depth_solution(ne,el_type,el_in,nplt,ndof_hb,phi_hb,1,hb,hb_val,hb_min,hb_max)
      ENDIF      
      
      
      
      IF (plot_bathy_option == 1) THEN
        PRINT("(A)"), "Writing bathymetry PostScript file..."
        CALL latex_all_labels(hb_min,hb_max,"bathymetry")         
        CALL write_psheader("bathy.ps",hb_unit)              
        CALL plot_contours(hb_unit,nplt,ntri,rect,ne,el_type,el_in,xyplt,hb_val,hb_min,hb_max)     
        IF (plot_bathy_mesh == 1) THEN
          CALL plot_mesh(hb_unit,ne,nverts,el_type,el_in,xy,ect)   
        ENDIF
        CALL write_all_axes(hb_unit)             
        CALL close_ps("bathy",hb_unit)
        CALL convert_ps("bathy",frmt,density,rm_ps)
      ENDIF
      
      IF (plot_mesh_option == 1) THEN          
        PRINT("(A)"), "Writing mesh PostScript file..."   
        CALL latex_all_labels()            
        CALL write_psheader("mesh.ps",mesh_unit)            
        CALL plot_mesh(mesh_unit,ne,nverts,el_type,el_in,xy,ect)   
        CALL write_all_axes(mesh_unit)          
        CALL close_ps("mesh",mesh_unit)         
        CALL convert_ps("mesh",frmt,density,0)
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
      
      
      

      Z_max = -1d10
      Z_min = 1d10
      vel_max = -1d10
      vel_min = 1d10
        
      DO snap = snap_start,snap_end
      
        IF ((plot_zeta_option == 1 .AND. cscale_zeta == "auto-all") .OR.  &
            (plot_vel_option == 1 .AND. cscale_vel == "auto-all")) THEN
              
          CALL evaluate_depth_solution(ne,el_type,el_in,nplt,ndof,phi,snap,Z,Z_val,Zsnap_min,Zsnap_max)
        
          IF (Zsnap_min < Z_min) THEN
            Z_min = Zsnap_min
          ENDIF
          
          IF (Zsnap_max > Z_max) THEN
            Z_max = Zsnap_max
          ENDIF
        ENDIF
          

        IF (plot_vel_option == 1 .AND. cscale_vel == "auto-all") THEN
          CALL evaluate_velocity_solution(ne,el_type,el_in,nplt,ndof,phi,snap,Qx,Qy,Z_val,hb_val,vel_val,velsnap_min,velsnap_max)
        
          IF (velsnap_min < vel_min) THEN
            vel_min = velsnap_min          
          ENDIF
         
          IF (velsnap_max > vel_max) THEN
            vel_max = velsnap_max
          ENDIF
        ENDIF
      ENDDO
      
      
      IF (cscale_zeta == "file") THEN
        READ(cscale_zeta_unit,*) start_snap,end_snap
        num_cscale_zeta_vals = end_snap - start_snap + 1
        ALLOCATE(cscale_zeta_vals(num_cscale_zeta_vals,3))
        DO i = 1,num_cscale_zeta_vals
          READ(cscale_zeta_unit,*) cscale_zeta_vals(i,1),cscale_zeta_vals(i,2),cscale_zeta_vals(i,3)
        ENDDO
        CLOSE(cscale_zeta_unit)
      ENDIF
      
      IF (cscale_vel == "file") THEN
        READ(cscale_vel_unit,*) start_snap,end_snap
        num_cscale_vel_vals = end_snap - start_snap + 1
        ALLOCATE(cscale_vel_vals(num_cscale_vel_vals,3))
        DO i = 1,num_cscale_vel_vals
          READ(cscale_vel_unit,*) cscale_vel_vals(i,1),cscale_vel_vals(i,2),cscale_vel_vals(i,3)
        ENDDO
        CLOSE(cscale_vel_unit)
      ENDIF


      
      PRINT("(A)"), " "
      
      
      IF (plot_zeta_option == 1) THEN   
        OPEN(unit=cscale_zeta_unit,file="zeta.cscale.out")  
        WRITE(cscale_zeta_unit,"(2I5)") snap_start-1,snap_end-1
      ENDIF
      IF (plot_vel_option == 1) THEN     
        OPEN(unit=cscale_vel_unit,file="vel.cscale.out") 
        WRITE(cscale_vel_unit,"(2I5)") snap_start-1,snap_end-1
      ENDIF
      
      
      snap_char = "0000"
      
      DO snap = snap_start,snap_end
      
        
        t_snap = t(snap)      
        PRINT("(A)"), "---------------------------------------------"
        PRINT("(A,I5,A,I5,A,F20.5)"), "Time snap: ", snap-1,"/",snap_end," t = ", t_snap        
        WRITE(snap_char,"(I4.4)") snap-1


        IF (plot_zeta_option == 1 .or. plot_vel_option == 1) THEN
          PRINT("(A)"), "  Evaluating zeta at additional plotting points..."
          CALL evaluate_depth_solution(ne,el_type,el_in,nplt,ndof,phi,snap,Z,Z_val,Zsnap_min,Zsnap_max)
        ENDIF
        
        IF (plot_zeta_option == 1) THEN        
          PRINT("(A)"), "  Writing zeta PostScript file..."   
          IF (cscale_zeta == "auto-snap") THEN
            Z_min = Zsnap_min
            Z_max = Zsnap_max            
          ELSE IF (cscale_zeta == "auto-all") THEN
            ! use existing values
          ELSE IF (cscale_zeta == "file") THEN
!             DO i = 1,num_cscale_zeta_vals
!               IF (ABS(cscale_zeta_vals(i,1)-t_snap) < 1d-8) THEN
!                 Z_min = cscale_zeta_vals(i,2)
!                 Z_max = cscale_zeta_vals(i,3)              
!                 EXIT
!               ENDIF
!             ENDDO            

            Z_min = cscale_zeta_vals(snap-1,2)
            Z_max = cscale_zeta_vals(snap-1,3)              
              
          ELSE
            Z_min = cscale_zeta_min
            Z_max = cscale_zeta_max
          ENDIF
          
          PRINT("(2(A,F10.5))"), "  Z_min = ", Z_min, "  Z_max = ", Z_max
          WRITE(cscale_zeta_unit,"(3(e24.17,1x))") t_snap,Z_min,Z_max
          
          CALL latex_all_labels(Z_min,Z_max,"surface elevation",t_snap,t_start,t_end)          
          CALL write_psheader("zeta_"//snap_char//".ps",Z_unit)            
          CALL plot_contours(Z_unit,nplt,ntri,rect,ne,el_type,el_in,xyplt,Z_val,Z_min,Z_max)      
          IF (plot_zeta_mesh == 1) THEN
            CALL plot_mesh(Z_unit,ne,nverts,el_type,el_in,xy,ect)
          ENDIF      
          CALL write_all_axes(Z_unit,t_snap,t_start,t_end)               
          CALL close_ps("zeta_"//snap_char,Z_unit)
          CALL convert_ps("zeta_"//snap_char,frmt,density,rm_ps)
        ENDIF        
      
      
      
        PRINT("(A)"), " "
        
      
      
      
        IF (plot_vel_option == 1) THEN
          PRINT("(A)"), "  Evaluating velocity at additional plotting points..."
          CALL evaluate_velocity_solution(ne,el_type,el_in,nplt,ndof,phi,snap,Qx,Qy,Z_val,hb_val,vel_val,velsnap_min,velsnap_max)
        ENDIF       
        
        IF (plot_vel_option == 1) THEN    
        
          PRINT("(A)"), "  Writing velocity PostScript file..."   
          IF (cscale_vel == "auto-snap") THEN
            vel_min = velsnap_min
            vel_max = velsnap_max            
          ELSE IF (cscale_vel == "auto-all") THEN
            ! use existing values
          ELSE IF (cscale_vel == "file") THEN
!             DO i = 1,num_cscale_vel_vals
!               IF (ABS(cscale_vel_vals(i,1)-t_snap) < 1d-8) THEN
!                 vel_min = cscale_vel_vals(i,2)
!                 vel_max = cscale_vel_vals(i,3)              
!                 EXIT
!               ENDIF
!             ENDDO        

           vel_min = cscale_vel_vals(snap-1,2)
           vel_max = cscale_vel_vals(snap-1,3)
            
          ELSE
            vel_min = cscale_vel_min
            vel_max = cscale_vel_max
          ENDIF
          
          PRINT("(2(A,F10.5))"), "  vel_min = ", vel_min, "  vel_max = ", vel_max
          WRITE(cscale_vel_unit,"(3(e24.17,1x))") t_snap,vel_min,vel_max      
          
          CALL latex_all_labels(vel_min,vel_max,"velocity",t_snap,t_start,t_end)                 
          CALL write_psheader("vel_"//snap_char//".ps",vel_unit)        
          CALL plot_contours(vel_unit,nplt,ntri,rect,ne,el_type,el_in,xyplt,vel_val,vel_min,vel_max)      
          IF (plot_vel_mesh == 1) THEN
            CALL plot_mesh(vel_unit,ne,nverts,el_type,el_in,xy,ect)
          ENDIF      
          CALL write_all_axes(vel_unit,t_snap,t_start,t_end)          
          CALL close_ps("vel_"//snap_char,vel_unit)   
          CALL convert_ps("vel_"//snap_char,frmt,density,rm_ps)
        ENDIF
            
      ENDDO
      
!       WRITE(start_num,"(I3)") snap_start
!       WRITE(nframes,"(I3)") snap_end-snap_start + 1
      
      IF (plot_zeta_option == 1 .and. make_movie == 1) THEN
        CALL SYSTEM("ffmpeg -i zeta_%04d."//frmt//" -y zeta.mp4")
      ENDIF          
      
      IF (plot_vel_option == 1 .and. make_movie == 1) THEN
       CALL SYSTEM("ffmpeg -i vel_%04d."//frmt//" -y vel.mp4")
      ENDIF
      
  

      END PROGRAM plot_grid_ps
      
