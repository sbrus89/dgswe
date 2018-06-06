      MODULE initialize
      
      USE globals, ONLY: rp
      USE plot_globals, ONLY: solution_type,plot_type

      IMPLICIT NONE

      CONTAINS
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE read_dgsweinp(s,input_path,substitute_path,replace_path,sub_path)
      
      USE globals, ONLY: deg2rad      
      USE read_dginp, ONLY: read_input,p,ctp,hbp,h0,sphi0,slam0,sta_opt,tf, &
                            grid_file,curve_file,bathy_file,stations_file,out_direc
      
      IMPLICIT NONE
      
      TYPE(solution_type), INTENT(INOUT) :: s
      CHARACTER(*), INTENT(IN) :: input_path
      INTEGER, INTENT(IN) :: substitute_path
      CHARACTER(*), INTENT(IN) :: replace_path
      CHARACTER(*), INTENT(IN) :: sub_path
      
      CALL read_input(0,input_path)
      
      IF (substitute_path == 1) THEN
        CALL substitute_partial_path(grid_file,replace_path,sub_path)
        CALL substitute_partial_path(curve_file,replace_path,sub_path)    
        CALL substitute_partial_path(stations_file,replace_path,sub_path)       
      ENDIF      
      
      
      sphi0 = sphi0*deg2rad  ! These are still being used in labels_mod.F90
      slam0 = slam0*deg2rad  ! to determine if axis labels should be in lat/lon
      
      s%p = p
      s%ctp = ctp 
      s%hbp = hbp
      s%tf = tf
      s%grid_file = grid_file
      s%curve_file = curve_file
      s%bathy_file = bathy_file
      s%out_direc = TRIM(ADJUSTL(input_path)) // TRIM(ADJUSTL(out_direc))
      s%stations_file = stations_file
      s%sta_opt = sta_opt
      s%sphi0 = sphi0
      s%slam0 = slam0
      s%h0 = h0
      
      RETURN      
      END SUBROUTINE
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE sizes(s)
      
      IMPLICIT NONE
      
      TYPE(solution_type), INTENT(INOUT) :: s      
      
      s%nel_type = 4
      
      s%ndof(1) = (s%p+1)*(s%p+2)/2
      s%ndof(2) = (s%p+1)**2
      s%ndof(3) = s%ndof(1)
      s%ndof(4) = s%ndof(2)      
      s%mndof = maxval(s%ndof)
      
      
      s%nverts(1) = 3
      s%nverts(2) = 4
      s%nverts(3) = 3
      s%nverts(4) = 4
      
      s%np(1) = 1
      s%np(2) = 1
      s%np(3) = s%ctp
      s%np(4) = s%ctp  
      s%mnp = maxval(s%np)+1

      s%nnds(1) = 3
      s%nnds(2) = 4
      s%nnds(3) = (s%ctp+1)*(s%ctp+2)/2
      s%nnds(4) = (s%ctp+1)*(s%ctp+1)      
      s%mnnds = maxval(s%nnds)       
      
      RETURN
      END SUBROUTINE sizes
      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE read_grid(s)
      
      USE grid_file_mod, ONLY: read_header,read_coords,read_connectivity,init_element_coordinates, &
                               read_open_boundaries,read_flow_boundaries,grid_size,read_curve_file, &
                               print_grid_info
      USE transformation, ONLY: element_area
      
      IMPLICIT NONE
      
      TYPE(solution_type), INTENT(INOUT) :: s
      LOGICAL :: cb_file_exists      
      
      PRINT("(A)"), s%grid_file
      CALL read_header(0,s%grid_file,s%grid_name,s%ne,s%nn)        
      CALL read_coords(s%nn,s%xy,s%depth)
      CALL read_connectivity(s%ne,s%ect,s%el_type) 
      CALL init_element_coordinates(s%ne,s%ctp,s%el_type,s%nverts,s%xy,s%ect,s%elxy)                  
      CALL read_open_boundaries(s%nope,s%neta,s%obseg,s%obnds)            
      CALL read_flow_boundaries(s%nbou,s%nvel,s%fbseg,s%fbnds)     
      CALL grid_size(s%ne,s%el_type,s%ect,s%xy,s%el_size)
      CALL read_curve_file(0,s%curve_file,s%ctp,s%nbou,s%xy,s%bndxy,cb_file_exists)  
      CALL print_grid_info(s%grid_file,s%grid_name,s%ne,s%nn)                 
      
      RETURN
      END SUBROUTINE read_grid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE connectivity(s)
      
      USE edge_connectivity_mod      
      
      IMPLICIT NONE
      
      TYPE(solution_type), INTENT(INOUT) :: s
      INTEGER :: nred
      
      CALL elements_per_node(s%ne,s%nn,s%nverts,s%el_type,s%ect,s%nepn,s%mnepn,s%epn)       
      CALL find_edge_pairs(s%ne,s%nverts,s%el_type,s%ect,s%nepn,s%epn,s%ned,s%ged2el,s%ged2nn,s%ged2led)      
      CALL find_interior_edges(s%ned,s%ged2el,s%nied,s%iedn,s%ed_type,s%recv_edge)      
      CALL find_open_edges(s%nope,s%obseg,s%obnds,s%ged2nn,s%nobed,s%obedn,s%ed_type,s%recv_edge)            
      CALL find_flow_edges(s%nbou,s%fbseg,s%fbnds,s%ged2nn,s%nnfbed,s%nfbedn,s%nfbednn,s%nfbed,s%fbedn,s%recv_edge,s%ed_type)     
      nred = 0
      CALL print_connect_info(s%mnepn,s%ned,s%nied,s%nobed,s%nfbed,s%nnfbed,nred)

      
      
      
      RETURN
      END SUBROUTINE connectivity
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE curvilinear(s)
      
      USE curvilinear_nodes_mod
      USE transformation     
      
      IMPLICIT NONE
      
      
      TYPE(solution_type), INTENT(INOUT) :: s            
      
      PRINT("(A)"), "Calculating curved boundary information..."
      CALL shape_functions_linear_at_ctp(s%nel_type,s%np,s%psiv)                   
      CALL eval_coordinates_curved(s%ctp,s%nnds,s%nverts,s%el_type,s%xy,s%ect,s%fbseg,s%fbnds, &
                                   s%nnfbed,s%nfbedn,s%nfbednn,s%ged2el,s%ged2led, &
                                   s%psiv,s%bndxy,s%elxy)  
      CALL element_area(s%ne,s%nel_type,s%np,s%el_type,s%elxy,s%el_area)     
      
      RETURN
      END SUBROUTINE curvilinear
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE find_output_type(sol)
      
      IMPLICIT NONE
      
      TYPE(solution_type), INTENT(INOUT) :: sol
      
      LOGICAL :: zeta_file_exists
      LOGICAL :: vel_file_exists
      LOGICAL :: bathy_file_exists
      INTEGER :: found
      
      found = 0
      
      ! Check for adcirc output
      INQUIRE(FILE=sol%out_direc//"fort.63", EXIST=zeta_file_exists)      
      INQUIRE(FILE=sol%out_direc//"fort.64", EXIST=vel_file_exists)
      IF (zeta_file_exists .or. vel_file_exists) THEN
        sol%output_type = "adcirc"
        found = 1
      ENDIF

      ! Check for DG-SWEM output
      INQUIRE(FILE=sol%out_direc//"DG.63", EXIST=zeta_file_exists)
      INQUIRE(FILE=sol%out_direc//"DG.64", EXIST=vel_file_exists)
      IF (zeta_file_exists .or. vel_file_exists) THEN
        sol%output_type = "dgswem"
        found = 1
      ENDIF

      ! Check for dgswe output
      INQUIRE(FILE=sol%out_direc//"Z.sol", EXIST=zeta_file_exists)
      INQUIRE(FILE=sol%out_direc//"Qx.sol", EXIST=vel_file_exists)
      INQUIRE(FILE=sol%out_direc//"hb.sol", EXIST=bathy_file_exists)      
      IF (zeta_file_exists .or. vel_file_exists .or. bathy_file_exists) THEN
        sol%output_type = "dgswe"
        found = 1
      ENDIF            
      
      ! Default to dgswe
      IF (found == 0) THEN
        sol%output_type = "dgswe"
      ENDIF
      
      RETURN
      END SUBROUTINE find_output_type
      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE read_solutions(zeta,vel,bathy,s)          
     
      USE read_write_output, ONLY: read_solution_full,read_fort6163,read_fort6264, &
                                   read_dg63,read_dg64
      USE grid_file_mod, ONLY: read_bathy_file
      USE bathymetry_interp_mod, ONLY: bathymetry_nodal2modal,dgswem_bathymetry_nodal2modal      
      
     
      IMPLICIT NONE
     
      TYPE(plot_type), INTENT(IN) :: zeta
      TYPE(plot_type), INTENT(IN) :: vel
      TYPE(plot_type), INTENT(IN) :: bathy
      TYPE(solution_type), INTENT(INOUT) :: s
     
     
      INTEGER :: el,nd,snap,nv
      REAL(rp) :: H
      LOGICAL :: file_exists
      
      REAL(rp), DIMENSION(:,:), ALLOCATABLE :: eta
      REAL(rp), DIMENSION(:,:), ALLOCATABLE :: uu2,vv2
      REAL(rp), DIMENSION(:,:), ALLOCATABLE :: hbm      
     
      IF (s%output_type == "adcirc") THEN       
      
        IF (zeta%plot_sol_option == 1 .or. vel%plot_sol_option == 1) THEN
          PRINT("(A)"), "Reading zeta solution..."          
          CALL read_fort6163(s%out_direc,"63",s%t,eta,s%nsnap_Z) 
          ALLOCATE(s%Z(4,s%ne,s%nsnap_Z))
          DO snap = 1,s%nsnap_Z
            DO el = 1,s%ne
              nv = s%nverts(s%el_type(el))            
              DO nd = 1,nv
                s%Z(nd,el,snap) = eta(s%ect(nd,el),snap)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        IF (bathy%plot_sol_option == 1 .or. vel%plot_sol_option == 1) THEN
          PRINT("(A)"), "Reading bathymetry solution..."              
          ALLOCATE(s%hb(4,s%ne,1))
          DO el = 1,s%ne
            nv = s%nverts(s%el_type(el))          
            DO nd = 1,nv
              s%hb(nd,el,1) = s%depth(s%ect(nd,el))
            ENDDO
          ENDDO
        ENDIF        
        IF (vel%plot_sol_option == 1) THEN
          PRINT("(A)"), "Reading u and v solutions..."       
          CALL read_fort6264(s%out_direc,"64",s%t,uu2,vv2,s%nsnap_Qx)
          s%nsnap_Qy = s%nsnap_Qx          
          ALLOCATE(s%Qx(4,s%ne,s%nsnap_Qx))
          ALLOCATE(s%Qy(4,s%ne,s%nsnap_Qx))        
          DO snap = 1,s%nsnap_Qx
            DO el = 1,s%ne
              nv = s%nverts(s%el_type(el))            
              DO nd = 1,nv
                H = s%Z(nd,el,snap)+s%hb(nd,el,1)
                s%Qx(nd,el,snap) = uu2(s%ect(nd,el),snap)*H
                s%Qy(nd,el,snap) = vv2(s%ect(nd,el),snap)*H
              ENDDO
            ENDDO
          ENDDO                
        ELSE 
          ALLOCATE(s%Qx(1,1,1),s%Qy(1,1,1))
        ENDIF             
      
      ELSE IF (s%output_type == "dgswem") THEN

        IF (zeta%plot_sol_option == 1 .or. vel%plot_sol_option == 1) THEN
          PRINT("(A)"), "Reading zeta solution..."          
          CALL read_DG63(s%out_direc,s%ne,s%t,s%Z,s%nsnap_Z) 
        ENDIF
        IF (vel%plot_sol_option == 1) THEN
          PRINT("(A)"), "Reading Qx and Qy solutions..."             
          CALL read_DG64(s%out_direc,s%ne,s%t,s%Qx,s%Qy,s%nsnap_Qx)        
          s%nsnap_Qy = s%nsnap_Qx
        ELSE 
          ALLOCATE(s%Qx(1,1,1),s%Qy(1,1,1))
        ENDIF  
        IF (bathy%plot_sol_option == 1 .or. vel%plot_sol_option == 1) THEN
          ALLOCATE(s%hb(3,s%ne,1))
          CALL dgswem_bathymetry_nodal2modal(s%ne,s%ect,s%depth,hbm)
          s%hb(:,:,1) = hbm(:,:)
        ENDIF  
     
      ELSE IF (s%output_type == "dgswe") THEN

        IF (zeta%plot_sol_option == 1 .or. vel%plot_sol_option == 1) THEN
          PRINT("(A)"), "Reading zeta solution..."          
          CALL read_solution_full(s%out_direc,"Z.sol","N",s%t,s%Z,s%nsnap_Z) 
        ENDIF
        IF (vel%plot_sol_option == 1) THEN
          PRINT("(A)"), "Reading Qx and Qy solutions..."       
          CALL read_solution_full(s%out_direc,"Qx.sol","N",s%t,s%Qx,s%nsnap_Qx)        
          CALL read_solution_full(s%out_direc,"Qy.sol","N",s%t,s%Qy,s%nsnap_Qy) 
        ELSE 
          ALLOCATE(s%Qx(1,1,1),s%Qy(1,1,1))
        ENDIF  
        IF (bathy%plot_sol_option == 1 .or. vel%plot_sol_option == 1) THEN
          INQUIRE(file=s%out_direc//"hb.sol",exist=file_exists)
          IF (file_exists .eqv. .true.) THEN
            PRINT("(A)"), "Reading bathymetry solution..."            
            CALL read_solution_full(s%out_direc,"hb.sol","N",s%t,s%hb,s%nsnap_hb)  
          ELSE
            CALL read_bathy_file(0,s%bathy_file,s%hbp,s%ne,s%el_type,s%nverts,s%depth,s%ect,s%elhb,file_exists)
            ALLOCATE(s%hb(s%mndof,s%ne,1))
            CALL bathymetry_nodal2modal(s%hbp,s%mndof,s%ne,s%el_type,s%elhb,hbm)
            s%hb(:,:,1) = hbm(:,:)
          ENDIF
        ENDIF    

     ENDIF           
     
     RETURN
     END SUBROUTINE read_solutions
     
     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     

      END MODULE initialize
