      MODULE read_solution
      
      USE plot_globals, ONLY: rp,plot_type,solution_type
      USE read_write_output, ONLY: read_solution_full,read_fort6163,read_fort6264, &
                                   read_dg63,read_dg64
      USE bathymetry_interp_mod, ONLY: bathymetry_nodal2modal,dgswem_bathymetry_nodal2modal
      USE grid_file_mod      
      
      IMPLICIT NONE
      
      CONTAINS
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
      SUBROUTINE read_solutions(bathy,zeta,vel,sol)
      
      IMPLICIT NONE
      
      TYPE(plot_type), INTENT(IN) :: bathy
      TYPE(plot_type), INTENT(IN) :: zeta
      TYPE(plot_type), INTENT(IN) :: vel     
      TYPE(solution_type), INTENT(INOUT) :: sol
      
      INTEGER :: snap
      INTEGER :: el
      INTEGER :: nd
      
      REAL(rp) :: H
      REAL(rp), DIMENSION(:,:), ALLOCATABLE :: eta
      REAL(rp), DIMENSION(:,:), ALLOCATABLE :: uu2,vv2
      REAL(rp), DIMENSION(:,:), ALLOCATABLE :: hbm
      
      LOGICAL :: file_exists
      

#ifdef adcirc        
      
      IF (zeta%plot_sol_option == 1 .or. vel%plot_sol_option == 1) THEN
        PRINT("(A)"), "Reading zeta solution..."          
        CALL read_fort6163(sol%out_direc,"63",sol%t,eta,sol%nsnap_Z) 
        ALLOCATE(sol%Z(3,sol%ne,sol%nsnap_Z))
        DO snap = 1,sol%nsnap_Z
          DO el = 1,sol%ne
            DO nd = 1,3
              sol%Z(nd,el,snap) = eta(sol%ect(nd,el),snap)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
      IF (bathy%plot_sol_option == 1 .or. vel%plot_sol_option == 1) THEN
        PRINT("(A)"), "Reading bathymetry solution..."       
        ALLOCATE(sol%hb(3,sol%ne,1))
        DO el = 1,sol%ne
          DO nd = 1,3
            sol%hb(nd,el,1) = sol%depth(sol%ect(nd,el))
          ENDDO
        ENDDO
      ENDIF        
      IF (vel%plot_sol_option == 1) THEN
        PRINT("(A)"), "Reading u and v solutions..."       
        CALL read_fort6264(sol%out_direc,"64",sol%t,uu2,vv2,sol%nsnap_Qx)
        ALLOCATE(sol%Qx(3,sol%ne,sol%nsnap_Qx))
        ALLOCATE(sol%Qy(3,sol%ne,sol%nsnap_Qx))        
        DO snap = 1,sol%nsnap_Qx
          DO el = 1,sol%ne
            DO nd = 1,3
              H = sol%Z(nd,el,snap)+sol%hb(nd,el,1)
              sol%Qx(nd,el,snap) = uu2(ect(nd,el),snap)*H
              sol%Qy(nd,el,snap) = vv2(ect(nd,el),snap)*H
            ENDDO
          ENDDO
        ENDDO                
      ELSE 
        ALLOCATE(sol%Qx(1,1,1),sol%Qy(1,1,1))
      ENDIF             
      
#elif dgswem

      IF (zeta%plot_sol_option == 1 .or. vel%plot_sol_option == 1) THEN
        PRINT("(A)"), "Reading zeta solution..."          
        CALL read_DG63(sol%out_direc,sol%ne,sol%t,sol%Z,sol%nsnap_Z) 
      ENDIF
      IF (vel%plot_sol_option == 1) THEN
        PRINT("(A)"), "Reading Qx and Qy solutions..."       
        CALL read_DG64(sol%out_direc,sol%ne,sol%t,sol%Qx,sol%Qy,sol%nsnap_Qx)        
      ELSE 
        ALLOCATE(sol%Qx(1,1,1),sol%Qy(1,1,1))
      ENDIF  
      IF (bathy%plot_sol_option == 1 .or. vel%plot_sol_option == 1) THEN
        ALLOCATE(sol%hb(3,sol%ne,1))
        CALL dgswem_bathymetry_nodal2modal(sol%ne,sol%ect,sol%depth,hbm)
        sol%hb(:,:,1) = hbm(:,:)
      ENDIF  

#else
 
      IF (zeta%plot_sol_option == 1 .or. vel%plot_sol_option == 1) THEN
        PRINT("(A)"), "Reading zeta solution..."          
        CALL read_solution_full(sol%out_direc,"Z.sol","N",sol%t,sol%Z,sol%nsnap_Z) 
      ENDIF
      IF (vel%plot_sol_option == 1) THEN
        PRINT("(A)"), "Reading Qx and Qy solutions..."       
        CALL read_solution_full(sol%out_direc,"Qx.sol","N",sol%t,sol%Qx,sol%nsnap_Qx)        
        CALL read_solution_full(sol%out_direc,"Qy.sol","N",sol%t,sol%Qy,sol%nsnap_Qy) 
      ELSE 
        ALLOCATE(sol%Qx(1,1,1),sol%Qy(1,1,1))
      ENDIF  
      IF (bathy%plot_sol_option == 1 .or. vel%plot_sol_option == 1) THEN
        INQUIRE(file="hb.sol",exist=file_exists)
        IF (file_exists == .true.) THEN
          PRINT("(A)"), "Reading bathymetry solution..."          
          CALL read_solution_full(sol%out_direc,"hb.sol","N",sol%t,sol%hb,sol%nsnap_hb)  
        ELSE
          CALL read_bathy_file(0,sol%bathy_file,sol%hbp,sol%ne,sol%el_type,sol%nverts,sol%depth,sol%ect,sol%elhb,file_exists)
          ALLOCATE(sol%hb(sol%mndof,sol%ne,1))
          CALL bathymetry_nodal2modal(sol%hbp,sol%mndof,sol%ne,sol%el_type,sol%elhb,hbm)
          sol%hb(:,:,1) = hbm(:,:)
        ENDIF
      ENDIF   
 
#endif

      END SUBROUTINE read_solutions
      
      END MODULE read_solution