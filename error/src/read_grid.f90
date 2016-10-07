
      SUBROUTINE read_grids()
      
      USE globals, ONLY: coarse,fine,base      
      USE allocation, ONLY: sizes
      
      IMPLICIT NONE
      
      INTEGER :: i,j
      
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", "           Coarse Grid Information           "
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", " "      
      
      CALL sizes(coarse)
      CALL read_grid(coarse)      
      
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", "            Fine Grid Information            "
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", " "     
     
      CALL sizes(fine)
      CALL read_grid(fine)
      
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", "            Base Grid Information            "
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", " "     
      
      CALL sizes(base)
      CALL read_grid(base)      
      
      RETURN
      END SUBROUTINE read_grids
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE read_grid(mesh)
      
      USE globals, ONLY: rp,solution,nverts
      USE grid_file_mod
      USE spherical_mod, ONLY: cpp_transformation

      IMPLICIT NONE

                 
      INTEGER :: el,j           
      INTEGER :: et,nv      
      INTEGER :: coord_sys
      REAL(rp) :: h0
      INTEGER :: alloc_status   
      TYPE(solution) :: mesh   
      LOGICAL :: cb_file_exists
      
      coord_sys = 1
      h0 = 0d0
      
      CALL read_header(0,mesh%grid_file,mesh%grid_name,mesh%ne,mesh%nn)
      
      CALL read_coords(mesh%nn,mesh%xy,mesh%depth,h0)  
      
!       CALL cpp_transformation(coord_sys,Erad,lambda0,phi0,mesh%nn,mesh%xy)      
      
      CALL read_connectivity(mesh%ne,mesh%ect,mesh%el_type)
      
      CALL init_element_coordinates(mesh%ne,mesh%ctp,mesh%el_type,nverts,mesh%xy,mesh%ect,mesh%elxy)     

      CALL read_open_boundaries(mesh%nope,mesh%neta,mesh%obseg,mesh%obnds) 

      CALL read_flow_boundaries(mesh%nbou,mesh%nvel,mesh%fbseg,mesh%fbnds) 
      
      CALL read_curve_file(0,mesh%curve_file,mesh%ctp,mesh%nbou,mesh%xy,mesh%bndxy,cb_file_exists)       
      
      CALL print_grid_info(mesh%grid_file,mesh%grid_name,mesh%ne,mesh%nn)      
      



      RETURN
      END SUBROUTINE read_grid