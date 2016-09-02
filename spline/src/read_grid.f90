      SUBROUTINE read_grid(mesh)
      
      USE globals, ONLY: rp,grid,nverts, &
                         Erad,lambda0,phi0          
      USE grid_file_mod
      USE spherical_mod, ONLY: cpp_transformation      

      IMPLICIT NONE
      INTEGER :: el,j
      INTEGER :: nv,et
      INTEGER :: alloc_status   
      INTEGER :: coord_sys
      TYPE(grid) :: mesh
      
!       coord_sys = 2      
!       Erad = 6378206.4d0 
!       lambda0 = 95d0
!       phi0 = 29d0

      coord_sys = 1      
      Erad = 1d0 
      lambda0 = 0d0
      phi0 = 0d0
      
      CALL read_header(0,mesh%grid_file,mesh%grid_name,mesh%ne,mesh%nn)

      CALL read_coords(mesh%nn,mesh%xy,mesh%depth)        
      
      CALL cpp_transformation(coord_sys,Erad,lambda0,phi0,mesh%nn,mesh%xy)        
      
      CALL read_connectivity(mesh%ne,mesh%ect,mesh%el_type)      

      CALL init_element_coordinates(mesh%ne,mesh%ctp,mesh%el_type,nverts,mesh%xy,mesh%ect,mesh%elxy)

      CALL read_open_boundaries(mesh%nope,mesh%neta,mesh%obseg,mesh%obnds)

      CALL read_flow_boundaries(mesh%nbou,mesh%nvel,mesh%fbseg,mesh%fbnds) 
      
      CALL print_grid_info(mesh%grid_file,mesh%grid_name,mesh%ne,mesh%nn)    
      
      
      
      
      
      ALLOCATE(mesh%elxy_spline((mesh%ctp+1)**2,mesh%ne,2))
      mesh%elxy_spline = mesh%elxy
      
      ALLOCATE(mesh%el_type_spline(mesh%ne))
      mesh%el_type_spline = mesh%el_type

      RETURN
      END SUBROUTINE read_grid