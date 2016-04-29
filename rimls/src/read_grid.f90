      SUBROUTINE read_grid(mesh)
      
      USE globals, ONLY: rp,grid,nverts, &
                         Erad,lambda0,phi0   
      USE grid_file_mod
      USE spherical_mod, ONLY: cpp_transformation

      IMPLICIT NONE

                 
      INTEGER :: el,j           
      INTEGER :: et,nv                        
      INTEGER :: coord_sys
      REAL(rp) :: h0
      INTEGER :: alloc_status   
      TYPE(grid) :: mesh      
      
      coord_sys = 1
      h0 = 0d0
      
      CALL read_header(mesh%grid_file,mesh%grid_name,mesh%ne,mesh%nn)
      
      CALL read_coords(mesh%nn,mesh%xy,mesh%depth,h0)  
      
      CALL cpp_transformation(coord_sys,Erad,lambda0,phi0,mesh%nn,mesh%xy)      
      
      CALL read_connectivity(mesh%ne,mesh%ect,mesh%el_type)
      
      CALL init_element_coordinates(mesh%ne,mesh%ctp,mesh%el_type,nverts,mesh%xy,mesh%ect,mesh%elxy)
      
      CALL init_element_depth(mesh%ne,mesh%hbp,mesh%el_type,nverts,mesh%depth,mesh%ect,mesh%elhb)
      
      CALL read_open_boundaries(mesh%nope,mesh%neta,mesh%obseg,mesh%obnds) 

      CALL read_flow_boundaries(mesh%nbou,mesh%nvel,mesh%fbseg,mesh%fbnds) 
      
      CALL print_grid_info(mesh%grid_file,mesh%grid_name,mesh%ne,mesh%nn)      
      
      
      CALL grid_size(mesh%ne,mesh%el_type,mesh%ect,mesh%xy,mesh%h)
      
      
      RETURN
      END SUBROUTINE read_grid