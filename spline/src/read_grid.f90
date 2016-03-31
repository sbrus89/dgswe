      SUBROUTINE read_grid(mesh)
      
      USE globals, ONLY: rp,grid,nverts,mnnds, &
                         Erad,lambda0,phi0          
      USE grid_file_mod
      USE spherical_mod, ONLY: cpp_transformation      

      IMPLICIT NONE
      INTEGER :: el,j
      INTEGER :: nv,et
      INTEGER :: alloc_status   
      INTEGER :: coord_sys
      TYPE(grid) :: mesh
      
      coord_sys = 1      
      Erad = 1d0
      lambda0 = 0d0
      phi0 = 0d0
      
      CALL read_header(mesh%grid_file,mesh%grid_name,mesh%ne,mesh%nn)

      CALL read_coords(mesh%nn,mesh%xy,mesh%depth)        
      
      CALL cpp_transformation(coord_sys,Erad,lambda0,phi0,mesh%nn,mesh%xy)        
      
      CALL read_connectivity(mesh%ne,mesh%ect,mesh%el_type)      

      
      

      ALLOCATE(mesh%elxy(mnnds,mesh%ne,2))     
      
      DO el = 1,mesh%ne
        et = mesh%el_type(el)
        nv = nverts(et)
        
        DO j = 1,nv
          mesh%elxy(j,el,1) = mesh%xy(1,mesh%ect(j,el))
          mesh%elxy(j,el,2) = mesh%xy(2,mesh%ect(j,el))
        ENDDO      
        
      ENDDO
      

      

      CALL read_open_boundaries(mesh%nope,mesh%neta,mesh%obseg,mesh%obnds)

      CALL read_flow_boundaries(mesh%nbou,mesh%nvel,mesh%fbseg,mesh%fbnds) 
      
      CALL print_grid_info(mesh%grid_file,mesh%grid_name,mesh%ne,mesh%nn)    

      RETURN
      END SUBROUTINE read_grid