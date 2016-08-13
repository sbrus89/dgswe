      SUBROUTINE read_grid(mesh)
      
      USE globals, ONLY: rp,grid,nverts
      USE grid_file_mod

      IMPLICIT NONE

      INTEGER :: el,j           
      INTEGER :: et,nv                        
      INTEGER :: coord_sys
      INTEGER :: bou,nd
      INTEGER :: myrank
      REAL(rp) :: h0
      INTEGER :: alloc_status   
      TYPE(grid) :: mesh      
  
      coord_sys = 1
      h0 = 0d0  
  
      CALL read_header(0,mesh%grid_file,mesh%grid_name,mesh%ne,mesh%nn)  
        
      CALL read_coords(mesh%nn,mesh%xy,mesh%depth,h0)  
      
      CALL read_connectivity(mesh%ne,mesh%ect,mesh%el_type)
            
      CALL init_element_coordinates(mesh%ne,mesh%ctp,mesh%el_type,nverts,mesh%xy,mesh%ect,mesh%elxy)
      
      CALL init_element_depth(mesh%ne,mesh%hbp,mesh%el_type,nverts,mesh%depth,mesh%ect,mesh%elhb)
      
      CALL read_open_boundaries(mesh%nope,mesh%neta,mesh%obseg,mesh%obnds) 

      CALL read_flow_boundaries(mesh%nbou,mesh%nvel,mesh%fbseg,mesh%fbnds)       
      
      CALL print_grid_info(mesh%grid_file,mesh%grid_name,mesh%ne,mesh%nn) 
      
      myrank = 0                 
      CALL read_curve_file(myrank,mesh%curve_file,mesh%ctp,mesh%nbou,mesh%xy,mesh%bndxy)       
      
      
      
      ALLOCATE(mesh%bnd_flag(mesh%nn))
      mesh%bnd_flag = 0
      
      DO bou = 1,mesh%nope
        DO nd = 1,mesh%obseg(bou)
          mesh%bnd_flag(nd) = 1
        ENDDO
      ENDDO
      
      DO bou = 1,mesh%nbou
        DO nd = 1,mesh%fbseg(1,bou)
          mesh%bnd_flag(nd) = 1
        ENDDO
      ENDDO      
      
      

      RETURN
      END SUBROUTINE read_grid
