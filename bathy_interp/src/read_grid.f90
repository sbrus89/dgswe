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
  
      CALL read_header(mesh%grid_file,mesh%grid_name,mesh%ne,mesh%nn)  
  

      
      CALL read_coords(mesh%nn,mesh%xy,mesh%depth,h0)  
      
      CALL read_connectivity(mesh%ne,mesh%ect,mesh%el_type)
      
      
      ALLOCATE(mesh%elxy(mesh%mnnds,mesh%ne,2))
      ALLOCATE(mesh%elhb(mesh%mnnds,mesh%ne))

      DO el = 1,mesh%ne       
        et = mesh%el_type(el)
        nv = nverts(et)
        DO j = 1,nv
          mesh%elxy(j,el,1) = mesh%xy(1,mesh%ect(j,el))
          mesh%elxy(j,el,2) = mesh%xy(2,mesh%ect(j,el))
          mesh%elhb(j,el)   = mesh%depth(mesh%ect(j,el))
        ENDDO      
        
      ENDDO      

      
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
