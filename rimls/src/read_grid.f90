      SUBROUTINE read_grid(mesh)
      
      USE globals, ONLY: rp,grid,ctp,mnnds,nverts, &
                         Erad,lambda0,phi0
      USE allocation, ONLY: grid_alloc    
      USE grid_file_mod
      USE spherical_mod, ONLY: cpp_transformation

      IMPLICIT NONE
      INTEGER :: i,j,k,el,n,nd
      INTEGER :: ne,nn,nbseg
      INTEGER :: btype
      INTEGER :: nvert
      INTEGER :: alloc_status   
      INTEGER, ALLOCATABLE, DIMENSION(:) :: vflag     
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: vxy_temp 
      TYPE(grid) :: mesh
                        
      INTEGER :: et,nv                        
      INTEGER :: coord_sys
      REAL(rp) :: h0
      
      coord_sys = 1
      h0 = 0d0
      
      CALL read_header(mesh%grid_file,mesh%grid_name,mesh%ne,mesh%nn)
      
      CALL read_coords(mesh%nn,mesh%xy,mesh%depth,h0)  
      
      CALL cpp_transformation(coord_sys,Erad,lambda0,phi0,mesh%nn,mesh%xy)      
      
      CALL read_connectivity(mesh%ne,mesh%ect,mesh%el_type)
      
      ALLOCATE(mesh%elxy(mnnds,mesh%ne,2))
      ALLOCATE(mesh%elhb(mnnds,mesh%ne))
      
      ne = mesh%ne
      nn = mesh%nn        

      DO el = 1,ne       
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
      
      RETURN
      END SUBROUTINE read_grid