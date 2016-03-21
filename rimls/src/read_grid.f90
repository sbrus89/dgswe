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
      
      ALLOCATE(mesh%vct(4,mesh%ne))
      ALLOCATE(mesh%nelnds(mesh%ne))
      ALLOCATE(mesh%elxy(mnnds,mesh%ne,2))
      ALLOCATE(mesh%elhb(mnnds,mesh%ne))
      ALLOCATE(mesh%vxyn(mesh%nn))
      ALLOCATE(mesh%xyhv(1,mesh%nn,3))      
      
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
      
      
      ALLOCATE(vxy_temp(2,nn),vflag(nn))
      vflag = 0
      

      DO i = 1,ne
        nvert = nverts(mesh%el_type(i))
        DO j = 1,nvert
          mesh%vct(j,i) = mesh%ect(j,i)
        ENDDO
      ENDDO

      
      n = 0
      DO i = 1,ne
        nvert = nverts(mesh%el_type(i))
        DO j = 1,nvert
          nd = mesh%vct(j,i)
          IF (vflag(nd) == 0) THEN  
            n = n + 1            
            mesh%vxyn(n) = nd
            vxy_temp(1,n) = mesh%xy(1,nd)
            vxy_temp(2,n) = mesh%xy(2,nd)
            vflag(nd) = 1
          ENDIF
        ENDDO
      ENDDO         
      
      ALLOCATE(mesh%vxy(2,n))
      
      DO i = 1,n
        mesh%vxy(1,i) = vxy_temp(1,i)
        mesh%vxy(2,i) = vxy_temp(2,i)
      ENDDO
      
      mesh%mnelnds = maxval(mesh%nelnds)      
      
      DO i = 1,nn  
        mesh%xyhv(1,i,1) = mesh%xy(1,i)
        mesh%xyhv(1,i,2) = mesh%xy(2,i)
        mesh%xyhv(1,i,3) = mesh%depth(i)
      ENDDO
      

      CALL read_open_boundaries(mesh%nope,mesh%neta,mesh%obseg,mesh%obnds) 

      CALL read_flow_boundaries(mesh%nbou,mesh%nvel,mesh%fbseg,mesh%fbnds)

      
      CALL print_grid_info(mesh%grid_file,mesh%grid_name,mesh%ne,mesh%nn)      
      
      RETURN
      END SUBROUTINE read_grid