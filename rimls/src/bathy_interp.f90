      SUBROUTINE bathymetry_eval_nodes(mesh)

      USE globals, ONLY: rp,grid,nel_type,nverts
      USE bathymetry_interp_mod, ONLY: shape_functions_eltypes_at_hbp,bathy_coordinates               

      IMPLICIT NONE
      
      TYPE(grid) :: mesh
      INTEGER :: space
      INTEGER :: mnnds,el
      REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: psic

      space = 1
      CALL shape_functions_eltypes_at_hbp(space,nel_type,mesh%np,psic) 
      
      mnnds = MAXVAL(mesh%nnds)
!       ALLOCATE(mesh%elhb(mnnds,mesh%ne))
      ALLOCATE(mesh%elxyh(mnnds,mesh%ne,2))
      
      DO el = 1,mesh%ne
        CALL bathy_coordinates(el,mesh%nnds,nverts,mesh%el_type,mesh%elxy,psic,mesh%elxyh, &
                               mesh%depth,mesh%ect,mesh%elhb)
      ENDDO


      RETURN
      END SUBROUTINE bathymetry_eval_nodes
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!          


      SUBROUTINE bathymetry_base_nodes(mesh)
      
      USE globals, ONLY: rp,grid,nel_type,nverts,nrpt
      USE basis, ONLY: element_nodes
      USE shape_functions_mod, ONLY: shape_functions_area_eval        
      USE bathymetry_interp_mod, ONLY: shape_functions_eltypes_at_hbp,bathy_coordinates      
      
      IMPLICIT NONE
      
      TYPE(grid) :: mesh
      INTEGER :: space
      INTEGER :: hbp
      INTEGER :: el,et
      INTEGER :: extract
      INTEGER :: hbp_type      
      INTEGER :: mnnds,npts,nnd
      REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: psic
      REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: dpdr
      REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: dpds
      REAL(rp), DIMENSION(:), ALLOCATABLE :: r,s      
      
      
      hbp = mesh%np(5)      
      
      IF (nrpt == 0) THEN
        space = 1
        extract = 5 ! extract interior nodes
        mesh%np(5) = hbp + 2
        mesh%np(6) = hbp + 1
        mesh%mninds = (mesh%np(6)-1)**2
      ELSE
        space = -1
        extract = 0 ! return all nodes
        mesh%mninds = (hbp+1)**2
      ENDIF     
      

                  
      CALL shape_functions_eltypes_at_hbp(space,nel_type,mesh%np,psic,dpdr,dpds,extract,mesh%nnds)      

      mnnds = MAXVAL(mesh%nnds)      
!       ALLOCATE(mesh%elhb(mnnds,mesh%ne))  
      ALLOCATE(mesh%dhbdx(mnnds,mesh%ne),mesh%dhbdy(mnnds,mesh%ne))
      ALLOCATE(mesh%nhb(mnnds,mesh%ne,3))
      ALLOCATE(mesh%elxyh(mnnds,mesh%ne,2))
      ALLOCATE(r(mnnds),s(mnnds))       
      
      DO el = 1,mesh%ne
      
        IF (space < 0 ) THEN
          et = mesh%el_type(el)
          IF (mod(et,2) == 1) THEN
            hbp_type = 5    
          ELSE IF (mod(et,2) == 0) THEN
            hbp_type = 6
          ENDIF    
        
          CALL element_nodes(et,space,mesh%np(hbp_type),npts,r,s)       
          CALL shape_functions_area_eval(et,mesh%np(et),nnd,npts,r,s, &
               psic(:,:,et),dpdr(:,:,et),dpds(:,:,et))        
        ENDIF
      
        CALL bathy_coordinates(el,mesh%nnds,nverts,mesh%el_type,mesh%elxy,psic,mesh%elxyh, &
                               mesh%depth,mesh%ect,mesh%elhb, &
                               dpdr,dpds,mesh%dhbdx,mesh%dhbdy,mesh%nhb)      
                             
      ENDDO                             
      
      
      
      ALLOCATE(mesh%npts_interior(mesh%ne))
      mesh%tpts_interior = 0
      DO el = 1,mesh%ne
        et = mesh%el_type(el)
        
        IF (mod(et,2) == 1) THEN
          mesh%npts_interior(el) = mesh%nnds(5)
        ELSE IF (mod(et,2) == 0) THEN
          mesh%npts_interior(el) = mesh%nnds(6)
        ENDIF
        
        mesh%tpts_interior = mesh%tpts_interior + mesh%npts_interior(el)
      ENDDO      

        
!       !Cross product normals as a check
!       x1 = elxy(2,el,1) - elxy(1,el,1)
!       x2 = elxy(3,el,1) - elxy(1,el,1)
!       
!       y1 = elxy(2,el,2) - elxy(1,el,2)
!       y2 = elxy(3,el,2) - elxy(1,el,2)
!       
!       z1 = elhb(2,el) - elhb(1,el)
!       z2 = elhb(3,el) - elhb(1,el)
!       
!       nx =   y1*z2 - y2*z1
!       ny = -(x1*z2 - x2*z1)
!       nz =   x1*y2 - x2*y1
!       
!       nrm = sqrt(nx**2 + ny**2 + nz**2)
!       
!       nx = nx/nrm
!       ny = ny/nrm
!       nz = nz/nrm
!       
!       PRINT*, (nhb(i,el), i=1,3)
!       PRINT*, nx,ny,nz
!       PRINT*, " "         
      
      RETURN
      END SUBROUTINE bathymetry_base_nodes
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!          


      SUBROUTINE element_center_nodes(mesh)
      
      USE globals, ONLY: rp,grid,nel_type,nverts
      USE bathymetry_interp_mod, ONLY: shape_functions_eltypes_at_hbp,bathy_coordinates      
      
      IMPLICIT NONE
      
      TYPE(grid) :: mesh
      INTEGER :: space
      INTEGER :: hbp
      INTEGER :: el,et
      INTEGER :: np(6),nnds(6)
      INTEGER :: ext
      INTEGER :: mnnds
      REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: psic

      space = 0      
      
      np = mesh%np        
      np(5) = 3
      np(6) = 2
      ext = 5    ! extract interior node
      
      nnds = mesh%nnds
      
      CALL shape_functions_eltypes_at_hbp(space,nel_type,np,psic,ext=ext,nnds=nnds)      
      
      mnnds = MAXVAL(mesh%nnds)          
      ALLOCATE(mesh%elxy_center(mnnds,mesh%ne,2))  
      
      DO el = 1,mesh%ne
        CALL bathy_coordinates(el,nnds,nverts,mesh%el_type,mesh%elxy,psic,mesh%elxy_center)     
      ENDDO 
  
      
      RETURN
      END SUBROUTINE element_center_nodes      