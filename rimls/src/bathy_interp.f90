      SUBROUTINE bathymetry_eval_nodes(mesh)

      USE globals, ONLY: rp,grid,nel_type,nverts
      USE bathymetry_interp_mod, ONLY: shape_functions_eltypes_at_hbp,bathy_coordinates               

      IMPLICIT NONE
      
      TYPE(grid) :: mesh
      INTEGER :: space
      REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: psic

      space = 1
      CALL shape_functions_eltypes_at_hbp(space,nel_type,mesh%np,psic) 
      
      CALL bathy_coordinates(mesh%ne,mesh%nnds,nverts,mesh%el_type,mesh%elxy,psic,mesh%elxyh, &
                             mesh%depth,mesh%ect,mesh%elhb)


      RETURN
      END SUBROUTINE bathymetry_eval_nodes
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!          


      SUBROUTINE bathymetry_base_nodes(mesh)
      
      USE globals, ONLY: rp,grid,nel_type,nverts,nrpt
      USE bathymetry_interp_mod, ONLY: shape_functions_eltypes_at_hbp,bathy_coordinates      
      
      IMPLICIT NONE
      
      TYPE(grid) :: mesh
      INTEGER :: space
      INTEGER :: hbp
      INTEGER :: el,et
      CHARACTER(2) :: extract
      REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: psic
      REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: dpdr
      REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: dpds
      
      
      hbp = mesh%np(5)      
      
      IF (nrpt == 0) THEN
        space = 1
        extract = "IN"
        mesh%np(5) = hbp + 2
        mesh%np(6) = hbp + 1
        mesh%mninds = (mesh%np(6)-1)**2
      ELSE
        space = -1
        extract = "NO"
      ENDIF     
      

            
      
      CALL shape_functions_eltypes_at_hbp(space,nel_type,mesh%np,psic,dpdr,dpds,extract,mesh%nnds)      

      
      CALL bathy_coordinates(mesh%ne,mesh%nnds,nverts,mesh%el_type,mesh%elxy,psic,mesh%elxyh, &
                             mesh%depth,mesh%ect,mesh%elhb, &
                             dpdr,dpds,mesh%dhbdx,mesh%dhbdy,mesh%nhb)      
      
      
      
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