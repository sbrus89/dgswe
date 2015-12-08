      MODULE allocation

      USE globals, ONLY: grid
      
      CONTAINS
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
      SUBROUTINE sizes(mesh)
      
      IMPLICIT NONE   
      
      INTEGER :: hbp,ctp
      TYPE(grid) :: mesh
      
      ctp = mesh%ctp      
      hbp = mesh%hbp

      mesh%nverts(1) = 3
      mesh%nverts(2) = 4
      mesh%nverts(3) = 3
      mesh%nverts(4) = 4          
      
      mesh%np(1) = 1
      mesh%np(2) = 1
      mesh%np(3) = ctp
      mesh%np(4) = ctp   
      mesh%np(5) = hbp
      mesh%np(6) = hbp       
      
      mesh%nnds(1) = 3
      mesh%nnds(2) = 4
      mesh%nnds(3) = (ctp+1)*(ctp+2)/2
      mesh%nnds(4) = (ctp+1)*(ctp+1) 
      mesh%nnds(5) = (hbp+1)*(hbp+2)/2
      mesh%nnds(6) = (hbp+1)*(hbp+1)      
      mesh%mnnds = maxval(mesh%nnds)         

      RETURN
      END SUBROUTINE sizes      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      

      SUBROUTINE grid_alloc(stage,mesh)

      IMPLICIT NONE
      
      INTEGER :: stage
      TYPE(grid) :: mesh
      INTEGER :: i      
      INTEGER :: n 
      INTEGER :: alloc_status(11)

      alloc_status = 0
      
      IF (stage == 1) THEN
      
        n = 11
        
        ALLOCATE(mesh%ect(mesh%mnnds,mesh%ne),STAT = alloc_status(1))
        ALLOCATE(mesh%vct(4,mesh%ne),STAT = alloc_status(2))
        ALLOCATE(mesh%xy(2,mesh%nn),STAT = alloc_status(3))
        ALLOCATE(mesh%depth(mesh%nn),STAT = alloc_status(4))
        ALLOCATE(mesh%nelnds(mesh%ne),STAT = alloc_status(5))
        ALLOCATE(mesh%el_type(mesh%ne),STAT = alloc_status(6))
        ALLOCATE(mesh%elxy(mesh%mnnds,mesh%ne,2),STAT = alloc_status(7))
        ALLOCATE(mesh%elhb(mesh%mnnds,mesh%ne),STAT = alloc_status(8))
        ALLOCATE(mesh%hbxy(3,mesh%mnnds*mesh%ne),STAT = alloc_status(9))        
        ALLOCATE(mesh%vxyn(mesh%nn),STAT = alloc_status(10))
        ALLOCATE(mesh%vxy(2,mesh%nn),STAT = alloc_status(11))
      
      ELSE IF (stage == 2) THEN
      
        n = 3     
        
        ALLOCATE(mesh%obseg(mesh%nope),STAT = alloc_status(1))
        ALLOCATE(mesh%obnds(mesh%neta,mesh%nope),STAT = alloc_status(2))
        ALLOCATE(mesh%bnd_flag(mesh%nn),STAT = alloc_status(3))        
      
      ELSE IF (stage == 3) THEN
      
        n = 2
        
        ALLOCATE(mesh%fbseg(2,mesh%nbou),STAT = alloc_status(1))
        ALLOCATE(mesh%fbnds(mesh%nvel,mesh%nbou),STAT = alloc_status(2))        
      
      ENDIF
      
      DO i = 1,n
        IF (alloc_status(i) /= 0) THEN
          PRINT*, "Allocation error: grid_alloc"
          PRINT*, "Stage = ", stage
        ENDIF
      ENDDO

      RETURN
      END SUBROUTINE grid_alloc
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 



      END MODULE allocation

