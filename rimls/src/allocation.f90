      MODULE allocation

      USE globals, ONLY: grid,mnnds

      CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
      SUBROUTINE sizes()
      
      USE globals, ONLY: p,ctp,hbp,ndof,mndof,np,nnds,mnnds,nverts,mninds,order
      
      IMPLICIT NONE   
      
      ndof(1) = (p+1)*(p+2)/2
      ndof(2) = (p+1)*(p+1)
      ndof(3) = ndof(1)
      ndof(4) = ndof(2)
      mndof = maxval(ndof)
      
      nverts(1) = 3
      nverts(2) = 4
      nverts(3) = 3
      nverts(4) = 4      
      
      np(1) = 1
      np(2) = 1
      np(3) = ctp
      np(4) = ctp     
      np(5) = hbp
      np(6) = hbp
      
      nnds(1) = 3
      nnds(2) = 4
      nnds(3) = (ctp+1)*(ctp+2)/2
      nnds(4) = (ctp+1)*(ctp+1) 
      nnds(5) = (hbp+1)*(hbp+2)/2
      nnds(6) = (hbp+1)*(hbp+1)
      mnnds = maxval(nnds)       
      
      order(1) = 1
      order(2) = 2
      order(3) = 3
      order(4) = 4
      order(5) = 5
      order(6) = 6
      order(7) = 5
      order(8) = 6           
      
      mninds = (hbp-1)*(hbp-1)

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
      INTEGER :: alloc_status(10)

      alloc_status = 0
      
      IF (stage == 1) THEN
      
        n = 10
        
        ALLOCATE(mesh%ect(mnnds,mesh%ne),    STAT = alloc_status(1))
        ALLOCATE(mesh%vct(4,mesh%ne),        STAT = alloc_status(2))
        ALLOCATE(mesh%xy(2,mesh%nn),         STAT = alloc_status(3))
        ALLOCATE(mesh%depth(mesh%nn),        STAT = alloc_status(4))
        ALLOCATE(mesh%nelnds(mesh%ne),       STAT = alloc_status(5))
        ALLOCATE(mesh%el_type(mesh%ne),      STAT = alloc_status(6))
        ALLOCATE(mesh%elxy(mnnds,mesh%ne,2), STAT = alloc_status(7))
        ALLOCATE(mesh%elhb(mnnds,mesh%ne),   STAT = alloc_status(8))
        ALLOCATE(mesh%vxyn(mesh%nn),         STAT = alloc_status(9))
        ALLOCATE(mesh%xyhv(1,mesh%nn,3),     STAT = alloc_status(10))
      
      ELSE IF (stage == 2) THEN
      
        n = 2
        
        ALLOCATE(mesh%obseg(mesh%nope),           STAT = alloc_status(1))
        ALLOCATE(mesh%obnds(mesh%neta,mesh%nope), STAT = alloc_status(2))
        
      ELSE IF (stage == 3) THEN
        
        n = 2
        
        ALLOCATE(mesh%fbseg(2,mesh%nbou),         STAT = alloc_status(1))
        ALLOCATE(mesh%fbnds(mesh%nvel,mesh%nbou), STAT = alloc_status(2))
        
      END IF
      
      DO i = 1,n
        IF (alloc_status(i) /=0) THEN
          PRINT*, "Allocation error: grid_alloc"
          PRINT*, "Stage = ", stage
        ENDIF
      ENDDO
      
      RETURN
      END SUBROUTINE grid_alloc
  
      END MODULE allocation