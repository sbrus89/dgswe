      MODULE allocation

      USE globals, ONLY: grid,mnnds

      CONTAINS


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