      PROGRAM node_precision

      IMPLICIT NONE

      
      INTEGER, PARAMETER :: prec = kind(1d0)
      INTEGER :: i,j,pe,nd
      INTEGER :: gne,gnn,npe,tnn,mnn,lnn,ndof,lines,ne,nn
      INTEGER :: lname
      INTEGER, ALLOCATABLE, DIMENSION(:) :: lnd2gnd
      CHARACTER(100) :: grid_name
      CHARACTER(6) :: dirname
      REAL(prec) :: xmax,ymax,hmax
      REAL(prec) :: xdiff,ydiff,hdiff
      REAL(prec), ALLOCATABLE, DIMENSION(:,:) :: xyh_global      
      REAL(prec), ALLOCATABLE, DIMENSION(:,:) :: xyh_local      
      
      grid_name = TRIM(ADJUSTL("../grids/inlet1_dble.grd"))
      
      OPEN(UNIT=14,FILE=grid_name)
      READ(14,*)
      READ(14,*) gne,gnn
      
      ALLOCATE(xyh_global(3,gnn),xyh_local(3,gnn))
      
      DO i = 1,gnn
        READ(14,*) j, xyh_global(1,j), xyh_global(2,j), xyh_global(3,j)
      ENDDO      
      
      CLOSE(14)


      
      
      
      
      
      OPEN(UNIT=19,FILE='PE0000/fort.81')      
      READ(19,*) npe
      CLOSE(19)

      DO pe = 1,npe
      
        dirname = "PE0000"
        lname = 6
        
        WRITE(dirname(3:lname),"(I4.4)") pe-1
        
        PRINT*, dirname
        
        OPEN(UNIT=19,FILE=dirname(1:lname)//'/'//'fort.81',POSITION='rewind')
        
        READ(19,*)
        READ(19,*) tnn
        READ(19,*) mnn
        READ(19,*) lnn
        READ(19,*) ndof
        READ(19,*) lines
        
        IF(.not. ALLOCATED(lnd2gnd)) THEN
          ALLOCATE(lnd2gnd(mnn))
        ENDIF
        
        DO nd = 1,lnn
          READ(19,*) lnd2gnd(nd)
        ENDDO
        
        CLOSE(19)    
        
        
        
        
        
        
        
        OPEN(UNIT=14,FILE=dirname(1:lname)//'/'//'fort.14')
        READ(14,*)
        READ(14,*) ne,nn
        
        DO i = 1,nn
          READ(14,*) j, xyh_local(1,lnd2gnd(j)), xyh_local(2,lnd2gnd(j)), xyh_local(3,lnd2gnd(j))
        ENDDO      
      
        CLOSE(14)        
        
        
      ENDDO
      
      
      
      
      
      
      
      xmax = 0d0
      ymax = 0d0
      hmax = 0d0
      
      DO i = 1,gnn
      
        xdiff = ABS(xyh_local(1,i)-xyh_global(1,i))
        IF (xdiff > xmax) THEN
          xmax = xdiff
        ENDIF
        
        ydiff = ABS(xyh_local(2,i)-xyh_global(2,i))
        IF (ydiff > ymax) THEN
          ymax = ydiff
        ENDIF
        
        hdiff = ABS(xyh_local(3,i)-xyh_global(3,i))
        IF (hdiff > hmax) THEN
          hmax = hdiff
        ENDIF        
        
      ENDDO
      
      PRINT*, "xmax = ", xmax
      PRINT*, "ymax = ", ymax
      PRINT*, "hmax = ", hmax

      END PROGRAM node_precision
