      SUBROUTINE read_stations()

      USE globals, ONLY: pres,nbou,fbseg,fbnds,xy,nsta,xysta,ndsta, &
                         out_direc

      IMPLICIT NONE
      
      INTEGER, PARAMETER :: nista = 4 ! number of stations bewteen 101 boundary nodes
      REAL(pres) :: r(nista),h
      INTEGER :: sta,i
      INTEGER :: bou,nd
      INTEGER :: nnds,btype
      LOGICAL :: file_exists
      
      h = 2d0/real(nista+1,pres)
      DO i = 1,nista
        r(i) = -1d0 + h*real(i,pres)
      ENDDO
      
      file_exists = .FALSE.
      INQUIRE(FILE=trim(out_direc) //"stations.d",EXIST=file_exists)
      
      IF(file_exists) THEN
      
        OPEN(UNIT=15,FILE=trim(out_direc) //"stations.d")
        READ(15,*) nsta
        ALLOCATE(xysta(2,nsta))
        DO sta = 1,nsta
          READ(15,*) xysta(1,sta), xysta(2,sta)
        ENDDO
        CLOSE(15)      
        
      ELSE
      
        ! Use the nodes along 101 type boundaries
      
        nsta = 0
      
        DO bou = 1,nbou
          nnds = fbseg(1,bou)
          btype = fbseg(2,bou)
          IF (btype == 101) THEN
            DO nd = 1,nnds
              nsta = nsta + 1            
            ENDDO
          ENDIF
        ENDDO       
      
        ALLOCATE(xysta(2,(nista+1)*nsta),ndsta(nista*nsta))
      
        nsta = 0
        
!         DO bou = 1,nbou
!           nnds = fbseg(1,bou)
!           btype = fbseg(2,bou)
!           IF (btype == 101) THEN
!             DO nd = 1,nnds
!               nsta = nsta + 1            
!               xysta(1,nsta) = xy(1,fbnds(nd,bou))
!               xysta(2,nsta) = xy(2,fbnds(nd,bou))                
!             ENDDO  
! !             EXIT ! only do first one
!           ENDIF
!         ENDDO          
      
        DO bou = 1,nbou
          nnds = fbseg(1,bou)
          btype = fbseg(2,bou)
          IF (btype == 101) THEN
            DO nd = 1,nnds-1
              nsta = nsta + 1            
              xysta(1,nsta) = xy(1,fbnds(nd,bou))
              xysta(2,nsta) = xy(2,fbnds(nd,bou))             
              DO i = 1,nista
                nsta = nsta + 1            
                xysta(1,nsta) = .5d0*((1d0-r(i))*xy(1,fbnds(nd,bou)) + (1d0+r(i))*xy(1,fbnds(nd+1,bou)))
                xysta(2,nsta) = .5d0*((1d0-r(i))*xy(2,fbnds(nd,bou)) + (1d0+r(i))*xy(2,fbnds(nd+1,bou)))  
              ENDDO
            ENDDO
            nsta = nsta + 1            
            xysta(1,nsta) = xy(1,fbnds(nnds,bou))
            xysta(2,nsta) = xy(2,fbnds(nnds,bou))     
            EXIT ! only do first one
          ENDIF
        ENDDO  
        
        
      
        OPEN(UNIT=15,FILE=trim(out_direc) //"stations.d")
        WRITE(15,*) nsta
        DO sta = 1,nsta
          WRITE(15,"(2(e24.17,2x))") xysta(1,sta), xysta(2,sta)
        ENDDO
        CLOSE(15)
      
      ENDIF
        
        
        
      RETURN
      END SUBROUTINE read_stations