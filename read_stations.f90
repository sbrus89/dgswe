      SUBROUTINE read_stations()

      USE globals, ONLY: nbou,fbseg,fbnds,xy,nsta,xysta,ndsta,elsta, &
                         out_direc

      IMPLICIT NONE
      
      INTEGER :: sta
      INTEGER :: bou,nd
      INTEGER :: nnds,btype

      ! For now just use the nodes along 101 type boundaries
      
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
      
      ALLOCATE(xysta(2,nsta),ndsta(nsta),elsta(nsta))
      
      nsta = 0
      
      DO bou = 1,nbou
        nnds = fbseg(1,bou)
        btype = fbseg(2,bou)
        IF (btype == 101) THEN
          DO nd = 1,nnds
            nsta = nsta + 1            
            xysta(1,nsta) = xy(1,fbnds(nd,bou))
            xysta(2,nsta) = xy(2,fbnds(nd,bou))    
            ndsta(nsta) = fbnds(nd,bou)
          ENDDO
        ENDIF
      ENDDO  
      
      OPEN(UNIT=15,FILE=trim(out_direc) //"stations.d")
      WRITE(15,*) nsta
      DO sta = 1,nsta
        WRITE(15,*) xysta(1,sta), xysta(2,sta)
      ENDDO
      CLOSE(15)
        
      RETURN
      END SUBROUTINE read_stations