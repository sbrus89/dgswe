      SUBROUTINE read_stations()

      USE globals, ONLY: nbou,fbseg,fbnds,xy,nsta,xysta,ndsta

      IMPLICIT NONE
      
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
      
      ALLOCATE(xysta(2,nsta),ndsta(nsta))
      
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

      RETURN
      END SUBROUTINE read_stations