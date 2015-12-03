      SUBROUTINE read_bathy()

      USE globals, ONLY:rp,base

      IMPLICIT NONE
      
      INTEGER :: i,j,el,nd
      INTEGER :: ne,hbp,nnd
      
      
      OPEN(UNIT=11, FILE=trim(base%bathy_file))    
      
      READ(11,"(2(I7,1x))") ne,hbp
      
      DO i = 1,ne
      
        READ(11,"(2(I7),1x,60(e24.17,1x))") el,nnd,(base%elhb(j,el), j = 1,nnd)
      
      ENDDO
      
      CLOSE(11) 


      RETURN 
      END SUBROUTINE read_bathy
