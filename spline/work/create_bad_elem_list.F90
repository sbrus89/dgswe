      PROGRAM create_bad_elem_list
      
      IMPLICIT NONE
      
      INTEGER :: read_stat
      INTEGER :: nlines
      INTEGER :: el
      INTEGER :: i
      INTEGER :: nbad
      INTEGER :: bad_list(1000)
      INTEGER :: flag
      CHARACTER(100) :: line
      
      
      OPEN(UNIT=12,FILE='bad_elems.txt')
      READ(12,*) line
      
      nbad = 0
      nlines = 1
      bad_list = 0
      
      DO 
                
        IF (mod(nlines,3) == 1) THEN
          READ(12,*) el
          
          flag = 0
   check: DO i = 1,nbad
            IF (bad_list(i) == el) THEN
              flag = 1
              EXIT check
            ENDIF
          ENDDO check
          
          IF (flag == 0) THEN
            nbad = nbad + 1
            bad_list(nbad) = el
          ENDIF
          
        ELSE
          READ(12,*,IOSTAT=read_stat) line
        ENDIF               
              
        IF(read_stat < 0) THEN  
          EXIT 
        ENDIF
        
        nlines = nlines + 1
      
      ENDDO
      
      CLOSE(12)
      
      
      
      OPEN(UNIT=13,FILE='element.label')
      
      WRITE(13,*) nbad
      DO i = 1,nbad
        WRITE(13,*) bad_list(i)
      ENDDO
      
      
      CLOSE(13)
      
      
      END PROGRAM