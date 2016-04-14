      SUBROUTINE interp_check()

      USE globals, ONLY: rp,base,eval

      IMPLICIT NONE
      
      INTEGER :: el,nd,et,nnd
      REAL(rp) :: diff,max_diff


      IF (base%ne == eval%ne .AND. base%hbp == eval%hbp) THEN
      
      ELSE
        RETURN
      ENDIF
      
      
      max_diff = 0d0
       
      DO el = 1,base%ne
        et = base%el_type(el)
      
        IF (mod(et,2) == 1) THEN   
          nnd = base%nnds(5)
        ELSE IF (mod(et,2) == 0) THEN
          nnd = base%nnds(6)          
        ENDIF              
        
        DO nd = 1,nnd
          
          diff = ABS(base%elhb(nd,el) - eval%elhb(nd,el))           
          
          IF (diff > max_diff) THEN
            max_diff = diff
          ENDIF            
            
        ENDDO        
        
      ENDDO      
      
      PRINT*, " "
      PRINT*, "Max difference = ", max_diff


      END SUBROUTINE interp_check