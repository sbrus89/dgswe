      SUBROUTINE read_solution()

      USE globals, ONLY:pres,coarse,fine,lines

      IMPLICIT NONE
      
      INTEGER :: i,dof,el,line,et
      REAL(pres) :: tcoarse,tfine,hb
      
#ifndef adcirc      
      
      OPEN(UNIT=11, FILE=trim(coarse%out_direc) //"solution_H.d")
      OPEN(UNIT=12, FILE=trim(coarse%out_direc) //"solution_Qx.d")     
      OPEN(UNIT=13, FILE=trim(coarse%out_direc) //"solution_Qy.d")      
      
      READ(11,*) 
      READ(12,*)
      READ(13,*)
      
      OPEN(UNIT=21, FILE=trim(fine%out_direc) //"solution_H.d")
      OPEN(UNIT=22, FILE=trim(fine%out_direc) //"solution_Qx.d")     
      OPEN(UNIT=23, FILE=trim(fine%out_direc) //"solution_Qy.d")      
      
      READ(21,*) 
      READ(22,*)
      READ(23,*)
      

      
      DO line = 1,lines+1
        READ(11,*) tcoarse
        READ(12,*) tcoarse
        READ(13,*) tcoarse
        
        READ(21,*) tfine
        READ(22,*) tfine
        READ(23,*) tfine       
        
        DO dof = 1,coarse%mndof
          READ(11,*) (coarse%H(el,dof), el = 1,coarse%ne)
          READ(12,*) (coarse%Qx(el,dof), el = 1,coarse%ne)
          READ(13,*) (coarse%Qy(el,dof), el = 1,coarse%ne)
        ENDDO    
        
        DO dof = 1,fine%mndof
          READ(21,*) (fine%H(el,dof), el = 1,fine%ne)
          READ(22,*) (fine%Qx(el,dof), el = 1,fine%ne)
          READ(23,*) (fine%Qy(el,dof), el = 1,fine%ne)
        ENDDO         
        
      ENDDO
      
      CLOSE(11)
      CLOSE(12)
      CLOSE(13)
      
      CLOSE(21)
      CLOSE(22)
      CLOSE(23)  
      
      IF(tcoarse /= tfine) THEN
        PRINT("(A)"), "Warning: time snaps do not agree"
        PRINT("(A,F15.7)"), "tcoarse = ", tcoarse
        PRINT("(A,F15.7)"), "tfine = ", tfine
        PRINT*, " " 
      ENDIF
      
      
#else      
      
      OPEN(UNIT=11, FILE=trim(coarse%out_direc) //"DG.63")
      OPEN(UNIT=12, FILE=trim(coarse%out_direc) //"DG.64")     
      
      READ(11,*) 
      READ(11,*)
      READ(12,*)
      READ(12,*)
      
      OPEN(UNIT=21, FILE=trim(fine%out_direc) //"DG.63")
      OPEN(UNIT=22, FILE=trim(fine%out_direc) //"DG.64")     
      
      READ(21,*) 
      READ(21,*)
      READ(22,*)      
      READ(22,*)

      
      DO line = 1,lines
        READ(11,*) tcoarse
        READ(12,*) tcoarse
        
        READ(21,*) tfine
        READ(22,*) tfine    
        
        DO el = 1,coarse%ne
          et = coarse%el_type(el)
          DO dof = 1,coarse%ndof(et)
            READ(11,*) i,coarse%H(el,dof), hb
            READ(12,*) i,coarse%Qx(el,dof), coarse%Qy(el,dof)
          ENDDO    
        ENDDO
        
        DO el = 1,fine%ne
          et = fine%el_type(el)
          DO dof = 1,fine%ndof(et)
            READ(21,*) i,fine%H(el,dof), hb
            READ(22,*) i,fine%Qx(el,dof), fine%Qy(el,dof)
          ENDDO         
        ENDDO
        
      ENDDO
      
      CLOSE(11)
      CLOSE(12)
      
      CLOSE(21)
      CLOSE(22)
      
      IF(tcoarse /= tfine) THEN
        PRINT("(A)"), "Warning: time snaps do not agree"
      ENDIF      
      
#endif      

      RETURN 
      END SUBROUTINE read_solution
