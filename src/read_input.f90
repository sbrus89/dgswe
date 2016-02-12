      SUBROUTINE read_input()

      USE globals, ONLY: base,eval,ctp,out_direc,theta_tol,deform_tol,sig

      IMPLICIT NONE
      
      INTEGER, PARAMETER :: ninp = 7
      INTEGER :: inp_read,skipped
      CHARACTER(100) :: temp
      
      OPEN(unit=15,file='spline.inp')
      
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", "             Input Information               "
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", " "
      
      inp_read = 0
      skipped = 0
      DO WHILE (inp_read < ninp)
      
        READ(15,"(A100)") temp
                    
        IF ( INDEX(temp,"!") == 1 .or. INDEX(temp,"          ") == 1) THEN
!           PRINT*, "input ignored"
            skipped = skipped + 1
        ELSE
          inp_read = inp_read + 1
          SELECT CASE (inp_read)
            CASE (1)
              base%grid_file = TRIM(ADJUSTL(temp))
              PRINT*, "base%grid_file = ", base%grid_file
            CASE (2)
              eval%grid_file = TRIM(ADJUSTL(temp))
              PRINT*, "eval%grid_file = ", eval%grid_file
            CASE (3)
              READ(temp,*) ctp
              PRINT*, "ctp = ", ctp
            CASE (4)
              READ(temp,*) sig
              PRINT*, "sig = ", sig
            CASE (5)
              READ(temp,*) theta_tol
              PRINT*, "theta_tol = ", theta_tol          
            CASE (6)
              READ(temp,*) deform_tol
              PRINT*, "deform_tol = ", deform_tol                        
            CASE (7)
              out_direc = TRIM(temp)
              PRINT*, "out_direc = ", out_direc
          END SELECT
            
        ENDIF
      
      
      ENDDO
      
      PRINT*, " "
      PRINT*, "Lines skipped: ", skipped
      PRINT*, " "
      
      CLOSE(15)
      

      END SUBROUTINE  read_input