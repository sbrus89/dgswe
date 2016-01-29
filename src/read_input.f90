      SUBROUTINE read_input()

      USE globals, ONLY: grid_file,ctp,out_direc

      IMPLICIT NONE
      
      INTEGER, PARAMETER :: ninp = 3
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
              grid_file = TRIM(temp)
              PRINT*, "grid_file = ", grid_file
            CASE (2)
              READ(temp,*) ctp
              PRINT*, "ctp = ", ctp
            CASE (3)
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