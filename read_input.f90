      SUBROUTINE read_input()

      USE globals, ONLY: base,fine,p,ctp,Erad,lambda0,phi0,deg2rad,refinement
                         

      IMPLICIT NONE
      
      INTEGER, PARAMETER :: ninp = 6
      INTEGER :: inp_read,skipped
      CHARACTER(100) :: temp
      
      OPEN(unit=15,file='dgswe.inp')
      
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
              base%grid_file = TRIM(temp)
              PRINT*, "grid_file = ", base%grid_file
            CASE (2)
              fine%grid_file = TRIM(temp)
              PRINT*, "forcing_file = ", fine%grid_file
            CASE (3)
              READ(temp,*) p
              PRINT*, "p = ", p
            CASE (4)
              READ(temp,*) ctp
              PRINT*, "ctp = ", ctp
            CASE (5) 
              READ(temp,*) Erad
              PRINT*, "Erad = ", Erad
            CASE (6)
              READ(temp,*) lambda0,phi0
              PRINT*, "lambda0,phi0 = ", lambda0 , phi0
              lambda0 = lambda0*deg2rad
              phi0 = phi0*deg2rad
          END SELECT
            
        ENDIF
      
      
      ENDDO
      
      PRINT*, " "
      PRINT*, "Lines skipped: ", skipped
      PRINT*, " "
      
      CLOSE(15)
      
      IF (base%grid_file == fine%grid_file) THEN
        refinement = .false.
      ELSE 
        refinement = .true.
      ENDIF
      PRINT*, "Refinement grid = ", refinement
      

      END SUBROUTINE  read_input