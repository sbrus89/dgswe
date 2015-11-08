      SUBROUTINE read_input()

      USE globals, ONLY: base,eval,p,ctp,Erad,lambda0,phi0,deg2rad,refinement,r,sigma_n
                         

      IMPLICIT NONE
      
      INTEGER, PARAMETER :: ninp = 7
      INTEGER :: inp_read,skipped
      CHARACTER(100) :: temp
      
      OPEN(unit=15,file='rimls.inp')
      
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
              PRINT*, "base%grid_file = ", base%grid_file
            CASE (2)
              eval%grid_file = TRIM(temp)
              PRINT*, "eval%grid_file = ", eval%grid_file
            CASE (3)
              READ(temp,*) ctp
              PRINT*, "ctp = ", ctp
            CASE (4) 
              READ(temp,*) Erad
              PRINT*, "Erad = ", Erad
            CASE (5)
              READ(temp,*) lambda0,phi0
              PRINT*, "lambda0,phi0 = ", lambda0 , phi0
              lambda0 = lambda0*deg2rad
              phi0 = phi0*deg2rad
            CASE (6)
              READ(temp,*) r
              PRINT*, "r = ", r
            CASE (7)
              READ(temp,*) sigma_n
              PRINT*, "sigma_n = ", sigma_n
          END SELECT
            
        ENDIF
      
      
      ENDDO
      
      PRINT*, " "
      PRINT*, "Lines skipped: ", skipped
      PRINT*, " "
      
      CLOSE(15)
      
      IF (base%grid_file == eval%grid_file) THEN
        refinement = .false.
      ELSE 
        refinement = .true.
      ENDIF
      PRINT*, "Refinement grid = ", refinement
      

      END SUBROUTINE  read_input