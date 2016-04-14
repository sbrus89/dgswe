      SUBROUTINE read_input()

      USE globals, ONLY: base,eval

      IMPLICIT NONE
      
      INTEGER, PARAMETER :: ninp = 8
      INTEGER :: inp_read,skipped
      CHARACTER(100) :: temp
      LOGICAL :: file_exists
      
      INQUIRE(file='bathy.inp',exist=file_exists)
      IF (file_exists == .FALSE.) THEN
        PRINT*, "bathy.inp file does not exist"
        STOP
      ENDIF
      
      OPEN(unit=15,file='bathy.inp')
      
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
              PRINT("(A,A)"), "base grid_file = ", base%grid_file
            CASE (2)
              base%bathy_file = TRIM(temp)
              PRINT("(A,A)"), "base bathymetry file= ", base%bathy_file
            CASE (3)
              READ(temp,*) base%ctp
              PRINT("(A,I3)"), "base ctp = ", base%ctp
            CASE (4)
              READ(temp,*) base%hbp
              PRINT("(A,I3)"), "base hbp = ", base%hbp
              PRINT*, " "
            CASE (5)
              eval%grid_file = TRIM(temp)
              PRINT("(A,A)"), "eval grid_file = ", eval%grid_file
            CASE (6)
              eval%out_direc = TRIM(temp)
              PRINT("(A,A)"), "output directory = ", eval%out_direc
            CASE (7)
              READ(temp,*) eval%ctp
              PRINT("(A,I3)"), "eval ctp = ", eval%ctp
            CASE (8)
              READ(temp,*) eval%hbp
              PRINT("(A,I3)"), "eval hbp = ", eval%hbp                           
          END SELECT
            
        ENDIF
      
      
      ENDDO
      
      PRINT*, " "
      PRINT("(A,I5)"), "Lines skipped: ", skipped
      PRINT*, " "  
      
      
      CLOSE(15)
      
      eval%sol_name = "eval"
      base%sol_name = "base"
      

      END SUBROUTINE  read_input