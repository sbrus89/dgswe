      SUBROUTINE read_input()

      USE globals, ONLY: base,eval,Erad,lambda0,phi0,deg2rad, &
                         refinement,r,sigma_n,out_direc,lsp,nrpt,eps,curve_file
                         

      IMPLICIT NONE
      
      INTEGER, PARAMETER :: ninp = 15
      INTEGER :: inp_read,skipped
      CHARACTER(100) :: temp
      LOGICAL :: file_exists
      
      INQUIRE(file='rimls.inp',exist=file_exists)
      IF (file_exists == .FALSE.) THEN
        PRINT*, "rimls.inp file does not exist"
        STOP
      ENDIF      
      
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
              base%grid_file = TRIM(ADJUSTL(temp))
              PRINT*, "base%grid_file = ", base%grid_file
            CASE (2) 
              READ(temp,*) base%hbp
              PRINT*, "hbp = ", base%hbp              
            CASE (3)
              READ(temp,*) base%ctp
              PRINT*, "ctp = ", base%ctp              
            CASE (4)
              eval%grid_file = TRIM(ADJUSTL(temp))
              PRINT*, "eval%grid_file = ", eval%grid_file
            CASE (5) 
              READ(temp,*) eval%hbp
              PRINT*, "hbp = ", eval%hbp
            CASE (6)
              READ(temp,*) eval%ctp
              PRINT*, "ctp = ", eval%ctp              
            CASE (7)
              curve_file = TRIM(ADJUSTL(temp))
              PRINT*, "curve_file = ", curve_file
            CASE (8)
              READ(temp,*) lsp
              PRINT*, "lsp = ", lsp              
            CASE (9) 
              READ(temp,*) Erad
              PRINT*, "Erad = ", Erad
            CASE (10)
              READ(temp,*) lambda0,phi0
              PRINT*, "lambda0,phi0 = ", lambda0 , phi0
              lambda0 = lambda0*deg2rad
              phi0 = phi0*deg2rad
            CASE (11)
              READ(temp,*) r
              PRINT*, "r = ", r
            CASE (12)
              READ(temp,*) sigma_n
              PRINT*, "sigma_n = ", sigma_n
            CASE (13)
              out_direc = TRIM(ADJUSTL(temp))
              PRINT*, "out_direc = ", out_direc
            CASE (14)
              READ(temp,*) nrpt
              PRINT*, "nrpt = ", nrpt
            CASE (15)
              READ(temp,*) eps
              PRINT*, "eps = ", eps              
          
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