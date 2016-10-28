      SUBROUTINE read_input()

      USE globals, ONLY: base,eval,Erad,lambda0,phi0,deg2rad,basis_opt, &
                         refinement,r,sigma_n,out_direc,lsp,nrpt,eps
                         

      IMPLICIT NONE
      
      INTEGER, PARAMETER :: ninp = 17
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
              base%curve_file = TRIM(ADJUSTL(temp))
              PRINT*, "curve_file = ", base%curve_file              
            CASE (5)
              eval%grid_file = TRIM(ADJUSTL(temp))
              PRINT*, "eval%grid_file = ", eval%grid_file
            CASE (6) 
              READ(temp,*) eval%hbp
              PRINT*, "hbp = ", eval%hbp
            CASE (7)
              READ(temp,*) eval%ctp
              PRINT*, "ctp = ", eval%ctp              
            CASE (8)
              eval%curve_file = TRIM(ADJUSTL(temp))
              PRINT*, "curve_file = ", eval%curve_file
            CASE (9)
              READ(temp,*) lsp
              PRINT*, "lsp = ", lsp    
            CASE (10)
              READ(temp,*) basis_opt
              PRINT*, "basis_opt = ", basis_opt               
            CASE (11) 
              READ(temp,*) Erad
              PRINT*, "Erad = ", Erad
            CASE (12)
              READ(temp,*) lambda0,phi0
              PRINT*, "lambda0,phi0 = ", lambda0 , phi0
              lambda0 = lambda0*deg2rad
              phi0 = phi0*deg2rad
            CASE (13)
              READ(temp,*) r
              PRINT*, "r = ", r
            CASE (14)
              READ(temp,*) sigma_n
              PRINT*, "sigma_n = ", sigma_n
            CASE (15)
              out_direc = TRIM(ADJUSTL(temp))
              PRINT*, "out_direc = ", out_direc
            CASE (16)
              READ(temp,*) nrpt
              PRINT*, "nrpt = ", nrpt
            CASE (17)
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
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      SUBROUTINE write_input(file_unit)
      
      USE globals, ONLY: base,eval,Erad,lambda0,phi0,deg2rad,basis_opt, &
                         refinement,r,sigma_n,out_direc,lsp,nrpt,eps      
      
      IMPLICIT NONE
      
      INTEGER :: file_unit
      CHARACTER(40) :: sha1         
      
      WRITE(file_unit,"(A)") "" 
      WRITE(file_unit,"(A,A)")     "base grid file = ", base%grid_file
      WRITE(file_unit,"(A,I5)")    "hbp = ", base%hbp              
      WRITE(file_unit,"(A,I5)")    "ctp = ", base%ctp           
      WRITE(file_unit,"(A,A)")     "curve file = ", base%curve_file      
      WRITE(file_unit,"(A,A)")     "eval grid file = ", eval%grid_file
      WRITE(file_unit,"(A,I5)")    "hbp = ", eval%hbp
      WRITE(file_unit,"(A,I5)")    "ctp = ", eval%ctp              
      WRITE(file_unit,"(A,A)")     "curve file = ", eval%curve_file
      WRITE(file_unit,"(A,I5)")    "lsp = ", lsp    
      WRITE(file_unit,"(A,I5)")    "basis_opt = ", basis_opt               
      WRITE(file_unit,"(A,F11.3)") "Erad = ", Erad
      WRITE(file_unit,"(A,2(F11.3))") "lambda0,phi0 = ", lambda0 , phi0
      WRITE(file_unit,"(A,F11.3)") "r = ", r
      WRITE(file_unit,"(A,F11.3)") "sigma_n = ", sigma_n
      WRITE(file_unit,"(A,A)")     "out_direc = ", out_direc
      WRITE(file_unit,"(A,I5)")    "nrpt = ", nrpt
      WRITE(file_unit,"(A,F11.3)") "eps = ", eps      
      
      WRITE(file_unit,"(A)") ""       
      
      WRITE(file_unit,"(A)") "base grid file SHA: "//sha1(base%grid_file,"./")
      WRITE(file_unit,"(A)") "eval grid file SHA: "//sha1(eval%grid_file,"./")   
      WRITE(file_unit,"(A)") "base curve file SHA: "//sha1(base%curve_file,"./")      
      WRITE(file_unit,"(A)") "eval curve file SHA: "//sha1(eval%curve_file,"./")
                  
      WRITE(file_unit,"(A)") ""                 
      
      RETURN
      END SUBROUTINE write_input