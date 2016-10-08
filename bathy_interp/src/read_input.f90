      SUBROUTINE read_input()

      USE globals, ONLY: base,eval

      IMPLICIT NONE
      
      INTEGER, PARAMETER :: ninp = 10
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
              READ(temp,*) base%ctp
              PRINT("(A,I3)"), "base ctp = ", base%ctp               
            CASE (3)
              base%curve_file = TRIM(temp)
              PRINT("(A,A)"), "base curve file= ", base%curve_file     
            CASE (4)
              READ(temp,*) base%hbp
              PRINT("(A,I3)"), "base hbp = ", base%hbp               
            CASE (5)
              base%bathy_file = TRIM(temp)
              PRINT("(A,A)"), "base bathymetry file= ", base%bathy_file
            CASE (6)
              eval%grid_file = TRIM(temp)
              PRINT("(A,A)"), "eval grid_file = ", eval%grid_file
            CASE (7)
              READ(temp,*) eval%ctp
              PRINT("(A,I3)"), "eval ctp = ", eval%ctp                         
            CASE (8)
              eval%curve_file = TRIM(temp)
              PRINT("(A,A)"), "eval curve file= ", eval%curve_file                    
            CASE (9)
              READ(temp,*) eval%hbp
              PRINT("(A,I3)"), "eval hbp = ", eval%hbp           
            CASE (10)
              eval%out_direc = TRIM(temp)
              PRINT("(A,A)"), "output directory = ", eval%out_direc              
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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      SUBROUTINE write_input(file_unit)
      
      USE globals, ONLY: base,eval     
      
      IMPLICIT NONE
      
      INTEGER :: file_unit
      CHARACTER(40) :: sha1         
      
      WRITE(file_unit,"(A)") "" 
      WRITE(file_unit,"(A,A)")     "base grid file = ", base%grid_file
      WRITE(file_unit,"(A,I5)")    "base ctp = ", base%ctp              
      WRITE(file_unit,"(A,A)")     "base curve file = ", base%curve_file              
      WRITE(file_unit,"(A,I5)")    "base hbp = ", base%hbp
      WRITE(file_unit,"(A,A)")     "base bathy file = ", base%bathy_file              
      WRITE(file_unit,"(A,A)")     "eval grid file = ", eval%grid_file
      WRITE(file_unit,"(A,I5)")    "eval ctp = ", eval%ctp              
      WRITE(file_unit,"(A,A)")     "eval curve file = ", eval%curve_file              
      WRITE(file_unit,"(A,I5)")    "eval hbp = ", eval%hbp   
      WRITE(file_unit,"(A,A)")     "output directory = ", eval%out_direc  
      
      WRITE(file_unit,"(A)") ""       
      
      WRITE(file_unit,"(A)") "base grid file SHA: "//sha1(base%grid_file,"./")
      WRITE(file_unit,"(A)") "base curve file SHA: "//sha1(base%curve_file,"./")
      WRITE(file_unit,"(A)") "base bathy file SHA: "//sha1(base%bathy_file,"./")
      WRITE(file_unit,"(A)") "eval grid file SHA: "//sha1(eval%grid_file,"./")  
      WRITE(file_unit,"(A)") "eval curve file SHA: "//sha1(eval%curve_file,"./")   
                  
      WRITE(file_unit,"(A)") ""                 
      
      RETURN
      END SUBROUTINE write_input