      SUBROUTINE read_input()

      USE globals, ONLY: grid_file,forcing_file,out_direc, &
                         p,ctp,hbp,dt,tf,dramp,cf,lines,nblk,npart,mndof

      IMPLICIT NONE
      
      INTEGER, PARAMETER :: ninp = 12
      INTEGER :: inp_read,skipped
      CHARACTER(50) :: temp
      LOGICAL :: file_exists
      

      INQUIRE(FILE='dgswe.inp', EXIST = file_exists)
      IF(file_exists == .FALSE.) THEN
        PRINT*, "dgswe.inp file does not exist"
      ELSE
        OPEN(unit=15,file='dgswe.inp')
      ENDIF
      
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", "             Input Information               "
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", " "
      
      inp_read = 0
      skipped = 0
      DO WHILE (inp_read < ninp)
      
        READ(15,"(A50)") temp
                    
        IF ( INDEX(temp,"!") == 1 .or. INDEX(temp,"          ") == 1) THEN
!           PRINT*, "input ignored"
            skipped = skipped + 1
        ELSE
          inp_read = inp_read + 1
          SELECT CASE (inp_read)
            CASE (1)
              grid_file = TRIM(temp)
            CASE (2)
              forcing_file = TRIM(temp)
            CASE (3)
              READ(temp,*) p
            CASE (4) 
              READ(temp,*) ctp
            CASE (5) 
              READ(temp,*) hbp
            CASE (6)
              READ(temp,*) dt
            CASE (7)
              READ(temp,*) tf
            CASE (8)
              READ(temp,*) dramp
            CASE (9)
              READ(temp,*) cf
            CASE (10)
              READ(temp,*) lines
            CASE (11)
              out_direc = TRIM(temp)
            CASE (12)
              READ(temp,*) npart
          END SELECT
            
        ENDIF
      
      
      ENDDO
      
      mndof = (p+1)**2
      
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", "             Input Information               "
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", " "      
      PRINT*, "grid_file = ", grid_file    
      PRINT*, "forcing_file = ", forcing_file 
      PRINT*, "p = ", p              
      PRINT*, "ctp = ", ctp       
      PRINT*, "hbp = ", hbp             
      PRINT*, "dt = ", dt    
      PRINT*, "tf = ", tf        
      PRINT*, "dramp = ", dramp   
      PRINT*, "cf = ", cf  
      PRINT*, "lines = ", lines   
      PRINT*, "out_direc = ", out_direc    
      PRINT*, "npart = ", npart       
      
      PRINT*, " "
      PRINT*, "Lines skipped: ", skipped
      PRINT*, " "
      
      CLOSE(15)
      

      END SUBROUTINE  read_input