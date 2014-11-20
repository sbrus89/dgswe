      SUBROUTINE read_input()

      USE globals, ONLY: grid_file,forcing_file,p,ctp,dt,tf,dramp,cf,lines,out_direc,npart
      USE messenger2, ONLY: myrank,dirname,lname,finish,nthreads

      IMPLICIT NONE
      
      INTEGER, PARAMETER :: ninp = 11
      INTEGER :: inp_read,skipped
      CHARACTER(55) :: temp
      LOGICAL :: file_exists
      
!       OPEN(UNIT=16, FILE=dirname(1:lname)//'/'//'fort.16')      
      
      INQUIRE(FILE=dirname(1:lname)//'/'//'dgswe.inp', EXIST = file_exists)
      IF(file_exists == .FALSE.) THEN
        PRINT*, "dgswe.inp file does not exist"
        CALL finish()
      ENDIF      
      
      OPEN(UNIT=15,FILE=dirname(1:lname)//'/'//'dgswe.inp')      
      
      inp_read = 0
      skipped = 0
      DO WHILE (inp_read < ninp)
      
        READ(15,"(A55)") temp
                    
        IF ( INDEX(temp,"!") == 1 .or. INDEX(temp,"          ") == 1) THEN
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
              READ(temp,*) dt
            CASE (6)
              READ(temp,*) tf
              tf = tf*86400d0
            CASE (7)
              READ(temp,*) dramp
            CASE (8)
              READ(temp,*) cf
            CASE (9)
              READ(temp,*) lines
            CASE (10)
              out_direc = TRIM(temp)
            CASE (11)
              READ(temp,*) npart
          END SELECT
            
        ENDIF
      
      
      ENDDO
           
      CLOSE(15)
      
#ifdef openmp      
      IF (npart < nthreads) THEN
        npart = nthreads
      ENDIF  
#endif      
      
      IF (myrank == 0) THEN
      
        PRINT "(A)", "---------------------------------------------"
        PRINT "(A)", "             Input Information               "
        PRINT "(A)", "---------------------------------------------"
        PRINT "(A)", " "      
        PRINT*, "grid_file = ", grid_file    
        PRINT*, "forcing_file = ", forcing_file 
        PRINT*, "p = ", p              
        PRINT*, "ctp = ", ctp       
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
      
#ifdef rk22
        PRINT*, "RK22 timestepping"
#elif rk33
        PRINT*, "RK33 timestepping"
#else
        PRINT*, "Forward Euler timestepping"
#endif      
        PRINT*, " "

      ENDIF

      END SUBROUTINE  read_input