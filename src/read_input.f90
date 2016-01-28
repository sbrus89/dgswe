      SUBROUTINE read_input()

      USE globals, ONLY: grid_file,forcing_file,p,ctp,dt,tf,dramp,cf,lines,out_direc,npart

      IMPLICIT NONE
      
      INTEGER, PARAMETER :: ninp = 11
      INTEGER :: inp_read,skipped
      CHARACTER(55) :: temp
      
      OPEN(unit=15,file='spline.inp')
      
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", "             Input Information               "
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", " "
      
      inp_read = 0
      skipped = 0
      DO WHILE (inp_read < ninp)
      
        READ(15,"(A55)") temp
                    
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
              forcing_file = TRIM(temp)
              PRINT*, "forcing_file = ", forcing_file
            CASE (3)
              READ(temp,*) p
              PRINT*, "p = ", p
            CASE (4)
              READ(temp,*) ctp
              PRINT*, "ctp = ", ctp
            CASE (5)
              READ(temp,*) dt
              PRINT*, "dt = ", dt
            CASE (6)
              READ(temp,*) tf
              tf = tf*86400d0
              PRINT*, "tf = ", tf
            CASE (7)
              READ(temp,*) dramp
              PRINT*, "dramp = ", dramp
            CASE (8)
              READ(temp,*) cf
              PRINT*, "cf = ", cf
            CASE (9)
              READ(temp,*) lines
              PRINT*, "lines = ", lines
            CASE (10)
              out_direc = TRIM(temp)
              PRINT*, "out_direc = ", out_direc
            CASE (11)
              READ(temp,*) npart
              PRINT*, "npart = ", npart
          END SELECT
            
        ENDIF
      
      
      ENDDO
      
      PRINT*, " "
      PRINT*, "Lines skipped: ", skipped
      PRINT*, " "
      
      CLOSE(15)
      

      END SUBROUTINE  read_input