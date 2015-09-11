      SUBROUTINE write_output(init)
      
      USE globals, ONLY: t,mndof,ne, &
                         Hwrite,Zwrite,Qxwrite,Qywrite
      USE messenger2, ONLY: myrank
      USE read_dginp, ONLY: out_direc,grid_file,dt,tf

      IMPLICIT NONE

      INTEGER :: dof,el
      LOGICAL :: init

      IF (init) THEN
      
        ! Initialize files and write initial condition      
      
        IF (myrank == 0) THEN
          PRINT "(A)", "---------------------------------------------"
          PRINT "(A)", "               Time Stepping                 "
          PRINT "(A)", "---------------------------------------------"
          PRINT "(A)", " "

          PRINT "(A,e12.4)", "Time step: ",dt
          PRINT "(A,e12.4)", "Final time: ",tf

          PRINT "(A)", " "
        ENDIF      

        OPEN(unit=63,file=trim(out_direc) // 'solution_H.d')
        OPEN(unit=641,file=trim(out_direc) // 'solution_Qx.d')
        OPEN(unit=642,file=trim(out_direc) // 'solution_Qy.d')

        WRITE(63,"(A)") grid_file
        WRITE(641,"(A)") grid_file
        WRITE(642,"(A)") grid_file

      ELSE

        IF(myrank == 0) THEN
          PRINT("(A,e15.8)"), 't = ', t
        ENDIF
                  
      ENDIF



      WRITE(63,"(e24.17)") t
      DO dof = 1,mndof
        WRITE(63,"(16000(e24.17,1x))") (Zwrite(el,dof)%ptr, el = 1,ne)
      ENDDO

      WRITE(641,"(e24.17)") t
      DO dof = 1,mndof
        WRITE(641,"(16000(e24.17,1x))") (Qxwrite(el,dof)%ptr, el = 1,ne)
      ENDDO

      WRITE(642,"(e24.17)") t
      DO dof = 1,mndof
        WRITE(642,"(16000(e24.17,1x))") (Qywrite(el,dof)%ptr, el = 1,ne)
      ENDDO     

      RETURN
      END SUBROUTINE