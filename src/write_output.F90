      SUBROUTINE write_output()
      
      USE globals, ONLY: t,Hwrite,Qxwrite,Qywrite,mndof,ne
      USE messenger2, ONLY: myrank

      IMPLICIT NONE

      INTEGER :: dof,el

      IF(myrank == 0) THEN
        PRINT("(A,e15.8)"), 't = ', t
      ENDIF

      WRITE(63,"(e24.17)") t
      DO dof = 1,mndof
        WRITE(63,"(16000(e24.17,1x))") (Hwrite(el,dof)%ptr, el = 1,ne)
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