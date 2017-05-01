      MODULE stats
      
      USE globals, ONLY:rp

      CONTAINS
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        

      SUBROUTINE r_squared(n,y,f,r2)

      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: n
      REAL(rp), DIMENSION(:), INTENT(IN) :: y
      REAL(rp), DIMENSION(:), INTENT(IN) :: f
      REAL(rp) :: r2
      
      INTEGER :: i
      REAL(rp) :: ybar
      REAL(rp) :: sstot
      REAL(rp) :: ssres
      
      ybar = 0d0
      DO i = 1,n
        ybar = ybar + y(i)
      ENDDO
      ybar = ybar/real(n,rp)
      
      sstot = 0d0
      DO i = 1,n
        sstot = sstot + (y(i)-ybar)**2
      ENDDO

      ssres = 0d0
      DO i = 1,n
        ssres = ssres + (y(i)-f(i))**2
      ENDDO
      
      r2 = 1d0-ssres/sstot
      
      RETURN
      END SUBROUTINE r_squared

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE scatter_index()
      
      IMPLICIT NONE
      
      
      
      RETURN
      END SUBROUTINE scatter_index


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      
      END MODULE stats