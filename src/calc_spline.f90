      MODULE calc_spline
      
      
      
      CONTAINS
      
      SUBROUTINE cubic_spline(sig,n,a,b,c,d,dt)

      USE globals, ONLY: rp
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: n
      INTEGER :: i      
      
      REAL(rp), INTENT(IN) :: sig
      REAL(rp), INTENT(IN) , DIMENSION(n) :: a
      REAL(rp), INTENT(OUT), DIMENSION(n) :: b,c,d
      REAL(rp), INTENT(OUT) :: dt
      REAL(rp) :: mult
      REAL(rp), DIMENSION(n) :: Ml,Md,Mu,v      
      
      
      
      dt = 1d0/(real(n,rp)-1d0)



      IF (sig < 1d-12) THEN   ! No tension      
      
          ! Set up matrix 
          Ml(1) = 0d0
          Md(1) = 1d0
          Mu(1) = 0d0
          DO i = 2,n-1
            Ml(i) = dt
            Md(i) = 4d0*dt
            Mu(i) = dt
          ENDDO
          Ml(n) = 0d0
          Md(n) = 1d0
          Mu(n) = 0d0

          ! Set up RHS
          v(1) = 0d0
          DO i = 2,n-1
            v(i) = 3d0/dt*(a(i+1)-2d0*a(i)+a(i-1))
          ENDDO
          v(n) = 0d0
          
      ELSE   ! With tension

          ! Set up matrix 
          Ml(1) = 0d0
          Md(1) = 1d0
          Mu(1) = 0d0
          DO i = 2,n-1
            Ml(i) = (2d0/sig**2)*(1d0/dt-sig/sinh(sig*dt))
            Md(i) = (4d0/sig**2)*((sig*cosh(sig*dt))/sinh(sig*dt)-1d0/dt)
            Mu(i) = (2d0/sig**2)*(1d0/dt-sig/sinh(sig*dt))
          ENDDO
          Ml(n) = 0d0
          Md(n) = 1d0
          Mu(n) = 0d0

          ! Set up RHS
          v(1) = 0d0
          DO i = 2,n-1
            v(i) = 1d0/dt*(a(i+1)-2d0*a(i)+a(i-1))
          ENDDO
          v(n) = 0d0
          
      ENDIF
      
      
      

      ! Solve system for c coefficients, forward sweep
      DO i = 2,n
        mult = Ml(i)/Md(i-1)
        Md(i) = Md(i) - mult*Mu(i-1)
        v(i) = v(i) - mult*v(i-1)
      ENDDO

      ! Solve system for c coefficients, backward sweep
      c(n) = v(n)/Md(n)
      DO i = n-1,1,-1
        c(i) = (v(i) - Mu(i)*c(i+1))/Md(i)
      ENDDO

      ! Solve for other coefficients d and b
      DO i = 1,n-1
        d(i) = (c(i+1)-c(i))/(3d0*dt)
        b(i) = (a(i+1)-a(i))/dt - dt*(2d0*c(i)+c(i+1))/3d0
      ENDDO
      
      
      
      
      RETURN
      END SUBROUTINE cubic_spline
      
      END MODULE calc_spline