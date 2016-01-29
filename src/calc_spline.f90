      MODULE calc_spline
      
      USE globals, ONLY: rp      
      
      IMPLICIT NONE
      
      
      
      CONTAINS
      
      SUBROUTINE calc_cubic_spline(sig,n,a,b,c,d,dt)


      
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
      END SUBROUTINE calc_cubic_spline
                  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   

      SUBROUTINE eval_cubic_spline(t,ti,a,b,c,d,f,fp,fpp)
      
      IMPLICIT NONE
      
      REAL(rp), INTENT(IN) :: t,ti
      REAL(rp), INTENT(IN) :: a,b,c,d
      REAL(rp), INTENT(OUT) :: f
      REAL(rp), INTENT(OUT), OPTIONAL :: fp,fpp
      
      f = a + b*(t-ti) + c*(t-ti)**2 + d*(t-ti)**3
      
      IF (PRESENT(fp)) THEN
        fp = b + 2d0*c*(t-ti) + 3d0*d*(t-ti)**2
      ENDIF 
      
      IF (PRESENT(fpp)) THEN
        fpp = 2d0*c + 6d0*d*(t-ti)
      ENDIF
      
      
      
      RETURN
      END SUBROUTINE eval_cubic_spline


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      SUBROUTINE newton(t,ti,xm,ym,ax,bx,cx,dx,ay,by,cy,dy,x,y)
      
      IMPLICIT NONE      
      
      INTEGER :: it,maxit
      REAL(rp), INTENT(INOUT) :: t
      REAL(rp), INTENT(IN) :: ti
      REAL(rp), INTENT(IN) :: xm,ym
      REAL(rp), INTENT(IN) :: ax,bx,cx,dx,ay,by,cy,dy
      REAL(rp), INTENT(OUT) :: x,y
      REAL(rp) :: f,fp,fpp,g,gp,gpp
      REAL(rp) :: d,dp,tol
      
      tol = 1d-4
      maxit = 10000
      
iter: DO it = 1,maxit

        CALL eval_cubic_spline(t,ti,ax,bx,cx,dx,f,fp,fpp)               
        CALL eval_cubic_spline(t,ti,ay,by,cy,dy,g,gp,gpp)
        
        d = 2d0*(f-xm)*fp + 2d0*(g-ym)*gp
        dp = 2d0*fp**2 + 2d0*(f-xm)*fpp + 2d0*gp**2 + 2d0*(g-ym)*gpp
        
        t = t - d/dp
        
        IF (ABS(d) < tol) THEN
!           PRINT*, "iterations", it
!           PRINT*, d
          EXIT iter
        ENDIF
        
        
      ENDDO iter
      
      CALL eval_cubic_spline(t,ti,ax,bx,cx,dx,x)
      CALL eval_cubic_spline(t,ti,ay,by,cy,dy,y)
      
      IF (it >= maxit) THEN
        PRINT*, "MAX ITERATIONS EXCEEDED"
      ENDIF     

      
      RETURN 
      END SUBROUTINE newton      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      

      
      END MODULE calc_spline