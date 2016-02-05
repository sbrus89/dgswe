      MODULE calc_spline
      
      USE globals, ONLY: rp      
      
      IMPLICIT NONE
      
      
      
      CONTAINS
                  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE spline_init(num)
      
      USE globals, ONLY: rp,base,ax,bx,cx,dx,ay,by,cy,dy, &
                         nfbnds,fbnds,fbnds_xy, &
                         tree_xy,closest,srchdp
      USE kdtree2_module                     
      
      IMPLICIT NONE
      
      INTEGER :: num
      INTEGER :: nmax
      INTEGER :: nd,seg,i,j,skip
      INTEGER :: segtype

      
      ! Count no normal flow boundaries and nodes
      
      num = 0       ! number of no normal flow segments
      nmax = 0      ! max number of nodes in any no normal flow segment
      nfbnds = 0    ! number of total no normal flow boundary nodes
      DO seg = 1,base%nbou
        segtype = base%fbseg(2,seg)
        IF( segtype == 0 .OR. segtype == 10 .OR. segtype == 20  .OR. &   ! land boundaries
            segtype == 1 .OR. segtype == 11 .OR. segtype == 21 ) THEN    ! island boundaries
        
          num = num + 1
          
          IF (base%fbseg(1,seg) > nmax) THEN
            nmax = base%fbseg(1,seg)
          ENDIF
          
          DO j = 1,base%fbseg(1,seg)
            nfbnds = nfbnds + 1            
          ENDDO
          
          
        ENDIF
      ENDDO
      
      ! Create a list of all no normal flow boundary nodes and their coordinates
      ! used to create the k-d tree to find spline coefficients to evaluate eval grid points
      
      ALLOCATE(fbnds(nfbnds),fbnds_xy(2,nfbnds))
      
      nfbnds = 0
      DO seg = 1,base%nbou
        segtype = base%fbseg(2,seg)
        IF( segtype == 0 .OR. segtype == 10 .OR. segtype == 20  .OR. &   ! land boundaries
            segtype == 1 .OR. segtype == 11 .OR. segtype == 21 ) THEN    ! island boundaries
                            
          DO j = 1,base%fbseg(1,seg)
            nd = base%fbnds(j,seg)
            
            skip = 0           ! skip duplicate nodes
            DO i = 1,nfbnds
              IF (fbnds(i) == nd) THEN
                skip = 1
              ENDIF
            ENDDO
          
            IF (skip == 0) THEN
              nfbnds = nfbnds + 1            

              fbnds(nfbnds) = nd
              fbnds_xy(1:2,nfbnds) = base%xy(1:2,nd)
            ENDIF
          ENDDO
          
          
        ENDIF
      ENDDO
      
      
      
      tree_xy => kdtree2_create(fbnds_xy(1:2,1:nfbnds), rearrange=.true., sort=.true.)
      ALLOCATE(closest(srchdp))
      
      

      PRINT "(A)", " "
      PRINT "(A,I5)", "Total number of type 0 normal flow boundaries ",num
      PRINT "(A,I5)", "Max number of nodes in a flow boundary segment ",nmax
      PRINT "(A)", " "
      
      ALLOCATE(ax(nmax),cx(nmax),bx(nmax-1),dx(nmax-1))
      ALLOCATE(ay(nmax),cy(nmax),by(nmax-1),dy(nmax-1))      
      
      
      RETURN
      END SUBROUTINE spline_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      
      SUBROUTINE calc_cubic_spline(sig,n,a,b,c,d,dt)


      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: n
      INTEGER :: i      
      INTEGER :: info
      
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
          c(1) = 0d0
          DO i = 2,n-1
            c(i) = 3d0/dt*(a(i+1)-2d0*a(i)+a(i-1))
          ENDDO
          c(n) = 0d0
          
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
          c(1) = 0d0
          DO i = 2,n-1
            c(i) = 1d0/dt*(a(i+1)-2d0*a(i)+a(i-1))
          ENDDO
          c(n) = 0d0
          
      ENDIF
      
      ! Solve system for c coefficients
      CALL DGTSV(n,1,Ml,Md,Mu,c,n,info)      

      
!       ! Solve system for c coefficients, forward sweep
!       ! RHS is v 
!       DO i = 2,n
!         mult = Ml(i)/Md(i-1)
!         Md(i) = Md(i) - mult*Mu(i-1)
!         v(i) = v(i) - mult*v(i-1)
!       ENDDO
! 
!       ! Solve system for c coefficients, backward sweep
!       c(n) = v(n)/Md(n)
!       DO i = n-1,1,-1
!         c(i) = (v(i) - Mu(i)*c(i+1))/Md(i)
!       ENDDO
 
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
      
      SUBROUTINE newton(t,ti,xr,ax,bx,cx,dx,ay,by,cy,dy,x,y)
      
      IMPLICIT NONE      
      
      INTEGER :: it,maxit
      REAL(rp), INTENT(INOUT) :: t
      REAL(rp), INTENT(IN) :: ti
      REAL(rp), INTENT(IN) :: xr(2)
      REAL(rp), INTENT(IN) :: ax,bx,cx,dx,ay,by,cy,dy
      REAL(rp), INTENT(OUT) :: x,y
      REAL(rp) :: f,fp,fpp,g,gp,gpp
      REAL(rp) :: d,dp,tol
      
      tol = 1d-8
      maxit = 10000
      
iter: DO it = 1,maxit

        CALL eval_cubic_spline(t,ti,ax,bx,cx,dx,f,fp,fpp)               
        CALL eval_cubic_spline(t,ti,ay,by,cy,dy,g,gp,gpp)
        
        d = 2d0*(f-xr(1))*fp + 2d0*(g-xr(2))*gp
        dp = 2d0*fp**2 + 2d0*(f-xr(1))*fpp + 2d0*gp**2 + 2d0*(g-xr(2))*gpp
        
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
        PRINT*, "MAX ITERATIONS EXCEEDED, ERROR: ", ABS(d)        
      ENDIF     

      
      RETURN 
      END SUBROUTINE newton      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      

      
      END MODULE calc_spline