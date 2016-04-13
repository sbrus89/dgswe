      MODULE calc_spline
      
      USE globals, ONLY: rp      
      
      IMPLICIT NONE
      
      
      
      CONTAINS
                  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE spline_init(num,nmax)
      
      USE globals, ONLY: rp,base,ax,bx,cx,dx,ay,by,cy,dy,dt, &
                         nfbnds,fbnds,fbnds_xy,nfbnd2el,fbnd2el, &
                         nverts,ctp,rpts     
      USE basis, ONLY: lglpts         
      
      IMPLICIT NONE
      
      INTEGER :: num
      INTEGER :: nmax
      INTEGER :: nd,bou,i,j,skip
      INTEGER :: el,et,nv,ged,led,edt
      INTEGER :: btype

      
      ! Count no normal flow boundaries and nodes
      
      num = 0       ! number of no normal flow boundaries
      nmax = 0      ! max number of nodes in any no normal flow boundary
      nfbnds = 0    ! number of total no normal flow boundary nodes
      DO bou = 1,base%nbou
        btype = base%fbseg(2,bou)
        IF( btype == 0 .OR. btype == 10 .OR. btype == 20  .OR. &   ! land boundaries
            btype == 1 .OR. btype == 11 .OR. btype == 21 ) THEN    ! island boundaries
        
          num = num + 1
          
          IF (base%fbseg(1,bou) > nmax) THEN
            nmax = base%fbseg(1,bou)
          ENDIF
          
          DO j = 1,base%fbseg(1,bou)
            nfbnds = nfbnds + 1            
          ENDDO
          
          
        ENDIF
      ENDDO
      
      ! Create a list of all no normal flow boundary nodes and their coordinates
      ! used to create the k-d tree to find spline coefficients to evaluate eval grid points
      
      ALLOCATE(fbnds(nfbnds),fbnds_xy(2,nfbnds))
      
      nfbnds = 0
      DO bou = 1,base%nbou
        btype = base%fbseg(2,bou)
        IF( btype == 0 .OR. btype == 10 .OR. btype == 20  .OR. &   ! land boundaries
            btype == 1 .OR. btype == 11 .OR. btype == 21 ) THEN    ! island boundaries
                            
          DO j = 1,base%fbseg(1,bou)
            nd = base%fbnds(j,bou)
            
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
      
      
      ALLOCATE(nfbnd2el(nfbnds),fbnd2el(base%mnepn,nfbnds))
      
      nfbnd2el = 0
      
      DO i = 1,nfbnds
        nd = fbnds(i)
        
        DO j = 1,base%nepn(nd)
          el = base%epn(j,nd)
          et = base%el_type(el)
          nv = nverts(et)
          
          DO led = 1,nv
            ged = base%el2ged(el,led)
            edt = base%ed_type(ged)
            
            IF (edt == 10) THEN
              nfbnd2el(i) = nfbnd2el(i) + 1
              fbnd2el(nfbnd2el(i),i) = el
            ENDIF
          ENDDO
          
        ENDDO
      ENDDO      
      
      

      PRINT "(A)", " "
      PRINT "(A,I5)", "Total number of type 0 normal flow boundaries ",num
      PRINT "(A,I5)", "Max number of nodes in a flow boundary ",nmax
      PRINT "(A)", " "
      
      ALLOCATE(ax(nmax,num),cx(nmax,num),bx(nmax-1,num),dx(nmax-1,num))
      ALLOCATE(ay(nmax,num),cy(nmax,num),by(nmax-1,num),dy(nmax-1,num)) 
      ALLOCATE(dt(nmax,num))
      
      ALLOCATE(rpts(ctp+1))
      
!       DO j = 0,ctp
!         rpts(j) = -1d0 + real(j,rp)*2d0/real(ctp,rp)   
!       ENDDO
      
      CALL lglpts(ctp,rpts)
      
      
      RETURN
      END SUBROUTINE spline_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      
      SUBROUTINE calc_cubic_spline(coord,bou,n,sig,a,b,c,d,dt)

      USE globals, ONLY: base
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: coord
      INTEGER, INTENT(IN) :: bou
      INTEGER, INTENT(IN) :: n
      INTEGER :: i      
      INTEGER :: info
      
      REAL(rp), INTENT(IN) :: sig
      REAL(rp), INTENT(OUT) , DIMENSION(:) :: a
      REAL(rp), INTENT(OUT), DIMENSION(:) :: b,c,d
      REAL(rp), INTENT(OUT), DIMENSION(:) :: dt
      REAL(rp) :: mult
      REAL(rp) :: x1,y1,x2,y2
      REAL(rp), DIMENSION(n) :: Ml,Md,Mu,v                 
     
     
      DO i = 1,n-1
      
        dt(i) = 1d0/(real(n,rp)-1d0)
        
!         x1 = base%xy(1,base%fbnds(i,bou))
!         y1 = base%xy(2,base%fbnds(i,bou))
!         
!         x2 = base%xy(1,base%fbnds(i+1,bou))
!         y2 = base%xy(2,base%fbnds(i+1,bou))
!         
!         dt(i) = sqrt((x2-x1)**2 + (y2-y1)**2)
      
      ENDDO
      

      ! Load nodal boundary coordinates 
      DO i = 1,n
        a(i) = base%xy(coord,base%fbnds(i,bou))
      ENDDO      
     



      IF (sig < 1d-12) THEN   ! No tension      
      
          ! Set up matrix 
          Ml(1) = 0d0
          Md(1) = 1d0
          Mu(1) = 0d0
          DO i = 2,n-1     
            Ml(i) = dt(i-1)
            Md(i) = 2d0*(dt(i-1)+dt(i))
            Mu(i) = dt(i)
          ENDDO
          Ml(n) = 0d0
          Md(n) = 1d0
          Mu(n) = 0d0

          ! Set up RHS
          c(1) = 0d0
          DO i = 2,n-1
            c(i) = 3d0/dt(i)*(a(i+1)-a(i)) - 3d0/dt(i-1)*(a(i)-a(i-1)) 
          ENDDO
          c(n) = 0d0
          
      ELSE   ! With tension (had to multiply the LHS of Palucci's notes by 2 )

          ! Set up matrix 
          Ml(1) = 0d0
          Md(1) = 1d0
          Mu(1) = 0d0
          DO i = 2,n-1        
            Ml(i) = 2d0*(1d0/dt(i-1) - sig/sinh(sig*dt(i-1)))/sig**2
            Md(i) = 2d0*(sig*cosh(sig*dt(i-1))/sinh(sig*dt(i-1)) - 1d0/dt(i-1) + &
                     sig*cosh(sig*dt(i))/sinh(sig*dt(i)) - 1d0/dt(i))/sig**2
            Mu(i) = 2d0*(1d0/dt(i) - sig/sinh(sig*dt(i)))/sig**2            
          ENDDO
          Ml(n) = 0d0
          Md(n) = 1d0
          Mu(n) = 0d0

          ! Set up RHS
          c(1) = 0d0
          DO i = 2,n-1
            c(i) = (a(i+1)-a(i))/dt(i) - (a(i)-a(i-1))/dt(i-1)
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
!         d(i) = (c(i+1)-c(i))/(3d0*dt)
!         b(i) = (a(i+1)-a(i))/dt - dt*(2d0*c(i)+c(i+1))/3d0      
      
        d(i) = (c(i+1)-c(i))/(3d0*dt(i))
        b(i) = (a(i+1)-a(i))/dt(i) - dt(i)*(2d0*c(i)+c(i+1))/3d0
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
      
      SUBROUTINE newton(r,dt,ti,xr,ax,bx,cx,dx,ay,by,cy,dy,x,y)
      
      IMPLICIT NONE      
      
      INTEGER :: it,maxit
      REAL(rp), INTENT(INOUT) :: r
      REAL(rp), INTENT(IN) :: dt
      REAL(rp), INTENT(IN) :: ti
      REAL(rp), INTENT(IN) :: xr(2)
      REAL(rp), INTENT(IN) :: ax,bx,cx,dx,ay,by,cy,dy
      REAL(rp), INTENT(OUT) :: x,y
      REAL(rp) :: t
      REAL(rp) :: f,fp,fpp,g,gp,gpp
      REAL(rp) :: d,dp,tol
      
      tol = 1d-8
      maxit = 1000
      
      t = .5d0*dt*(r + 1d0) + ti          ! initial guess for iteration 
      
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
      
      r =  2d0/dt*(t-ti)-1d0
      
      IF (it >= maxit) THEN
        PRINT "(A,E28.16)", "MAX ITERATIONS EXCEEDED IN FINDING EVALUATION PARAMETER, ERROR: ", ABS(d)        
      ENDIF     

      
      RETURN 
      END SUBROUTINE newton      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      

      
      END MODULE calc_spline