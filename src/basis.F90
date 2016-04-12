      MODULE basis
      
      USE lapack_interfaces
      USE globals, ONLY: rp      

      CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE element_basis(et,p,ndfs,npts,r,s,phi,dpdr,dpds)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: p,et,npts
      INTEGER, INTENT(OUT) :: ndfs
      REAL(rp), DIMENSION(:), INTENT(IN) :: r,s
      REAL(rp), DIMENSION(:,:), INTENT(OUT) :: phi
      REAL(rp), DIMENSION(:,:), OPTIONAL, INTENT(OUT) :: dpdr,dpds
      
      INTEGER :: calc_deriv
      
      
      
      IF (PRESENT(dpdr) .AND. PRESENT(dpds)) THEN
        calc_deriv = 1
      ELSE
        calc_deriv = 0
      ENDIF      
      
      
      
      IF (mod(et,2) == 1) THEN
      
        IF (calc_deriv == 1) THEN
          CALL tri_basis(p,ndfs,npts,r,s,phi,dpdr,dpds)
        ELSE
          CALL tri_basis(p,ndfs,npts,r,s,phi)
        ENDIF
        
      ELSE IF (mod(et,2) == 0) THEN
      
        IF (calc_deriv == 1) THEN
          CALL quad_basis(p,ndfs,npts,r,s,phi,dpdr,dpds)
        ELSE
          CALL quad_basis(p,ndfs,npts,r,s,phi)
        ENDIF
        
      ENDIF
      
      
      
      END SUBROUTINE element_basis

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        SUBROUTINE tri_basis(p,ndfs,npts,r,s,phi,dpdr,dpds)

        IMPLICIT NONE
        INTEGER, INTENT(IN) :: p,npts
        INTEGER, INTENT(OUT) :: ndfs
        REAL(rp), INTENT(IN) :: r(npts),s(npts) 
        REAL(rp), DIMENSION(:,:), INTENT(OUT) :: phi
        REAL(rp), DIMENSION(:,:), OPTIONAL, INTENT(OUT) :: dpdr,dpds         
        
        INTEGER :: m,i,j,pt
        INTEGER :: calc_deriv        
        REAL(rp) :: dpda,dpdb,dadr,dads,ii,tmp       
        REAL(rp) :: a(npts),b(npts)
        REAL(rp) :: Pi(npts),Pj(npts)
        REAL(rp) :: dPi(npts),dPj(npts)

        
        ndfs = (p+1)*(p+2)/2
        
        IF (PRESENT(dpdr) .AND. PRESENT(dpds)) THEN
          calc_deriv = 1
        ELSE
          calc_deriv = 0
        ENDIF        
      
        ! Change quadrature points from r,s (master element) to a,b extended coordinates
        DO pt = 1,npts
          IF(abs(s(pt) - 1d0) > 1d-14) THEN
            a(pt) = 2d0*(1d0+r(pt))/(1d0-s(pt))-1d0 
          ELSE
            a(pt) = -1d0
          ENDIF
          b(pt) = s(pt)
        ENDDO
        
        ! Calculate basis function values and derivative values at area quadrature points
        m = 0
        DO i = 0,p
          ii = real(i,rp)
          DO j = 0,p-i

            m = m+1

            CALL jacobi(0    ,0,i,a,npts,Pi)
            CALL jacobi(2*i+1,0,j,b,npts,Pj)           

            ! Calculate function values
            DO pt = 1,npts 
              ! orthonormal
              ! phi(m,pt) = sqrt(2d0)*Pi(pt)*Pj(pt)*(1d0-b(pt))**i
              
              ! phi_1 = 1
              phi(m,pt) = 2d0*Pi(pt)*Pj(pt)*(1d0-b(pt))**i
            ENDDO

            IF (calc_deriv == 1) THEN
              CALL djacobi(0    ,0,i,a,npts,dPi)
              CALL djacobi(2*i+1,0,j,b,npts,dPj)                          

              ! Calculate derivative values
              DO pt = 1,npts
              
!                 dadr = 2d0/(1d0-s(pt))
!                 dads = 2d0*(1d0+r(pt))/(1d0-s(pt))**2d0
!                
!                 ! orthonormal
!                 ! dpda = sqrt(2d0)*dPi(pt)*Pj(pt)*(1d0-b(pt))**ii
!                 ! dpdb = sqrt(2d0)*Pi(pt)*(dPj(pt)*(1d0-b(pt))**ii - ii*(1d0-b(pt))**(ii-1d0)*Pj(pt))
!                 
!                 ! phi_1 = 1
!                 dpda = 2d0*dPi(pt)*Pj(pt)*(1d0-b(pt))**ii
!                 dpdb = 2d0*Pi(pt)*(dPj(pt)*(1d0-b(pt))**ii - ii*(1d0-b(pt))**(ii-1d0)*Pj(pt))
!               
!                 dpdr(m,pt) = dpda*dadr
!                 dpds(m,pt) = dpda*dads + dpdb

                ! correction for corner nodes 
                dpdr(m,pt) = dPi(pt)*Pj(pt)
                IF (i>0) THEN
                  dpdr(m,pt) = dpdr(m,pt)*2d0*(1d0-b(pt))**(i-1)
                ENDIF
                
                dpds(m,pt) = dPi(pt)*Pj(pt)*(1d0+a(pt))
                IF (i>0) THEN
                  dpds(m,pt) = dpds(m,pt)*(1d0-b(pt))**(i-1)
                ENDIF
                
                tmp = dPj(pt)*(1d0-b(pt))**i
                IF (i>0) THEN
                  tmp = tmp - ii*Pj(pt)*(1d0-b(pt))**(i-1)
                ENDIF
                dpds(m,pt) = dpds(m,pt) + Pi(pt)*tmp              

                dpdr(m,pt) = 2d0*dpdr(m,pt)
                dpds(m,pt) = 2d0*dpds(m,pt)                             
                
              ENDDO
            ENDIF

          ENDDO
        ENDDO  
        
        END SUBROUTINE tri_basis        
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
        
        SUBROUTINE quad_basis(p,ndfs,npts,r,s,phi,dpdr,dpds)
        
        IMPLICIT NONE
        
        INTEGER, INTENT(IN) :: p,npts
        INTEGER, INTENT(OUT) :: ndfs
        REAL(rp), INTENT(IN) :: r(npts),s(npts)
        REAL(rp), DIMENSION(:,:), INTENT(OUT) :: phi
        REAL(rp), DIMENSION(:,:), OPTIONAL, INTENT(OUT) :: dpdr,dpds
        
        INTEGER :: m,i,j,pt    
        REAL(rp) :: Pi(npts),Pj(npts)
        REAL(rp) :: dPi(npts),dPj(npts)  
        INTEGER :: calc_deriv
        
        ndfs = (p+1)**2
        
        IF (PRESENT(dpdr) .AND. PRESENT(dpds)) THEN
          calc_deriv = 1
        ELSE
          calc_deriv = 0
        ENDIF
        
        ! Calculate basis function values and derivative values at area quadrature points
        m = 0
        DO i = 0,p
          DO j = 0,p

            m = m+1

            CALL jacobi(0,0,i,r,npts,Pi)
            CALL jacobi(0,0,j,s,npts,Pj)              

            ! Calculate function values
            DO pt = 1,npts 
              phi(m,pt) = 2d0*Pi(pt)*Pj(pt)
            ENDDO

            IF (calc_deriv == 1) THEN
              CALL djacobi(0,0,i,r,npts,dPi)
              CALL djacobi(0,0,j,s,npts,dPj)
            
              ! Calculate derivative values
              DO pt = 1,npts
                dpdr(m,pt) = 2d0*dPi(pt)*Pj(pt)
                dpds(m,pt) = 2d0*Pi(pt)*dPj(pt)
              ENDDO
            ENDIF

          ENDDO
        ENDDO

        
        END SUBROUTINE quad_basis            
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE jacobi(alpha_i,beta_i,deg,x,npts,v)

        IMPLICIT NONE
        INTEGER, INTENT(IN) :: alpha_i,beta_i,deg,npts
        REAL(rp), INTENT(IN) :: x(npts)
        REAL(rp), INTENT(OUT) :: v(npts)
        INTEGER :: i,j,np1
        REAL(rp) :: pnm1(npts),pn(npts),pnp1(npts)
        REAL(rp) :: alpha,beta,an,anp1,bn,n

        alpha = real(alpha_i,rp)
        beta = real(beta_i,rp)

        ! Calculate constant P0
        DO i = 1,npts
          pnm1(i) = sqrt(2d0**(-alpha-beta-1d0)*real(fact(alpha_i+beta_i+1),rp) &
                         /(real(fact(alpha_i),rp)*real(fact(beta_i),rp))) 
          v(i) = pnm1(i)          
        ENDDO

        IF (deg == 0) THEN
          RETURN
        ENDIF

        ! Calculate linear P1
        DO i = 1,npts
          pn(i) = .5d0*pnm1(i)*sqrt((alpha+beta+3d0)/((alpha+1d0)*(beta+1d0)))*((alpha+beta+2d0)*x(i) + alpha-beta)
          v(i) = pn(i)
        ENDDO

        IF (deg == 1) THEN
          RETURN
        ENDIF

        ! Calculate a1 and b1
        n = 1d0
        an = (2d0/(2d0*n+alpha+beta))*sqrt((n*(n+alpha+beta)*(n+alpha)*(n+beta))/((2d0*n+alpha+beta-1d0)*(2d0*n+alpha+beta+1d0)))
        bn = -(alpha*alpha-beta*beta)/((2d0*n+alpha+beta)*(2d0*n+alpha+beta+2d0))

        ! Loop for Pn+1
        DO np1 = 2,deg 
          n = real(np1,rp)
          anp1 = (2d0/(2d0*n+alpha+beta)) &
                  *sqrt((n*(n+alpha+beta)*(n+alpha)*(n+beta))/((2d0*n+alpha+beta-1d0)*(2d0*n+alpha+beta+1d0)))
          DO i = 1,npts
            pnp1(i) = (1d0/anp1)*(x(i)*pn(i)-an*pnm1(i)-bn*pn(i))
            pnm1(i) = pn(i)
            pn(i) = pnp1(i)
            v(i) = pn(i)
          ENDDO
          an = anp1
          bn = -(alpha*alpha-beta*beta)/((2d0*n+alpha+beta)*(2d0*n+alpha+beta+2d0))         
        ENDDO

        RETURN
      END SUBROUTINE jacobi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE djacobi(alpha_i,beta_i,deg_i,x,npts,dP)

        IMPLICIT NONE
        INTEGER, INTENT(IN) :: alpha_i,beta_i,deg_i,npts
        REAL(rp), INTENT(IN) :: x(npts)
        REAL(rp), INTENT(OUT) :: dP(npts)
        REAL(rp):: v(npts),deg,alpha,beta
        INTEGER :: i,pt

        deg = real(deg_i,rp)
        alpha = real(alpha_i,rp)
        beta = real(beta_i,rp)
        
        IF (deg == 0) THEN

          DO pt = 1,npts
            dP(pt) = 0d0
          ENDDO

          RETURN
        ENDIF

        CALL jacobi(alpha_i+1,beta_i+1,deg_i-1,x,npts,v)

        DO pt = 1,npts

          dP(pt) = sqrt(deg*(deg+alpha+beta+1d0))*v(pt)

        ENDDO

        RETURN
      END SUBROUTINE djacobi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE adcirc_basis(p,ndof,npts,rv,sv,phi)
      
      IMPLICIT NONE
      
      INTEGER :: pt,dof,i
      INTEGER :: p,ndof,npts
      REAL(rp) :: rv(npts),sv(npts),s,r
      REAL(rp) :: phi(ndof*npts)
      
      DO pt = 1,npts
        r = rv(pt)
        s = sv(pt)
        phi(pt) = 1d0
        phi(npts + pt) = 0.5d0*(1d0+3d0*s)
        phi(2*npts + pt) = .5d0*(1d0+2d0*r+s)       
      ENDDO        
        
      IF (p > 1) THEN
        DO pt = 1,npts
          r = rv(pt)
          s = sv(pt)        
          phi(3*npts + pt) = -.5d0 + s + .5d0*(5d0*s**2)
          phi(4*npts + pt) = .25d0*(1d0 + 2d0*r + s)*(3d0 + 5d0*s)
          phi(5*npts + pt) = .25d0*(1d0 + 6d0*r**2 + 4d0*s + s**2 + 6d0*r*(1d0 + s))
        ENDDO
      ENDIF
        
      IF (p > 2) THEN
        DO pt = 1,npts 
          r = rv(pt)
          s = sv(pt)          
          phi(6*npts + pt) = .125*(-3d0 - 15d0*s + 15d0*s**2 + 35d0*s**3)
          phi(7*npts + pt) = .125*(1d0 + 2d0*r + s)*(1d0 + 18d0*s + 21d0*s**2)
          phi(8*npts + pt) = .125*(5d0 + 7d0*s)*(1d0 + 6d0*r**2 + 4d0*s + s**2 + 6d0*r*(1d0 + s))
          phi(9*npts + pt) = .125*(1d0 + 20d0*r**3 + 9d0*s + 9d0*s**2 + s**3 + 30d0*r**2*(1d0 + s) + 12d0*r*(1d0 + 3d0*s + s**2))
        ENDDO
      ENDIF

      IF (p > 3) THEN
        DO pt = 1,npts
          r = rv(pt)
          s = sv(pt)          
          phi(10*npts + pt) = (1d0/8d0)*(3d0 - 12d0*s - 42d0*s**2 + 28d0*s**3 + 63d0*s**4)
          phi(11*npts + pt) = (1d0/4d0)*(1d0 + 2d0*r + s)*(-2d0 + 21d0*s**2 + 21d0*s**3)
          phi(12*npts + pt) = (1d0/4d0)*(2d0 + 10d0*s + 9d0*s**2)*(1d0 + 6d0*r**2 + 4d0*s + s**2 + 6d0*r*(1d0 + s))
          phi(13*npts + pt) = (1d0/16d0)*(7d0 + 9d0*s)*(1d0 + 20d0*r**3 + 9d0*s + 9d0*s**2 + s**3 + 30d0*r**2*(1 + s) + 12d0*r*(1d0 + 3d0*s + s**2))
          phi(14*npts + pt) = (1d0/16d0)*(1d0 + 70d0*r**4 + 16d0*s + 36d0*s**2 + 16d0*s**3 + s**4 + 140d0*r**3*(1d0 + s) + 30d0*r**2*(3d0 + 8d0*s + 3d0*s**2) + 20d0*r*(1d0 + 6d0*s + 6d0*s**2 + s**3))
        ENDDO
      ENDIF
      
      RETURN
      END SUBROUTINE      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE simple_basis(p,ndfs,npts,r,s,phi)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: p
      INTEGER, INTENT(OUT) :: ndfs
      INTEGER, INTENT(IN) :: npts
      REAL(rp), DIMENSION(:), INTENT(IN) :: r
      REAL(rp), DIMENSION(:), INTENT(IN) :: s
      REAL(rp), DIMENSION(:,:), INTENT(OUT) :: phi
      
      INTEGER :: i,j,m,pt
      
      ndfs = (p+1)*(p+2)/2      
      
      DO pt = 1,npts
      
        m = 0
        DO i = 0,p
          DO j = 0,p-i
            m = m + 1            
            phi(m,pt) = r(pt)**i*s(pt)**j
          ENDDO
        ENDDO      
      
      ENDDO
      
      RETURN
      END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE linear(n,r,s,phil)     

        IMPLICIT NONE
        INTEGER :: n
        REAL(rp), DIMENSION(:,:), INTENT(OUT) :: phil
        REAL(rp), DIMENSION(:), INTENT(IN) :: r,s
        
        INTEGER :: pt

        
        DO pt = 1,n

          phil(1,pt) = -.5d0*(r(pt)+s(pt))
          phil(2,pt) = .5d0*(r(pt)+1d0)
          phil(3,pt) = .5d0*(s(pt)+1d0)  
          
        ENDDO

        RETURN
      END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      FUNCTION fact(n) RESULT (n_fact)
        IMPLICIT NONE
        INTEGER :: n_fact
        INTEGER, INTENT(IN) :: n
        INTEGER :: k

        n_fact = 1
        DO k = 2,n
          n_fact = k*n_fact
        ENDDO

      END FUNCTION fact
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    

      SUBROUTINE element_nodes(et,space,p,n,r,s,ext)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: et,space,p
      INTEGER, INTENT(OUT) :: n
      REAL(rp), DIMENSION(:), INTENT(OUT) :: r,s
      INTEGER, INTENT(IN), OPTIONAL :: ext
      
      INTEGER :: nv
      
      IF (mod(et,2) == 1) THEN
        nv = 3
        CALL tri_nodes(1,p,n,r,s)
      ELSE IF  (mod(et,2) == 0) THEN
        nv = 4
        CALL quad_nodes(1,p,n,r,s)
      ENDIF
      
      IF (PRESENT(ext)) THEN
        CALL nodes_extract(ext,nv,p,n,r,s)
      ENDIF
      
      END SUBROUTINE element_nodes
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
      SUBROUTINE tri_nodes(space,p,n,r,s)
!c     Calculates nodal set for well-behaved interpolation
!c     p (input) : the order of the basis
!c     n (input): the number of points in the triangle
!c     r(np) (output) ; nodal r-coordinates for master element
!c     s(np) (output) : nodal s-coordinates for master element 

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: p,space 
      INTEGER, INTENT(OUT) :: n
      REAL(rp), DIMENSION(:), INTENT(OUT) :: r,s      
      
      INTEGER :: i,j,m
      REAL(rp) :: ii,jj,a,dx,tol
      REAL(rp) :: aopt(15)
      REAL(rp) :: lgl(p+1),xeq(p+1)           
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: x,y
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: l1,l2,l3      
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: b1,b2,b3
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: wx,wy,w1,w2,w3
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: w1e,w2e,w3e
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: var1,var2,var3
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: w1mat,w2mat,w3mat     
     
     

!c     Define optimal alpha values
      DATA aopt/0d0,0d0,1.4152d0,0.1001d0,0.2751d0,0.9800d0,1.0999d0 &
               ,1.2832d0,1.3648d0,1.4773d0,1.4959d0,1.5743d0,1.5770d0 &
               ,1.6223d0,1.6258d0/

      IF(p.lt.16) THEN
        a = aopt(p)
      ELSE
        a = 5d0/3d0
      ENDIF
      
      n = (p+1)*(p+2)/2
      
      ALLOCATE(x(n),y(n))            
      ALLOCATE(l1(n),l2(n),l3(n))
      ALLOCATE(b1(n),b2(n),b3(n))
      ALLOCATE(wx(n),wy(n),w1(n),w2(n),w3(n))
      ALLOCATE(w1e(n),w2e(n),w3e(n))
      ALLOCATE(var1(n),var2(n),var3(n))
      ALLOCATE(w1mat(n,p+1),w2mat(n,p+1),w3mat(n,p+1))
      
      IF (space == -1) THEN
        CALL tri_nodes_random(n,r,s)
        RETURN
      ENDIF

!c     Create equally distributed nodes on equalateral triangle

      m = 0

! Edge 3 points
      i = 1
      ii = dble(i)
      DO j = 1,p+1
        jj = dble(j)
        m = m+1
        l1(m) = (ii-1d0)/dble(p)
        l3(m) = (jj-1d0)/dble(p)
      ENDDO

! Edge 1 points
      DO i = 2,p+1
        ii = dble(i)
        j = p+2-i 
        jj = dble(j)
        m = m+1
        l1(m) = (ii-1d0)/dble(p)
        l3(m) = (jj-1d0)/dble(p)
      ENDDO


! Edge 2 points
      DO i = p,2,-1
        ii = dble(i)
        j = 1
        jj = dble(j)
        m = m+1
        l1(m) = (ii-1d0)/dble(p)
        l3(m) = (jj-1d0)/dble(p)
      ENDDO

! Internal points
      DO i = 2,p
        ii = dble(i)
        DO j = 2,p+1-i 
          jj = dble(j)
          m = m+1
          l1(m) = (ii-1d0)/dble(p)
          l3(m) = (jj-1d0)/dble(p)
        ENDDO
      ENDDO

!       do i = 1,p+1
!         ii = dble(i)
!         do j = 1,p+2-i 
!           jj = dble(j)
!           m = m+1
!           l1(m) = (ii-1d0)/dble(p)
!           l3(m) = (jj-1d0)/dble(p)
!         enddo
!       enddo

      DO i = 1,m
        l2(i) = 1d0-l1(i)-l3(i)
        x(i) = -l2(i)+l3(i)
        y(i) = (-l2(i)-l3(i)+2d0*l1(i))/dsqrt(3d0)
      ENDDO

      IF(space.eq.0)then
        CALL xytors(n,x,y,r,s) ;
        RETURN
      ENDIF

!c     Compute blending function
      DO i = 1,n
        b1(i) = 4d0*l2(i)*l3(i)
        b2(i) = 4d0*l3(i)*l1(i)
        b3(i) = 4d0*l2(i)*l1(i)
      ENDDO

!c     Compute 1-D Legendre-Gauss-Lobotto points on (-1,1)
      CALL lglpts(p,lgl)

!       do i = 1,p+1
!         write(40,*) lgl(i)
!       enddo

!c     Compute equadistant points on (-1,1)
      dx = 2d0/dble(p)
      xeq(1)=-1d0
      DO i = 2,p+1
        xeq(i) = xeq(i-1)+dx
      ENDDO
      xeq(p+1) = 1d0

!c     Compute warping function arugments
      DO i = 1,n
        w1e(i) = l3(i)-l2(i)
        w2e(i) = l1(i)-l3(i)
        w3e(i) = l2(i)-l1(i)
      ENDDO

!c     Elvaluate Lagrange polynomials based on equa-spaced nodes
      CALL lagrange(p+1,n,xeq,w1e,w1mat)
      CALL lagrange(p+1,n,xeq,w2e,w2mat)
      CALL lagrange(p+1,n,xeq,w3e,w3mat)

!       do j = 1,n
!         write(50,41) (w1mat(j,i), i = 1,p+1)
!       enddo
!  41   format(16000(e24.17,1x))

!c     Compute 1 dimensional mapping functions for edges

      tol = 1d0 - 1.0d-10

      DO i = 1,n
        IF(dabs(w1e(i)).lt.tol) THEN
          var1(i) = 1d0
        ELSE
          var1(i) = 0d0
        ENDIF
        IF(dabs(w2e(i)).lt.tol) THEN
          var2(i) = 1d0
        ELSE
          var2(i) = 0d0
        ENDIF
        IF(dabs(w3e(i)).lt.tol) THEN
          var3(i) = 1d0
        ELSE
          var3(i) = 0d0
        ENDIF
      ENDDO

      DO j = 1,n
        w1(j) = 0d0
        w2(j) = 0d0
        w3(j) = 0d0
      ENDDO

      DO j = 1,n
        DO i = 1,p+1
         w1(j) = w1(j) + (lgl(i)-xeq(i))*w1mat(j,i)
         w2(j) = w2(j) + (lgl(i)-xeq(i))*w2mat(j,i)
         w3(j) = w3(j) + (lgl(i)-xeq(i))*w3mat(j,i)
        ENDDO
        w1(j) = w1(j)/(1d0-(var1(j)*w1e(j))**2)+w1(j)*(var1(j)-1d0)
        w2(j) = w2(j)/(1d0-(var2(j)*w2e(j))**2)+w2(j)*(var2(j)-1d0)
        w3(j) = w3(j)/(1d0-(var3(j)*w3e(j))**2)+w3(j)*(var3(j)-1d0)
      ENDDO

!c     Apply gerneralized warping function to equally distributed nodes on equalateral triangle
      DO i = 1,n
        x(i) = x(i)+(1d0+(a*l1(i))**2d0)*b1(i)*w1(i) &
                   -.5d0*(1d0+(a*l2(i))**2d0)*b2(i)*w2(i) &
                   -.5d0*(1d0+(a*l3(i))**2d0)*b3(i)*w3(i)
        y(i) = y(i)+dsqrt(3d0)/2d0*(1d0+(a*l2(i))**2d0)*b2(i)*w2(i) &
                   -dsqrt(3d0)/2d0*(1d0+(a*l3(i))**2d0)*b3(i)*w3(i)
      ENDDO

      CALL xytors(n,x,y,r,s)
      
!         print*,' '
!         print*, 'Triangle Points'
!         do i = 1,n
!           print 21, r(i),s(i)
!         enddo
!         print*, ' '

 21   format(28(f10.4,1x))
 61   format(16000(e24.17,1x))      

      RETURN
      END SUBROUTINE
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      

      SUBROUTINE quad_nodes(space,p,n,r,s)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: space,p
      INTEGER, INTENT(OUT) :: n      
      REAL(rp), DIMENSION(:), INTENT(OUT) :: r,s
      

      INTEGER :: i,k,j,nl,m
      REAL(rp) :: xi(p+1)

      
      
      n = (p+1)**2
      
      ! Get 1-D LGL points 
      
      IF (space == 1) THEN
        CALL lglpts(p,xi)
      ELSE IF (space == -1) THEN
        CALL quad_nodes_random(n,r,s)
        RETURN
      ELSE
        xi(1) = -1d0
        DO i = 1,p-1
          xi(i+1) = xi(i) + 2d0/real(p,rp)
        ENDDO
        xi(p+1) = 1d0
      ENDIF
        
      ! Do tensor product, ordering the nodes counter-clockwise

      ! Find number of loops around refence quad element, excluding middle points
      IF (p <= 2)THEN
        nl = 1
      ELSE
        IF(mod(p,2) == 1) THEN
          nl = p-1
        ELSE IF (mod(p,2) == 0) THEN
          nl = p-2
        ENDIF
      ENDIF
        
      m = 1
      DO k = 1,nl ! loop over number of loops
         
        ! Edge 4
        DO i = k,p+1 - (k-1)
          j = k
          r(m) = xi(i)
          s(m) = xi(j)
          
          m = m+1
        ENDDO
        
        ! Edge 1
        DO j = 2 + (k-1),p+1 - (k-1)
          i = p+1 - (k-1)
          r(m) = xi(i)
          s(m) = xi(j)
          
          m = m+1
        ENDDO

        ! Edge 2
        DO i = p - (k-1),1 + (k-1),-1
          j = p+1 - (k-1)
          r(m) = xi(i)
          s(m) = xi(j)
          
          m = m+1      
        ENDDO

        ! Edge 3
        DO j = p - (k-1),2 + (k-1),-1
          i = 1 + (k-1)
          r(m) = xi(i)
          s(m) = xi(j)

          m = m+1
        ENDDO    
          
      ENDDO
        
      ! middle point
      IF (mod(p+1,2) == 1) THEN
        i = p/2 + 1
        r(m) = xi(i)
        s(m) = xi(i)
          
      ENDIF   
                
!         PRINT*, ' '
!         PRINT*, 'Quadrilateral Points'
!         DO m = 1,n
!           PRINT("(2(f10.4))"), r(m),s(m)
!         ENDDO         
      
      RETURN
      END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      subroutine lglpts(n,r)
!c     Returns the n+1 nth order Legendre-Gauss-Lebotto points on an iterval of (-1,1) in r.
!c     n (input) : order 
!c     r(n+1) (output) : Legendre-Gauss-Lebotto points

      implicit none
      integer n,nn,i
      real(rp) ii,h
      integer info
      real(rp) r(n+1),d(n+1),e(n)
      real(rp) z(n+1,n+1),work(1,2*n-2)

      if(n.eq.1) then
        r(1) = -1d0
        r(2) = 1d0
        go to 10
      endif
     
      nn = n - 2

      if(nn.eq.0) then
        r(1) = -1d0
        r(2) = 0d0
        r(3) = 1d0
        go to 10
      endif

      do i = 1,nn+1
        d(i) = 0d0
      enddo
      do i = 1,nn
        ii = dble(i)
        h = 2d0*(ii-1.0)+2d0
        e(i) = 2d0/(h+2d0)*dsqrt(ii*(ii+2d0)*(ii+1d0)*(ii+1d0)/(h+1d0)/(h+3d0))
      enddo
      
      call dsteqr('i',nn+1,d,e,z,nn+1,work,info)

      r(1) = -1d0
      do i = 1,nn+1   
        r(i+1) = d(i)
      enddo
      r(n+1) = 1d0
   
   
 10   continue 
 
!       print*,' '
!       print*, 'LGL Points'
!       do i = 1,n+1
!         print 21, r(i)
!       enddo
!       print*, ' '

 21   format(28(f10.4,1x))      

      return 
      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      subroutine lagrange(nn,ne,xn,xe,pmat)
!c     Evaluates the Lagrange polynomial associated with nodes nn at points xe
!c     nn (input) : number of nodes
!c     ne (input) : number of evaluation points
!c     xn(nn) (input) : vector of nodal values
!c     xe(ne) (input) : vector of evaluation points
!c     pmat(nn,ne) (output) : matrix with structure:
!c     [L0(x1)  L1(x1)  L2(x1)]
!c     [L0(x2)  L1(x2)  L2(x2)]
!c     [L0(x3)  L1(x3)  L3(x3)]

      implicit none
      integer nn,ne,n,i,j
      real(rp) p,xn(nn),xe(ne),pmat(ne,nn)

      do n = 1,ne
        do i = 1,nn
           p = 1d0
           do j = 1,i-1
             p = p*(xe(n)-xn(j))/(xn(i)-xn(j))
           enddo
           do  j = i+1,nn
             p = p*(xe(n)-xn(j))/(xn(i)-xn(j))
           enddo
           pmat(n,i) = p
        enddo
      enddo

      return
      end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      subroutine xytors(np,x,y,r,s)  
      
      implicit none
      integer np,i
      real(rp) x(np),y(np)
      real(rp) r(np),s(np)
      real(rp) l1(np),l2(np),l3(np)

!c     Convert from x,y to r,s
      do i = 1,np
        l1(i) = (dsqrt(3d0)*y(i)+1d0)/3d0
        l2(i) = (-3d0*x(i)-dsqrt(3d0)*y(i)+2d0)/6d0
        l3(i) = (3d0*x(i)-dsqrt(3d0)*y(i)+2d0)/6d0
        r(i) = -l2(i)+l3(i)-l1(i)
        s(i) = -l2(i)-l3(i)+l1(i)
      enddo

      return
      end subroutine
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      SUBROUTINE tri_nodes_random(n,r,s)            
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: n
      REAL(rp), DIMENSION(:), INTENT(OUT) :: r
      REAL(rp), DIMENSION(:), INTENT(OUT) :: s
      
      INTEGER :: i
      REAL(rp) :: rnd
      
      CALL RANDOM_SEED(put=(/3/))  
      
      DO i = 1,n
        CALL RANDOM_NUMBER(rnd)        
        r(i) = 2d0*rnd - 1d0
        
        CALL RANDOM_NUMBER(rnd)         
        s(i) = (r(i)-1d0)*rnd - 1d0
      ENDDO
      
      RETURN
      END SUBROUTINE tri_nodes_random      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      SUBROUTINE quad_nodes_random(n,r,s)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: n
      REAL(rp), DIMENSION(:), INTENT(OUT) :: r
      REAL(rp), DIMENSION(:), INTENT(OUT) :: s 
      
      INTEGER :: i
      REAL(rp) :: rnd
      
      CALL RANDOM_SEED(put=(/3/))      
      
      DO i = 1,n
        CALL RANDOM_NUMBER(rnd)       
        r(i) = 2d0*rnd-1d0
        
        CALL RANDOM_NUMBER(rnd)         
        s(i) = 2d0*rnd-1d0
      ENDDO      
      
      RETURN
      END SUBROUTINE quad_nodes_random

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

     SUBROUTINE nodes_extract(ext,nv,p,n,r,s)
     
     IMPLICIT NONE
     
     INTEGER, INTENT(IN) :: ext
     INTEGER, INTENT(IN) :: p
     INTEGER, INTENT(INOUT) :: n
     REAL(rp), DIMENSION(:), INTENT(INOUT) :: r
     REAL(rp), DIMENSION(:), INTENT(INOUT) :: s
     
     INTEGER :: i,j
     INTEGER :: nv,led
     INTEGER :: n_ext
     REAL(rp) :: r_ext(n),s_ext(n)

     
  
         
     IF (ext <= 0) THEN          
     
       RETURN
       
     ELSE IF (ext > nv) THEN
       
       n_ext = 0     
       DO i = nv*(p-1)+nv+1,n
         n_ext = n_ext + 1
           
         r_ext(n_ext) = r(i)
         s_ext(n_ext) = s(i)           
       ENDDO       
              
     ELSE 
     
       led = ext     
       n_ext = p+1     
       DO i = 1,n_ext
         j = mod(led,nv)*p + i
         
         IF (j == nv*p + 1) THEN
           j = 1
         ENDIF 
           
         r_ext(i) = r(j)
         s_ext(i) = s(j)           
       ENDDO             
         
     ENDIF
     
     
     r = 0
     s = 0 
     DO i = 1,n_ext
       r(i) = r_ext(i)
       s(i) = s_ext(i)       
     ENDDO
     n = n_ext
     
     RETURN
     END SUBROUTINE nodes_extract
     

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      END MODULE basis
