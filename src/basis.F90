      MODULE basis

      CONTAINS

      SUBROUTINE area_basis()

        USE globals, ONLY: rp,nel_type,nqpta,mnqpta,wpta,qpta,ndof,mndof, &
                           phia,dpdr,dpds                           
        USE allocation, ONLY: alloc_basis_arrays             
        USE messenger2, ONLY: myrank        
        USE read_dginp, ONLY: p

        IMPLICIT NONE
        INTEGER :: i,j,pt,et,dof,ndf
        REAL(rp) :: qint
        REAL(rp) :: mm(mndof,mndof)
  
        CALL alloc_basis_arrays()
      
        DO et = 1,nel_type        
          

          CALL element_basis(et,p,ndf,nqpta(et),qpta(:,1,et),qpta(:,2,et),phia(:,:,et),dpdr(:,:,et),dpds(:,:,et))         
          
!           PRINT "(A)", 'Basis functions at quadrature points'
!           DO i = 1,ndf
!             PRINT "(100(F10.3))", (phia(i,j,et),j=1,nqpta(et))
!           ENDDO
!           PRINT "(A)", ' '            
            
          
          ! Compute mass matrix (constant*indentity matrix)
          DO i = 1,ndf
            DO j = 1,ndf
              mm(i,j) = 0d0
              DO pt = 1,nqpta(et)
                mm(i,j) = mm(i,j) + wpta(pt,et)*phia(i,pt,et)*phia(j,pt,et)
              ENDDO
            ENDDO
          ENDDO
        
          
          
        ENDDO                
        
        CALL modal2nodal()
        
        IF (myrank == 0) THEN
          DO et = 1,nel_type
          
            PRINT "(A)", "---------------------------------------------"
            PRINT "(A)", "         Basis Function Information          "
            PRINT "(A)", "---------------------------------------------"
            PRINT "(A)", " "

            PRINT "(A,I5)", "Polynomial order:",p           
          
            PRINT "(A,I5)", "Number of degrees of freedom:",ndf
            PRINT "(A)", " "        

            PRINT "(A)", 'Mass matrix'
            DO i = 1,ndf
              PRINT "(100(F10.3))", (mm(i,j),j=1,ndof(et))
            ENDDO
            PRINT "(A)", ' '

!             PRINT "(A)", 'Basis functions at quadrature points'
!             DO i = 1,ndf
!               PRINT "(100(F10.3))", (phia(i,j,et),j=1,nqpta(et))
!             ENDDO
!             PRINT "(A)", ' '             
            
          ENDDO
        ENDIF


        RETURN
      END SUBROUTINE area_basis
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE edge_basis()

        USE globals, ONLY: rp,ndof,nverts,mndof,nqpte,mnqpte,qpte,wpte,phie,phie_int,nel_type
        USE read_dginp, ONLY: p

        IMPLICIT NONE

        INTEGER :: dof,pt,et,i
        INTEGER :: tpts,ndf
        REAL(rp) :: phi(mndof*4*mnqpte)       

        DO et = 1,nel_type
          
          tpts = nverts(et)*nqpte(et)
          CALL element_basis(et,p,ndf,tpts,qpte(:,1,et),qpte(:,2,et),phie(:,:,et))
          
          DO pt = 1,tpts
            DO dof = 1,ndf
              phie_int(dof,pt,et) = phie(dof,pt,et)*wpte(pt,et)
            ENDDO
          ENDDO
          
        ENDDO        


        RETURN
      END SUBROUTINE edge_basis      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE element_basis(et,p,ndfs,npts,r,s,phi,dpdr,dpds)
      
      USE globals, ONLY: rp
      
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
        
        USE globals, ONLY: rp
        
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: p,npts
        INTEGER, INTENT(OUT) :: ndfs
        REAL(rp), INTENT(IN) :: r(npts),s(npts) 
        REAL(rp), DIMENSION(:,:), INTENT(OUT) :: phi
        REAL(rp), DIMENSION(:,:), OPTIONAL, INTENT(OUT) :: dpdr,dpds         
        
        INTEGER :: m,i,j,pt
        INTEGER :: calc_deriv        
        REAL(rp) :: dpda,dpdb,dadr,dads,ii       
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
!               phi(m,pt) = sqrt(2d0)*Pi(pt)*Pj(pt)*(1d0-b(pt))**i
              phi(m,pt) = 2d0*Pi(pt)*Pj(pt)*(1d0-b(pt))**i
            ENDDO

            IF (calc_deriv == 1) THEN
              CALL djacobi(0    ,0,i,a,npts,dPi)
              CALL djacobi(2*i+1,0,j,b,npts,dPj)                          

              ! Calculate derivative values
              DO pt = 1,npts
                dadr = 2d0/(1d0-s(pt))
                dads = 2d0*(1d0+r(pt))/(1d0-s(pt))**2d0
              
!                 dpda = sqrt(2d0)*dPi(pt)*Pj(pt)*(1d0-b(pt))**ii
!                 dpdb = sqrt(2d0)*Pi(pt)*(dPj(pt)*(1d0-b(pt))**ii - ii*(1d0-b(pt))**(ii-1d0)*Pj(pt))
                dpda = 2d0*dPi(pt)*Pj(pt)*(1d0-b(pt))**ii
                dpdb = 2d0*Pi(pt)*(dPj(pt)*(1d0-b(pt))**ii - ii*(1d0-b(pt))**(ii-1d0)*Pj(pt))
              
                dpdr(m,pt) = dpda*dadr
                dpds(m,pt) = dpda*dads + dpdb
              ENDDO
            ENDIF

          ENDDO
        ENDDO  
        
        END SUBROUTINE tri_basis        
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
        
        SUBROUTINE quad_basis(p,ndfs,npts,r,s,phi,dpdr,dpds)
        
        USE globals, ONLY: rp
        
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

        USE globals, ONLY: rp

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

        USE globals, ONLY: rp

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

      SUBROUTINE linear(phil)
      
        USE globals, ONLY: rp,nqpta,qpta

        IMPLICIT NONE
        REAL(rp),INTENT(OUT) :: phil(3,nqpta(1))
        INTEGER :: pt
        REAL(rp) :: r,s
        
        DO pt = 1,nqpta(1)
          r = qpta(pt,1,1)
          s = qpta(pt,2,1)

          phil(1,pt) = -.5d0*(r+s)
          phil(2,pt) = .5d0*(r+1d0)
          phil(3,pt) = .5d0*(s+1d0)  
          
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

      SUBROUTINE element_nodes(et,space,p,n,r,s)
      
      USE globals, ONLY: rp
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: et,space,p
      INTEGER, INTENT(OUT) :: n
      REAL(rp), DIMENSION(:), INTENT(OUT) :: r,s
      
      IF (mod(et,2) == 1) THEN
        CALL tri_nodes(1,p,n,r,s)
      ELSE IF  (mod(et,2) == 0) THEN
        CALL quad_nodes(1,p,n,r,s)
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

      USE globals, ONLY: rp
      USE messenger2, ONLY: myrank

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
      
!       open(unit=60,file=DIRNAME//'/'//'tri.d')

!       IF (myrank == 0) THEN
!         print*,' '
!         print*, 'Triangle Points'
!         do i = 1,n
! !           write(60,61) r(i),s(i),x(i),y(i)
!           print 21, r(i),s(i)
!         enddo
!         print*, ' '
!       ENDIF

 21   format(28(f10.4,1x))
 61   format(16000(e24.17,1x))      

      RETURN
      END SUBROUTINE
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      

      SUBROUTINE quad_nodes(space,p,n,r,s)
      
      USE globals, ONLY: rp
      USE messenger2, ONLY: myrank
      
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
      
!       IF (myrank == 0) THEN           
!         PRINT*, ' '
!         PRINT*, 'Quadrilateral Points'
!         DO m = 1,n
!           PRINT("(2(f10.4))"), r(m),s(m)
!         ENDDO         
!       ENDIF
      
      RETURN
      END SUBROUTINE


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine lglpts(n,r)
!c     Returns the n+1 nth order Legendre-Gauss-Lebotto points on an iterval of (-1,1) in r.
!c     n (input) : order 
!c     r(n+1) (output) : Legendre-Gauss-Lebotto points

      use globals, only: rp

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

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

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

      use globals, only: rp

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

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	
      subroutine xytors(np,x,y,r,s)
      
      use globals, only: rp      
      
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

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      END MODULE basis
