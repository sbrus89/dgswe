      MODULE basis

      CONTAINS

      SUBROUTINE area_basis()

        USE globals, ONLY: pres,nel_type,nqpta,nqpte,wpta,qpta,p,ndof, &
                           phia,phia_int,phia_int_init,phil,dpdx,dpdy,dpdx_init,dpdy_init,dpdr,dpds, &
                           ect,xy,ne,area

        IMPLICIT NONE
        INTEGER :: i,j,pt,et
        INTEGER :: alloc_status
        INTEGER :: mndof,mnqpta
        REAL(pres) :: ml2(3,ndof(1)),mml(3,3)
        REAL(pres) :: qint
        
        mndof = maxval(ndof)
        mnqpta = maxval(nqpta)

        ALLOCATE(phia(mndof,mnqpta,nel_type),phia_int(ne,mndof*mnqpta),phia_int_init(ne,mndof*mnqpta),STAT = alloc_status)
        IF(alloc_status /= 0) THEN
          PRINT*, 'Allocation error: phia,phia_int'
        ENDIF 
        
        ALLOCATE(dpdr(mndof,mnqpta,nel_type),dpds(mndof,mnqpta,nel_type),STAT = alloc_status)
        IF(alloc_status /= 0) THEN
          PRINT*, 'Allocation error: phia,phia_int'
        ENDIF         
        
        ALLOCATE(phil(3,mnqpta,nel_type),STAT = alloc_status)
        IF(alloc_status /= 0) THEN
          PRINT*, 'Allocation error: phil'
        ENDIF

        ALLOCATE(dpdx(ne,mndof*mnqpta),dpdy(ne,mndof*mnqpta),STAT = alloc_status)
        IF(alloc_status /= 0) THEN
          PRINT*, 'Allocation error: dpdx,dpdy'
        ENDIF 
        
        ALLOCATE(dpdx_init(ne,mndof*mnqpta),dpdy_init(ne,mndof*mnqpta),STAT = alloc_status)
        IF(alloc_status /= 0) THEN
          PRINT*, 'Allocation error: dpdx_init,dpdy_init'
        ENDIF         

        PRINT "(A)", "---------------------------------------------"
        PRINT "(A)", "         Basis Function Information          "
        PRINT "(A)", "---------------------------------------------"
        PRINT "(A)", " "

        PRINT "(A,I5)", "Polynomial order:",p 

      
        DO et = 1,nel_type
        
          IF (mod(et,2) == 1) THEN
            CALL tri_area_basis(et,p,ndof(et),nqpta(et))
          ELSE IF (mod(et,2) == 0) THEN
            CALL quad_area_basis(et,p,ndof(et),nqpta(et))
          ENDIF
          
        ENDDO


        ! Calculate linear nodal basis functions
        CALL linear(phil)

        ! Compute RHS L2 projection matrix
        DO i = 1,3
          DO j = 1,ndof(1)
            qint = 0d0
            DO pt = 1,nqpta(1)
              qint = qint + wpta(pt,1)*phil(i,pt,1)*phia(j,pt,1)
            ENDDO
            ml2(i,j) = qint
          ENDDO
        ENDDO

        ! Compute linear mass matrix
        DO i = 1,3
          DO j = 1,3
            qint = 0d0
            DO pt = 1,nqpta(1)
              qint = qint + wpta(pt,1)*phil(i,pt,1)*phil(j,pt,1)
            ENDDO
            mml(i,j) = qint
          ENDDO
        ENDDO
        
        CALL modal2nodal()

        OPEN(unit=10,file='../output/projection.d')

        WRITE(10,*) ndof(1)
        DO i = 1,3
          WRITE(10,"(160(e24.17,1x))") (ml2(i,j),j=1,ndof(1))
        ENDDO

        DO i = 1,3
          WRITE(10,"(3(e24.17,1x))") (mml(i,j),j=1,3)
        ENDDO

        RETURN
      END SUBROUTINE area_basis
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE edge_basis()

        USE globals, ONLY: pres,ndof,nqpte,qpte,wpte,p,phi,phie,phie_int,nqpta,nel_type

        IMPLICIT NONE
        
        INTEGER :: alloc_status
        INTEGER :: mndof,mnqpte
        INTEGER :: et,i

        mndof = maxval(ndof)
        mnqpte = maxval(nqpte)        

        ALLOCATE(phie(mndof,4*mnqpte,nel_type),phie_int(mndof,4*mnqpte,nel_type),STAT = alloc_status)
        IF(alloc_status /= 0) THEN
          PRINT*, 'Allocation error: phie,phie_int'
        ENDIF 
        

        DO et = 1,nel_type
        
          IF (mod(et,2) == 1) THEN
            CALL tri_edge_basis(et,p,ndof(et),nqpte(et))
          ELSE IF (mod(et,2) == 0) THEN
            CALL quad_edge_basis(et,p,ndof(et),nqpte(et))
          ENDIF
          
        ENDDO        


        RETURN
      END SUBROUTINE edge_basis      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE tri_area_basis(et,p,ndof,nqpta)
        
        USE globals, ONLY: pres,qpta,wpta,phia,phia_int,dpdr,dpds
        
        IMPLICIT NONE
        
        INTEGER :: m,i,j,pt,et
        INTEGER :: p,nqpta
        INTEGER :: ndof
        REAL(pres) :: dpda,dpdb,dadr,dads,ii       
        REAL(pres) :: r(nqpta),s(nqpta),a(nqpta),b(nqpta)
        REAL(pres) :: Pi(nqpta),Pj(nqpta)
        REAL(pres) :: dPi(nqpta),dPj(nqpta)
        REAL(pres) :: mm(ndof,ndof)
        REAL(pres) :: ml2(3,ndof),mml(3,3)
        REAL(pres) :: qint
        
      
        ! Change quadrature points from r,s (master element) to a,b extended coordinates
        DO pt = 1,nqpta
          r(pt) = qpta(pt,1,et)
          s(pt) = qpta(pt,2,et)
          
          a(pt) = 2d0*(1d0+r(pt))/(1d0-s(pt))-1d0 
          b(pt) = s(pt)
        ENDDO
        
        ! Calculate basis function values and derivative values at area quadrature points
        m = 0
        DO i = 0,p
          DO j = 0,p-i

            m = m+1

            Pi = 0d0
            dPi = 0d0
            Pj = 0d0
            dPj = 0d0

            CALL jacobi(0,0,i,a,nqpta,Pi)
            CALL djacobi(0,0,i,a,nqpta,dPi)
            CALL jacobi(2*i+1,0,j,b,nqpta,Pj)
            CALL djacobi(2*i+1,0,j,b,nqpta,dPj)

            ! Calculate function values
            DO pt = 1,nqpta 
!               phia(m,pt) = sqrt(2d0)*Pi(pt)*Pj(pt)*(1d0-b(pt))**i
              phia(m,pt,et) = 2d0*Pi(pt)*Pj(pt)*(1d0-b(pt))**i
            ENDDO

            ii = real(i,pres)

            ! Calculate derivative values
            DO pt = 1,nqpta
              dadr = 2d0/(1d0-s(pt))
!               dpda = sqrt(2d0)*dPi(pt)*Pj(pt)*(1d0-b(pt))**ii
              dpda = 2d0*dPi(pt)*Pj(pt)*(1d0-b(pt))**ii
              dpdr(m,pt,et) = dpda*dadr

              dads = 2d0*(1d0+r(pt))/(1d0-s(pt))**2d0
!               dpdb = sqrt(2d0)*Pi(pt)*(dPj(pt)*(1d0-b(pt))**ii - ii*(1d0-b(pt))**(ii-1d0)*Pj(pt))
              dpdb = 2d0*Pi(pt)*(dPj(pt)*(1d0-b(pt))**ii - ii*(1d0-b(pt))**(ii-1d0)*Pj(pt))
              dpds(m,pt,et) = dpda*dads + dpdb
            ENDDO

          ENDDO
        ENDDO

        ! Compute mass matrix (indentity matrix)
        DO i = 1,ndof
          DO j = 1,ndof
            qint = 0d0
            DO pt = 1,nqpta
              qint = qint + wpta(pt,et)*phia(i,pt,et)*phia(j,pt,et)
            ENDDO
            mm(i,j) = qint
          ENDDO
        ENDDO
        
        PRINT "(A,I5)", "Number of degrees of freedom:",ndof
        PRINT "(A)", " "        

        PRINT "(A)", 'Mass matrix'
        DO i = 1,ndof
          PRINT "(100(F10.3))", (mm(i,j),j=1,ndof)
        ENDDO
        PRINT "(A)", ' '

        PRINT "(A)", 'Basis functions at quadrature points'
        DO i = 1,m
          PRINT "(100(F10.3))", (phia(i,j,et),j=1,nqpta)
        ENDDO
        PRINT "(A)", ' '     
        
        END SUBROUTINE tri_area_basis
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE quad_area_basis(et,p,ndof,nqpta)
        
        USE globals, ONLY: pres,qpta,wpta,phia,phia_int,dpdr,dpds
        
        IMPLICIT NONE
        
        INTEGER :: m,i,j,pt,et
        INTEGER :: p,nqpta
        INTEGER :: ndof
        REAL(pres) :: dpda,dpdb,dadr,dads      
        REAL(pres) :: r(nqpta),s(nqpta),a(nqpta),b(nqpta)
        REAL(pres) :: Pi(nqpta),Pj(nqpta)
        REAL(pres) :: dPi(nqpta),dPj(nqpta)
        REAL(pres) :: mm(ndof,ndof)
        REAL(pres) :: ml2(3,ndof),mml(3,3)
        REAL(pres) :: qint
        
      
        ! Change quadrature points from r,s (master element) to a,b extended coordinates
        DO pt = 1,nqpta
          r(pt) = qpta(pt,1,et)
          s(pt) = qpta(pt,2,et)
        ENDDO
        
        ! Calculate basis function values and derivative values at area quadrature points
        m = 0
        DO i = 0,p
          DO j = 0,p

            m = m+1

            Pi = 0d0
            dPi = 0d0
            Pj = 0d0
            dPj = 0d0

            CALL jacobi(0,0,i,r,nqpta,Pi)
            CALL djacobi(0,0,i,r,nqpta,dPi)
            CALL jacobi(0,0,j,s,nqpta,Pj)
            CALL djacobi(0,0,j,s,nqpta,dPj)

            ! Calculate function values
            DO pt = 1,nqpta 
              phia(m,pt,et) = 2d0*Pi(pt)*Pj(pt)
            ENDDO

            ! Calculate derivative values
            DO pt = 1,nqpta
              dpdr(m,pt,et) = 2d0*dPi(pt)*Pj(pt)
              dpds(m,pt,et) = 2d0*Pi(pt)*dPj(pt)
            ENDDO

          ENDDO
        ENDDO

        ! Compute mass matrix (indentity matrix)
        DO i = 1,ndof
          DO j = 1,ndof
            qint = 0d0
            DO pt = 1,nqpta
              qint = qint + wpta(pt,et)*phia(i,pt,et)*phia(j,pt,et)
            ENDDO
            mm(i,j) = qint
          ENDDO
        ENDDO
        
        PRINT "(A,I5)", "Number of degrees of freedom:",ndof
        PRINT "(A)", " "        

        PRINT "(A)", 'Mass matrix'
        DO i = 1,ndof
          PRINT "(100(F10.3))", (mm(i,j),j=1,ndof)
        ENDDO
        PRINT "(A)", ' '

        PRINT "(A)", 'Basis functions at quadrature points'
        DO i = 1,m
          PRINT "(100(F10.3))", (phia(i,j,et),j=1,nqpta)
        ENDDO
        PRINT "(A)", ' '     
        
        END SUBROUTINE quad_area_basis        
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
      SUBROUTINE tri_edge_basis(et,p,ndof,nqpte)
      
        USE globals, ONLY: pres,wpte,qpte,phie,phie_int
      
        IMPLICIT NONE
        INTEGER :: led,pt,m,i,j,ind
        INTEGER :: ndof,nqpte,et,p
        REAL(pres) :: r(nqpte),s(nqpte),a(nqpte),b(nqpte)
        REAL(pres) :: Pi(nqpte),Pj(nqpte)   
      
        DO led = 1,3

          ! Determine r,s coordinates of quadrature points along local edge in master element
          SELECT CASE(led)
            CASE(1)
              DO pt = 1,nqpte
                r(pt) = -qpte(pt,et)
                s(pt) =  qpte(pt,et)
!                 PRINT*, r(pt),s(pt)
              ENDDO
            CASE(2)
              DO pt = 1,nqpte
                r(pt) = -1d0
                s(pt) = -qpte(pt,et)
!                 PRINT*, r(pt),s(pt)
              ENDDO 
            CASE(3)
              DO pt = 1,nqpte
                r(pt) = qpte(pt,et)
                s(pt) = -1d0
!                 PRINT*, r(pt),s(pt)
              ENDDO 
          END SELECT

          ! Change quadrature points from r,s (master element) to a,b extended coordinates
          DO pt = 1,nqpte
            a(pt) = 2d0*(1d0+r(pt))/(1d0-s(pt)) - 1d0
            b(pt) = s(pt)
          ENDDO

          ! Calculate basis function values edge quadrature points
          m = 0
          DO i = 0,p
            DO j = 0,p-i

              m = m+1

              Pi = 0d0
              Pj = 0d0
              CALL jacobi(0,0,i,a,nqpte,Pi)
              CALL jacobi(2*i+1,0,j,b,nqpte,Pj)
            
              ! Calculate function values
              DO pt = 1,nqpte 
                ind = (led-1)*nqpte+pt
!                 phie(m,ind) = sqrt(2d0)*Pi(pt)*Pj(pt)*(1d0-b(pt))**i
!                 phie_int(m,ind) = phie(m,ind)*wpte(pt)
                phie(m,ind,et) = 2d0*Pi(pt)*Pj(pt)*(1d0-b(pt))**i
                phie_int(m,ind,et) = phie(m,ind,et)*wpte(pt,et)
              ENDDO

            ENDDO
          ENDDO

        ENDDO
      
      
      END SUBROUTINE tri_edge_basis
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
      SUBROUTINE quad_edge_basis(et,p,ndof,nqpte)
      
        USE globals, ONLY: pres,wpte,qpte,phie,phie_int
      
        IMPLICIT NONE
        INTEGER :: led,pt,m,i,j,ind
        INTEGER :: ndof,nqpte,et,p
        REAL(pres) :: r(nqpte),s(nqpte)
        REAL(pres) :: Pi(nqpte),Pj(nqpte)   
      
        DO led = 1,4

          ! Determine r,s coordinates of quadrature points along local edge in master element
          SELECT CASE(led)
            CASE(1)
              DO pt = 1,nqpte
                r(pt) = 1d0
                s(pt) = qpte(pt,et)
              ENDDO
            CASE(2)
              DO pt = 1,nqpte
                r(pt) = -qpte(pt,et)
                s(pt) = 1d0
              ENDDO 
            CASE(3)
              DO pt = 1,nqpte
                r(pt) = -1d0
                s(pt) = -qpte(pt,et)
              ENDDO 
            CASE(4) 
              DO pt = 1,nqpte
                r(pt) = qpte(pt,et)
                s(pt) = -1d0
              ENDDO
          END SELECT

          ! Calculate basis function values edge quadrature points
          m = 0
          DO i = 0,p
            DO j = 0,p

              m = m+1

              Pi = 0d0
              Pj = 0d0
              CALL jacobi(0,0,i,r,nqpte,Pi)
              CALL jacobi(0,0,j,s,nqpte,Pj)
            
              ! Calculate function values
              DO pt = 1,nqpte 
                ind = (led-1)*nqpte+pt

                phie(m,ind,et) = 2d0*Pi(pt)*Pj(pt)
                phie_int(m,ind,et) = phie(m,ind,et)*wpte(pt,et)
              ENDDO

            ENDDO
          ENDDO
          

        ENDDO
      
      
      END SUBROUTINE quad_edge_basis      
      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE jacobi(alpha_i,beta_i,deg,x,npts,v)

        USE globals, ONLY: pres,ndof

        IMPLICIT NONE
        INTEGER, INTENT(IN) :: alpha_i,beta_i,deg,npts
        REAL(pres), INTENT(IN) :: x(npts)
        REAL(pres), INTENT(OUT) :: v(npts)
        INTEGER :: i,j,np1
        REAL(pres) :: pnm1(npts),pn(npts),pnp1(npts)
        REAL(pres) :: alpha,beta,an,anp1,bn,n

        alpha = real(alpha_i,pres)
        beta = real(beta_i,pres)

        ! Calculate constant P0
        DO i = 1,npts
          pnm1(i) = sqrt(2d0**(-alpha-beta-1d0)*real(fact(alpha_i+beta_i+1),pres) &
                         /(real(fact(alpha_i),pres)*real(fact(beta_i),pres))) 
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
          n = real(np1,pres)
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

        USE globals, ONLY: pres,ndof

        IMPLICIT NONE
        INTEGER, INTENT(IN) :: alpha_i,beta_i,deg_i,npts
        REAL(pres), INTENT(IN) :: x(npts)
        REAL(pres), INTENT(OUT) :: dP(npts)
        REAL(pres):: v(npts),deg,alpha,beta
        INTEGER :: i,pt

        deg = real(deg_i,pres)
        alpha = real(alpha_i,pres)
        beta = real(beta_i,pres)
        
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
      
        USE globals, ONLY: pres,nqpta,qpta

        IMPLICIT NONE
        REAL(pres),INTENT(OUT) :: phil(3,nqpta(1))
        INTEGER :: pt
        REAL(pres) :: r,s
        
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
      
      subroutine tri_nodes(space,p,np,r,s)
!c     Calculates nodal set for well-behaved interpolation
!c     p (input) : the order of the basis
!c     np (input): the number of points in the triangle
!c     r(np) (output) ; nodal r-coordinates for master element
!c     s(np) (output) : nodal s-coordinates for master element 

      use globals, only: pres

      implicit none
      integer p,np,i,j,m,space 
      real(pres) ii,jj,a,dx,tol
      real(pres) r(np),s(np),x(np),y(np),l1(np),l2(np),l3(np)
      real(pres) aopt(15),b1(np),b2(np),b3(np)
      real(pres) wx(np),wy(np),w1(np),w2(np),w3(np)
      real(pres) w1e(np),w2e(np),w3e(np)
      real(pres) w1mat(np,p+1),w2mat(np,p+1),w3mat(np,p+1)
      real(pres) var1(np),var2(np),var3(np)
      real(pres) lgl(p+1),xeq(p+1)
      !c DW
      

!c     Define optimal alpha values
      data aopt/0d0,0d0,1.4152d0,0.1001d0,0.2751d0,0.9800d0,1.0999d0 &
               ,1.2832d0,1.3648d0,1.4773d0,1.4959d0,1.5743d0,1.5770d0 &
               ,1.6223d0,1.6258d0/

      if(p.lt.16) then
        a = aopt(p)
      else
        a = 5d0/3d0
      endif

!c     Create equally distributed nodes on equalateral triangle

      m = 0

! Edge 3 points
      i = 1
      ii = dble(i)
      do j = 1,p+1
        jj = dble(j)
        m = m+1
        l1(m) = (ii-1d0)/dble(p)
        l3(m) = (jj-1d0)/dble(p)
      enddo

! Edge 1 points
      do i = 2,p+1
        ii = dble(i)
        j = p+2-i 
        jj = dble(j)
        m = m+1
        l1(m) = (ii-1d0)/dble(p)
        l3(m) = (jj-1d0)/dble(p)
      enddo


! Edge 2 points
      do i = p,2,-1
        ii = dble(i)
        j = 1
        jj = dble(j)
        m = m+1
        l1(m) = (ii-1d0)/dble(p)
        l3(m) = (jj-1d0)/dble(p)
      enddo

! Internal points
      do i = 2,p
        ii = dble(i)
        do j = 2,p+1-i 
          jj = dble(j)
          m = m+1
          l1(m) = (ii-1d0)/dble(p)
          l3(m) = (jj-1d0)/dble(p)
        enddo
      enddo

!       do i = 1,p+1
!         ii = dble(i)
!         do j = 1,p+2-i 
!           jj = dble(j)
!           m = m+1
!           l1(m) = (ii-1d0)/dble(p)
!           l3(m) = (jj-1d0)/dble(p)
!         enddo
!       enddo

      do i = 1,m
        l2(i) = 1d0-l1(i)-l3(i)
        x(i) = -l2(i)+l3(i)
        y(i) = (-l2(i)-l3(i)+2d0*l1(i))/dsqrt(3d0)
      enddo

      if(space.eq.0)then
        call xytors(np,x,y,r,s) ;
        return
      endif

!c     Compute blending function
      do i = 1,np
        b1(i) = 4d0*l2(i)*l3(i)
        b2(i) = 4d0*l3(i)*l1(i)
        b3(i) = 4d0*l2(i)*l1(i)
      enddo

!c     Compute 1-D Legendre-Gauss-Lobotto points on (-1,1)
      call lglpts(p,lgl)

      do i = 1,p+1
        write(40,*) lgl(i)
      enddo

!c     Compute equadistant points on (-1,1)
      dx = 2d0/dble(p)
      xeq(1)=-1d0
      do i = 2,p+1
        xeq(i) = xeq(i-1)+dx
      enddo
      xeq(p+1) = 1d0

!c     Compute warping function arugments
      do i = 1,np
        w1e(i) = l3(i)-l2(i)
        w2e(i) = l1(i)-l3(i)
        w3e(i) = l2(i)-l1(i)
      enddo

!c     Elvaluate Lagrange polynomials based on equa-spaced nodes
      call lagrange(p+1,np,xeq,w1e,w1mat)
      call lagrange(p+1,np,xeq,w2e,w2mat)
      call lagrange(p+1,np,xeq,w3e,w3mat)

      do j = 1,np
        write(50,41) (w1mat(j,i), i = 1,p+1)
      enddo

 41   format(16000(e24.17,1x))

!c     Compute 1 dimensional mapping functions for edges

      tol = 1d0 - 1.0d-10

      do i = 1,np
        if(dabs(w1e(i)).lt.tol) then
          var1(i) = 1d0
        else
          var1(i) = 0d0
        endif
        if(dabs(w2e(i)).lt.tol) then
          var2(i) = 1d0
        else
          var2(i) = 0d0
        endif
        if(dabs(w3e(i)).lt.tol) then
          var3(i) = 1d0
        else
          var3(i) = 0d0
        endif
      enddo

      do j = 1,np
        w1(j) = 0d0
        w2(j) = 0d0
        w3(j) = 0d0
      enddo

      do j = 1,np
        do i = 1,p+1
         w1(j) = w1(j) + (lgl(i)-xeq(i))*w1mat(j,i)
         w2(j) = w2(j) + (lgl(i)-xeq(i))*w2mat(j,i)
         w3(j) = w3(j) + (lgl(i)-xeq(i))*w3mat(j,i)
        enddo
        w1(j) = w1(j)/(1d0-(var1(j)*w1e(j))**2)+w1(j)*(var1(j)-1d0)
        w2(j) = w2(j)/(1d0-(var2(j)*w2e(j))**2)+w2(j)*(var2(j)-1d0)
        w3(j) = w3(j)/(1d0-(var3(j)*w3e(j))**2)+w3(j)*(var3(j)-1d0)
      enddo

!c     Apply gerneralized warping function to equally distributed nodes on equalateral triangle
      do i = 1,np
        x(i) = x(i)+(1d0+(a*l1(i))**2d0)*b1(i)*w1(i) &
                   -.5d0*(1d0+(a*l2(i))**2d0)*b2(i)*w2(i) &
                   -.5d0*(1d0+(a*l3(i))**2d0)*b3(i)*w3(i)
        y(i) = y(i)+dsqrt(3d0)/2d0*(1d0+(a*l2(i))**2d0)*b2(i)*w2(i) &
                   -dsqrt(3d0)/2d0*(1d0+(a*l3(i))**2d0)*b3(i)*w3(i)
      enddo

      call xytors(np,x,y,r,s)
      
!       open(unit=60,file=DIRNAME//'/'//'tri.d')

      print*,' '
      print*, 'Triangle Points'
      do i = 1,np
!         write(60,61) r(i),s(i),x(i),y(i)
        print 21, r(i),s(i)
      enddo
      print*, ' '

 21   format(28(f10.4,1x))
 61   format(16000(e24.17,1x))      

      return
      end subroutine
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       

      SUBROUTINE quad_nodes(space,np,nnds,r,s)
      
      USE globals, ONLY: pres
      
      IMPLICIT NONE
      
      INTEGER :: space,np,nnds
      INTEGER :: i,k,j,n,nnp
      REAL(pres) :: xi(np+1)
      REAL(pres) :: r(nnds),s(nnds)
      
      
      ! Get 1-D LGL points 
      
      IF (space == 1) THEN
        CALL lglpts(np,xi)
      ELSE
        xi(1) = -1d0
        DO i = 1,np-1
          xi(i+1) = xi(i) + 2d0/real(np,pres)
        ENDDO
        xi(np+1) = 1d0
      ENDIF
        
      ! Do tensor product, ordering the nodes counter-clockwise

      ! Find number of loops around refence quad element, excluding middle points
      IF (np <= 2)THEN
        nnp = 1
      ELSE
        IF(mod(np,2) == 1) THEN
          nnp = np-1
        ELSE IF (mod(np,2) == 0) THEN
          nnp = np-2
        ENDIF
      ENDIF
        
      n = 1
      DO k = 1,nnp ! loop over number of loops
         
        ! Edge 4
        DO i = k,np+1 - (k-1)
          j = k
          r(n) = xi(i)
          s(n) = xi(j)
          
          n = n+1
        ENDDO
        
        ! Edge 1
        DO j = 2 + (k-1),np+1 - (k-1)
          i = np+1 - (k-1)
          r(n) = xi(i)
          s(n) = xi(j)
          
          n = n+1
        ENDDO

        ! Edge 2
        DO i = np - (k-1),1 + (k-1),-1
          j = np+1 - (k-1)
          r(n) = xi(i)
          s(n) = xi(j)
          
          n = n+1      
        ENDDO

        ! Edge 3
        DO j = np - (k-1),2 + (k-1),-1
          i = 1 + (k-1)
          r(n) = xi(i)
          s(n) = xi(j)

          n = n+1
        ENDDO    
          
      ENDDO
        
      ! middle point
      IF (mod(np+1,2) == 1) THEN
        i = np/2 + 1
        r(n) = xi(i)
        s(n) = xi(i)
          
      ENDIF   
      
      PRINT*, ' '
      PRINT*, 'Quadrilateral Points'
      DO n = 1,nnds
        PRINT("(2(f10.4))"), r(n),s(n)
      ENDDO         
      
      RETURN
      END SUBROUTINE


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine lglpts(n,r)
!c     Returns the n+1 nth order Legendre-Gauss-Lebotto points on an iterval of (-1,1) in r.
!c     n (input) : order 
!c     r(n+1) (output) : Legendre-Gauss-Lebotto points

      use globals, only: pres

      implicit none
      integer n,nn,i
      real(pres) ii,h
      integer info
      real(pres) r(n+1),d(n+1),e(n)
      real(pres) z(n+1,n+1),work(1,2*n-2)

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

      use globals, only: pres

      implicit none
      integer nn,ne,n,i,j
      real(pres) p,xn(nn),xe(ne),pmat(ne,nn)

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
      
      use globals, only: pres      
      
      implicit none
      integer np,i
      real(pres) x(np),y(np)
      real(pres) r(np),s(np)
      real(pres) l1(np),l2(np),l3(np)

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
