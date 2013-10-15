      MODULE basis

      CONTAINS

      SUBROUTINE area_basis()

        USE globals, ONLY: pres,nqpta,wpta,qpta,p,ndof,phia,phia_int,phil,dpdx,dpdy,ect,xy,ne,area

        IMPLICIT NONE
        INTEGER :: m,i,j,pt,el,dof,ind,jj
        INTEGER :: alloc_status
        REAL(pres) :: dpda,dpdb,dadr,dads,ii
        REAL(pres) :: dpdr(ndof,nqpta),dpds(ndof,nqpta)
        REAL(pres) :: x1,x2,x3,y1,y2,y3
        REAL(pres) :: r(nqpta),s(nqpta),a(nqpta),b(nqpta)
        REAL(pres) :: Pi(nqpta),Pj(nqpta)
        REAL(pres) :: dPi(nqpta),dPj(nqpta)
        REAL(pres) :: mm(ndof,ndof)
        REAL(pres) :: ml2(3,ndof),mml(3,3)
        REAL(pres) :: qint

        ALLOCATE(phia(ndof,nqpta),phia_int(ndof,nqpta),STAT = alloc_status)
        IF(alloc_status /= 0) THEN
          PRINT*, 'Allocation error: phia,phia_int'
        ENDIF 

        ALLOCATE(phil(3,nqpta),STAT = alloc_status)
        IF(alloc_status /= 0) THEN
          PRINT*, 'Allocation error: phil'
        ENDIF

        ALLOCATE(dpdx(ne,ndof*nqpta),dpdy(ne,ndof*nqpta),STAT = alloc_status)
        IF(alloc_status /= 0) THEN
          PRINT*, 'Allocation error: dpdx,dpdy'
        ENDIF 

      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", "         Basis Function Information          "
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", " "

      PRINT "(A,I5)", "Polynomial order:",p 
      PRINT "(A,I5)", "Number of degrees of freedom:",ndof
      PRINT "(A)", " "

        ! Change quadrature points from r,s (master element) to a,b extended coordinates
        DO pt = 1,nqpta
          r(pt) = qpta(pt,1)
          s(pt) = qpta(pt,2)
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
!               phia_int(m,pt) = phia(m,pt)*wpta(pt)
              phia(m,pt) = 2d0*Pi(pt)*Pj(pt)*(1d0-b(pt))**i
              phia_int(m,pt) = .5d0*phia(m,pt)*wpta(pt)
            ENDDO

            ii = real(i,pres)

            ! Calculate derivative values
            DO pt = 1,nqpta
              dadr = 2d0/(1d0-s(pt))
!               dpda = sqrt(2d0)*dPi(pt)*Pj(pt)*(1d0-b(pt))**ii
              dpda = 2d0*dPi(pt)*Pj(pt)*(1d0-b(pt))**ii
              dpdr(m,pt) = dpda*dadr

              dads = 2d0*(1d0+r(pt))/(1d0-s(pt))**2d0
!               dpdb = sqrt(2d0)*Pi(pt)*(dPj(pt)*(1d0-b(pt))**ii - ii*(1d0-b(pt))**(ii-1d0)*Pj(pt))
              dpdb = 2d0*Pi(pt)*(dPj(pt)*(1d0-b(pt))**ii - ii*(1d0-b(pt))**(ii-1d0)*Pj(pt))
              dpds(m,pt) = dpda*dads + dpdb
            ENDDO

          ENDDO
        ENDDO

        ! Compute mass matrix (indentity matrix)
        DO i = 1,ndof
          DO j = 1,ndof
            qint = 0d0
            DO pt = 1,nqpta
              qint = qint + wpta(pt)*phia(i,pt)*phia(j,pt)
            ENDDO
            mm(i,j) = qint
          ENDDO
        ENDDO

        PRINT "(A)", 'Mass matrix'
        DO i = 1,ndof
          PRINT "(100(F10.3))", (mm(i,j),j=1,ndof)
        ENDDO
        PRINT "(A)", ' '

        PRINT "(A)", 'Basis functions at quadrature points'
        DO i = 1,m
          PRINT "(100(F10.3))", (phia(i,j),j=1,nqpta)
        ENDDO
        PRINT "(A)", ' '

!             DO J = 0,P
!                DO I = 0,J
!                   JJ = J - I
! 
!                   print*, ((2.D0*I+1.D0)*(2.D0*JJ+2.D0*I+2.D0)/4.D0)
!                ENDDO
!             ENDDO

        ! Derivative coordinate transformation
        DO dof = 1,ndof
          DO pt = 1,nqpta
            ind = (dof-1)*nqpta + pt
            DO el = 1,ne

              x1 = xy(1,ect(1,el))
              x2 = xy(1,ect(2,el))
              x3 = xy(1,ect(3,el)) 

              y1 = xy(2,ect(1,el))
              y2 = xy(2,ect(2,el))
              y3 = xy(2,ect(3,el))

!               dpdx(el,ind) = wpta(pt)*(dpdr(dof,pt)*(y3 - y1) + dpds(dof,pt)*(y1 - y2))/area(el)
!               dpdy(el,ind) = wpta(pt)*(dpdr(dof,pt)*(x1 - x3) + dpds(dof,pt)*(x2 - x1))/area(el)

              dpdx(el,ind) = .5d0*wpta(pt)*(dpdr(dof,pt)*(y3 - y1) + dpds(dof,pt)*(y1 - y2))/area(el)
              dpdy(el,ind) = .5d0*wpta(pt)*(dpdr(dof,pt)*(x1 - x3) + dpds(dof,pt)*(x2 - x1))/area(el)
            ENDDO
          ENDDO
        ENDDO

        ! Calculate linear nodal basis functions
        CALL linear(phil)

        ! Compute RHS L2 projection matrix
        DO i = 1,3
          DO j = 1,ndof
            qint = 0d0
            DO pt = 1,nqpta
              qint = qint + wpta(pt)*phil(i,pt)*phia(j,pt)
            ENDDO
            ml2(i,j) = qint
          ENDDO
        ENDDO

        ! Compute linear mass matrix
        DO i = 1,3
          DO j = 1,3
            qint = 0d0
            DO pt = 1,nqpta
              qint = qint + wpta(pt)*phil(i,pt)*phil(j,pt)
            ENDDO
            mml(i,j) = qint
          ENDDO
        ENDDO

        OPEN(unit=10,file='../output/projection.d')

        WRITE(10,*) ndof
        DO i = 1,3
          WRITE(10,"(160(e24.17,1x))") (ml2(i,j),j=1,ndof)
        ENDDO

        DO i = 1,3
          WRITE(10,"(3(e24.17,1x))") (mml(i,j),j=1,3)
        ENDDO

        RETURN
      END SUBROUTINE area_basis

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE edge_basis()

        USE globals, ONLY: pres,ndof,nqpte,qpte,wpte,p,phie,phie_int

        IMPLICIT NONE
        INTEGER :: led,pt,m,i,j,ind
        INTEGER :: alloc_status
        REAL(pres) :: r(nqpte),s(nqpte),a(nqpte),b(nqpte)
        REAL(pres) :: Pi(nqpte),Pj(nqpte)

        ALLOCATE(phie(ndof,3*nqpte),phie_int(ndof,3*nqpte),STAT = alloc_status)
        IF(alloc_status /= 0) THEN
          PRINT*, 'Allocation error: phie,phie_int'
        ENDIF 

        DO led = 1,3

          ! Determine r,s coordinates of quadrature points along local edge in master element
          SELECT CASE(led)
            CASE(1)
              DO pt = 1,nqpte
                r(pt) = -qpte(pt)
                s(pt) =  qpte(pt)
!                 PRINT*, r(pt),s(pt)
              ENDDO
            CASE(2)
              DO pt = 1,nqpte
                r(pt) = -1d0
                s(pt) = -qpte(pt)
!                 PRINT*, r(pt),s(pt)
              ENDDO 
            CASE(3)
              DO pt = 1,nqpte
                r(pt) = qpte(pt)
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
                phie(m,ind) = 2d0*Pi(pt)*Pj(pt)*(1d0-b(pt))**i
                phie_int(m,ind) = .5d0*phie(m,ind)*wpte(pt)
              ENDDO

            ENDDO
          ENDDO

        ENDDO

        RETURN
      END SUBROUTINE edge_basis

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
        REAL(pres),INTENT(OUT) :: phil(3,nqpta)
        INTEGER :: pt
        REAL(pres) :: r,s
        
        DO pt = 1,nqpta
          r = qpta(pt,1)
          s = qpta(pt,2)

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

      END MODULE basis
