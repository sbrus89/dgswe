      MODULE evaluate
      
      IMPLICIT NONE
      
      
      
      CONTAINS
      
      
      SUBROUTINE vandermonde()
      
      USE globals, ONLY: pres,ctp,np,nnds,mnnds,nel_type,V,ipiv
      USE basis, ONLY: tri_nodes,tri_basis,quad_nodes,quad_basis
      
      IMPLICIT NONE
      INTEGER :: et,pt,dof,i,n
      REAL(pres) :: r(mnnds),s(mnnds)
      REAL(pres) :: phi(mnnds*mnnds)
      INTEGER :: info      
      
      ALLOCATE(V(mnnds,mnnds,nel_type))
      ALLOCATE(ipiv(mnnds,nel_type))
      
      DO et = 1,nel_type
        n = nnds(et)
        IF (mod(et,2) == 1) THEN
          PRINT*, np(et)
          CALL tri_nodes(1,np(et),n,r,s)
          CALL tri_basis(np(et),n,n,r,s,phi)       
        ELSE IF (mod(et,2) == 0) THEN
          CALL quad_nodes(1,np(et),n,r,s)
          CALL quad_basis(np(et),n,n,r,s,phi)
        ENDIF
        
        DO pt = 1,n
          DO dof = 1,n
            i = (dof-1)*n + pt
            V(dof,pt,et) = phi(i)
          ENDDO
        ENDDO
        
        CALL DGETRF(n,n,V(1:n,1:n,et),n,ipiv(1:n,et),info)
!         DO pt = 1,n
!             PRINT("(100(e15.5))"), (V(dof,pt,et), dof = 1,n)
!         ENDDO        
!         PRINT*, " "
        
      ENDDO
      
      
      
      RETURN      
      END SUBROUTINE vandermonde
 
 
      SUBROUTINE newton(x,y,eln,r,s)

      USE globals, ONLY: pres,el_type,np,nnds,mnnds,V,elxy,ipiv
      USE basis, ONLY: tri_basis,quad_basis

      IMPLICIT NONE
      INTEGER :: it,eln,et,p,n,i
      INTEGER :: info
      INTEGER :: maxit
      REAL(pres) :: tol
      REAL(pres) :: x,y
      REAL(pres) :: r(1),s(1)
      REAL(pres) :: f,g
      REAL(pres) :: dfdr,dfds,dgdr,dgds,jac
      REAL(pres) :: phi(mnnds),dpdr(mnnds),dpds(mnnds)
      REAL(pres) :: l(mnnds,1),dldr(mnnds,1),dlds(mnnds,1)
        
      tol = 1d-12
      maxit = 10000
      info = 0
      
      et = el_type(eln)
      p = np(et)  
      n = nnds(et)
        
      IF (mod(et,2) == 1) THEN
        r(1) = -1d0/3d0
        s(1) = -1d0/3d0
      ELSE IF (mod(et,2) == 0) THEN
        r(1) = 1d0
        s(1) = 1d0
      ENDIF

      DO it = 1,maxit
           
        IF (mod(et,2) == 1) THEN
          CALL tri_basis(p,n,1,r,s,phi,dpdr,dpds)
        ELSE IF (mod(et,2) == 0) THEN
          CALL quad_basis(p,n,1,r,s,phi,dpdr,dpds)
        ENDIF
        
        DO i = 1,n
          l(i,1) = phi(i)
          dldr(i,1) = dpdr(i)
          dlds(i,1) = dpds(i)
        ENDDO
        
        CALL DGETRS("N",n,1,V(1:n,1:n,et),n,ipiv(1:n,et),l,n,info)
!         IF (info /= 0 ) PRINT*, "LAPACK ERROR"  
        CALL DGETRS("N",n,1,V(1:n,1:n,et),n,ipiv(1:n,et),dldr,n,info)
!         IF (info /= 0 ) PRINT*, "LAPACK ERROR"      
        CALL DGETRS("N",n,1,V(1:n,1:n,et),n,ipiv(1:n,et),dlds,n,info)
!         IF (info /= 0 ) PRINT*, "LAPACK ERROR"      
        
        dfdr = 0d0
        dfds = 0d0
        dgdr = 0d0
        dgds = 0d0
        f = 0d0
        g = 0d0        
        
        DO i = 1,n
          dfdr = dfdr + dldr(i,1)*elxy(i,eln,1)
          dfds = dfds + dlds(i,1)*elxy(i,eln,1)
          dgdr = dgdr + dldr(i,1)*elxy(i,eln,2)
          dgds = dgds + dlds(i,1)*elxy(i,eln,2)
          f = f + l(i,1)*elxy(i,eln,1)
          g = g + l(i,1)*elxy(i,eln,2)
        ENDDO
        
        jac = dfdr*dgds - dgdr*dfds
        
        f = f - x
        g = g - y
        
        r(1) = r(1) - (1d0/jac)*( dgds*f - dfds*g)
        s(1) = s(1) - (1d0/jac)*(-dgdr*f + dfdr*g)         
        
        IF (ABS(f) < tol .AND. ABS(g) < tol) THEN
          EXIT
        ENDIF        
       
        
      ENDDO
      
      IF (it >= maxit) THEN
        PRINT("(A)"), "   MAX ITERATIONS EXCEEDED"
        PRINT("(2(A,F20.15))"), "   r = ",r(1), "   s = ", s(1)
      ELSE       
        PRINT("(A,I7)"), "   iterations: ",it
        PRINT("(2(A,F20.15))"), "   r = ",r(1), "   s = ", s(1)
      ENDIF

      RETURN
      END SUBROUTINE newton
      

      
      
      END MODULE evaluate