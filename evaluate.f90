      MODULE evaluate
      
      IMPLICIT NONE
      
      
      
      CONTAINS
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
            
      SUBROUTINE vandermonde(sol)
      
      USE globals, ONLY: pres,nel_type,solution
      USE basis, ONLY: tri_nodes,tri_basis,quad_nodes,quad_basis
      
      IMPLICIT NONE
      TYPE(solution) :: sol
      
      INTEGER :: et,pt,dof,i,n
      REAL(pres) :: r(sol%mnnds),s(sol%mnnds)
      REAL(pres) :: phi(sol%mnnds*sol%mnnds)
      INTEGER :: info    
      

      
      ALLOCATE(sol%V(sol%mnnds,sol%mnnds,nel_type))
      ALLOCATE(sol%ipiv(sol%mnnds,nel_type))
      
      DO et = 1,nel_type
        n = sol%nnds(et)
        IF (mod(et,2) == 1) THEN
          CALL tri_nodes(1,sol%np(et),n,r,s)
          CALL tri_basis(sol%np(et),n,n,r,s,phi)       
        ELSE IF (mod(et,2) == 0) THEN
          CALL quad_nodes(1,sol%np(et),n,r,s)
          CALL quad_basis(sol%np(et),n,n,r,s,phi)
        ENDIF
        
        DO pt = 1,n
          DO dof = 1,n
            i = (dof-1)*n + pt
            sol%V(dof,pt,et) = phi(i)
          ENDDO
        ENDDO
        
!         CALL DGETRF(n,n,V(1:n,1:n,et),n,ipiv(1:n,et),info)
        CALL DGETRF(n,n,sol%V(1,1,et),sol%mnnds,sol%ipiv(1,et),info)
        
!         DO pt = 1,n
!             PRINT("(100(e15.5))"), (sol%V(dof,pt,et), dof = 1,n)
!         ENDDO        
!         PRINT*, " "       
        
      ENDDO
      
      
      
      RETURN      
      END SUBROUTINE vandermonde
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
      
      SUBROUTINE re_vert(sol)
      
      USE globals, ONLY: pres,solution,nel_type
      USE basis, ONLY: tri_nodes, quad_nodes
      
      IMPLICIT NONE
      
      INTEGER :: et,pt
      INTEGER :: n,p
      REAL(pres) :: r(4),s(4)
      TYPE(solution) :: sol
      
      DO et = 1,nel_type
        IF (mod(et,2) == 1) THEN
          n = sol%nnds(1)
          p = sol%np(1)
          CALL tri_nodes(1,p,n,r,s)
        ELSE IF (mod(et,2) == 0) THEN
          n = sol%nnds(2)
          p = sol%np(2)
          CALL quad_nodes(1,p,n,r,s)
        ENDIF
        
        DO pt = 1,n
          sol%rsre(1,pt,et) = r(pt)
          sol%rsre(2,pt,et) = s(pt)
        ENDDO
        
!         PRINT("(4(F15.5))"), (sol%rsre(1,pt,et), pt = 1,n)
!         PRINT("(4(F15.5))"), (sol%rsre(2,pt,et), pt = 1,n)        
!         PRINT*, " "
        
      ENDDO      
      
      RETURN
      END SUBROUTINE re_vert

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     

      SUBROUTINE l_eval()
      
      USE globals, ONLY: pres,nel_type,nqpta,qpta,mnqpta,fine
      USE basis, ONLY: tri_basis,quad_basis
      
      IMPLICIT NONE
      
      INTEGER :: p,n,et,m,i,pt,mnnds,npt
      INTEGER :: info
      REAL(pres) :: r(mnqpta),s(mnqpta)
      REAL(pres) :: phi(fine%mnnds*mnqpta)
      
      mnnds = fine%mnnds
      
      ALLOCATE(fine%l(mnnds,mnqpta,nel_type))
      
      
      DO et = 1,nel_type
        p = fine%p
        n = fine%nnds(et)
        npt = nqpta(et)
        
        DO pt = 1,npt
          r(pt) = qpta(pt,1,et) 
          s(pt) = qpta(pt,2,et)
        ENDDO
        
        IF (mod(et,2) == 1) THEN
          CALL tri_basis(p,n,npt,r,s,phi)
        ELSE IF (mod(et,2) == 0) THEN
          CALL quad_basis(p,n,npt,r,s,phi)
        ENDIF
        
        DO pt = 1,npt
          DO m = 1,n
            i = (m-1)*npt + pt
            fine%l(m,pt,et) = phi(i)      
          ENDDO
        ENDDO
        

        CALL DGETRS("N",n,npt,fine%V(1,1,et),mnnds,fine%ipiv(1,et),fine%l(1,1,et),mnnds,info)           
      
      ENDDO
      
      RETURN
      END SUBROUTINE l_eval
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
 
 
      SUBROUTINE newton(x,y,eln,r,s,hb)

      USE globals, ONLY: pres,coarse
      USE basis, ONLY: tri_basis,quad_basis

      IMPLICIT NONE
      INTEGER :: it,eln,et,p,n,i,mnnds
      INTEGER :: info
      INTEGER :: maxit
      REAL(pres) :: tol
      REAL(pres) :: x,y
      REAL(pres) :: r(1),s(1),hb
      REAL(pres) :: f,g,error
      REAL(pres) :: dfdr,dfds,dgdr,dgds,jac
      REAL(pres) :: phi(coarse%mnnds),dpdr(coarse%mnnds),dpds(coarse%mnnds)
      REAL(pres) :: l(coarse%mnnds,3)
        
      tol = 1d-10
      maxit = 1000
      info = 0
      
      et = coarse%el_type(eln)
      p = coarse%np(et)  
      n = coarse%nnds(et)
      mnnds = coarse%mnnds
        
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
          l(i,2) = dpdr(i)
          l(i,3) = dpds(i)          
        ENDDO
        
        CALL DGETRS("N",n,3,coarse%V(1,1,et),mnnds,coarse%ipiv(1,et),l,mnnds,info)
!         IF (info /= 0 ) PRINT*, "LAPACK ERROR"      
        
        dfdr = 0d0
        dfds = 0d0
        dgdr = 0d0
        dgds = 0d0
        f = 0d0
        g = 0d0        
        
        DO i = 1,n
          
          dfdr = dfdr + l(i,2)*coarse%elxy(i,eln,1)
          dfds = dfds + l(i,3)*coarse%elxy(i,eln,1)
          dgdr = dgdr + l(i,2)*coarse%elxy(i,eln,2)
          dgds = dgds + l(i,3)*coarse%elxy(i,eln,2)
          
          f = f + l(i,1)*coarse%elxy(i,eln,1)
          g = g + l(i,1)*coarse%elxy(i,eln,2)
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
      
      error = max(abs(f),abs(g))
      IF (it >= maxit) THEN
        PRINT("(A,E22.15)"), "   MAX ITERATIONS EXCEEDED, error = ",error
        PRINT("(2(A,F20.15))"), "   r = ",r(1), "   s = ", s(1)
      ELSE       
        PRINT("(A,I7,A,E22.15)"), "   iterations: ",it, "  error = ",error
        PRINT("(2(A,F20.15))"), "   r = ",r(1), "   s = ", s(1)
      ENDIF
      
      hb = 0d0
      DO i = 1,n
        hb = hb + l(i,1)*coarse%elhb(i,eln)
      ENDDO

      RETURN
      END SUBROUTINE newton
!       

      
      
      END MODULE evaluate