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

      SUBROUTINE function_eval()
      
      USE globals, ONLY: pres,nel_type,nqpta,qpta,mnqpta,fine
      USE basis, ONLY: tri_basis,quad_basis,adcirc_basis
      
      IMPLICIT NONE
      
      INTEGER :: p,n,et,m,i,pt,mnnds,mndof,npt
      INTEGER :: info
      REAL(pres) :: r(mnqpta),s(mnqpta)
      REAL(pres), ALLOCATABLE, DIMENSION(:) :: phi,dpdr,dpds
      
      mnnds = fine%mnnds
      mndof = fine%mndof
      
      ALLOCATE(phi(max(mnnds,mndof)*mnqpta))
      ALLOCATE(dpdr(max(mnnds,mndof)*mnqpta))
      ALLOCATE(dpds(max(mnnds,mndof)*mnqpta))      
      
      ! evaluate shape functions
      
      ALLOCATE(fine%l(mnnds,mnqpta,nel_type))   
      ALLOCATE(fine%dldr(mnnds,mnqpta,nel_type))
      ALLOCATE(fine%dlds(mnnds,mnqpta,nel_type))
      
      DO et = 1,nel_type
        p = fine%np(et)     ! transformation order
        n = fine%nnds(et) ! transformation nodes
        npt = nqpta(et)
        
        DO pt = 1,npt
          r(pt) = qpta(pt,1,et) 
          s(pt) = qpta(pt,2,et)
        ENDDO
        
        IF (mod(et,2) == 1) THEN
          CALL tri_basis(p,n,npt,r,s,phi,dpdr,dpds)
        ELSE IF (mod(et,2) == 0) THEN
          CALL quad_basis(p,n,npt,r,s,phi,dpdr,dpds)
        ENDIF
        
        DO pt = 1,npt
          DO m = 1,n
            i = (m-1)*npt + pt
            fine%l(m,pt,et) = phi(i) 
            fine%dldr(m,pt,et) = dpdr(i)
            fine%dlds(m,pt,et) = dpds(i)            
          ENDDO
        ENDDO
        
        CALL DGETRS("N",n,npt,fine%V(1,1,et),mnnds,fine%ipiv(1,et),fine%l(1,1,et),mnnds,info)  
        CALL DGETRS("N",n,npt,fine%V(1,1,et),mnnds,fine%ipiv(1,et),fine%dldr(1,1,et),mnnds,info)      
        CALL DGETRS("N",n,npt,fine%V(1,1,et),mnnds,fine%ipiv(1,et),fine%dlds(1,1,et),mnnds,info)              
      
      ENDDO
      
      
      ! evaluate basis functions
      
      ALLOCATE(fine%phi(mndof,mnqpta,nel_type))  
      
      DO et = 1,nel_type
        p = fine%p        ! solution order
        n = fine%ndof(et) ! solution degrees of freedom
        npt = nqpta(et) 
        
        DO pt = 1,npt
          r(pt) = qpta(pt,1,et) 
          s(pt) = qpta(pt,2,et)
        ENDDO
        
        IF (mod(et,2) == 1) THEN
          CALL tri_basis(p,n,npt,r,s,phi)
!           CALL adcirc_basis(p,n,npt,r,s,phi)
        ELSE IF (mod(et,2) == 0) THEN
          CALL quad_basis(p,n,npt,r,s,phi)
        ENDIF
        
        DO pt = 1,npt
          DO m = 1,n
            i = (m-1)*npt + pt
            fine%phi(m,pt,et) = phi(i)      
          ENDDO
        ENDDO
      
      ENDDO      
      
      RETURN
      END SUBROUTINE function_eval
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      SUBROUTINE detJ_eval()
      
      USE globals, ONLY: pres,fine,mnqpta,nqpta
      
      IMPLICIT NONE
      
      INTEGER :: el,pt,nd
      INTEGER :: et
      REAL(pres) :: x,y
      REAL(pres) :: dxdr,dxds,dydr,dyds
      
      ALLOCATE(fine%detJ(mnqpta,fine%ne))
      
      DO el = 1,fine%ne
        et = fine%el_type(el)
        
        DO pt = 1,nqpta(et)
          dxdr = 0d0
          dxds = 0d0
          dydr = 0d0
          dyds = 0d0
          
          DO nd = 1,fine%nnds(et)
            x = fine%elxy(nd,el,1)
            y = fine%elxy(nd,el,2)
          
            dxdr = dxdr + fine%dldr(nd,pt,et)*x
            dxds = dxds + fine%dlds(nd,pt,et)*x
            dydr = dydr + fine%dldr(nd,pt,et)*y
            dyds = dyds + fine%dlds(nd,pt,et)*y
                       
          ENDDO
            
          fine%detJ(pt,el) = dxdr*dyds - dxds*dydr        
          
        ENDDO
      
      ENDDO
      
      RETURN
      END SUBROUTINE detJ_eval      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
 
 
      SUBROUTINE newton(x,y,npt,eln,r,s,hb)

      USE globals, ONLY: pres,coarse
      USE basis, ONLY: tri_basis,quad_basis

      IMPLICIT NONE
      INTEGER :: eln,npt
      INTEGER :: it,i,m,pt
      INTEGER et,p,n,mnnds     
      INTEGER :: info
      INTEGER :: maxit
      REAL(pres) :: tol,error
      REAL(pres), DIMENSION(npt) :: x,y
      REAL(pres), DIMENSION(npt) :: r,s,hb
      REAL(pres), DIMENSION(npt) :: f,g
      REAL(pres), DIMENSION(npt) :: dfdr,dfds,dgdr,dgds,jac
      REAL(pres), DIMENSION(npt*coarse%mnnds) :: phi,dpdr,dpds
      REAL(pres) :: l(coarse%mnnds,3*npt)
        
      tol = 1d-10
      maxit = 1000
      info = 0
      
      et = coarse%el_type(eln)
      p = coarse%np(et)  
      n = coarse%nnds(et)
      mnnds = coarse%mnnds
        
      IF (mod(et,2) == 1) THEN
        DO pt = 1,npt      
          r(pt) = -1d0/3d0
          s(pt) = -1d0/3d0
        ENDDO
      ELSE IF (mod(et,2) == 0) THEN
        DO pt = 1,npt
          r(pt) = 0d0
          s(pt) = 0d0
        ENDDO
      ENDIF

      DO it = 1,maxit
           
        IF (mod(et,2) == 1) THEN
          CALL tri_basis(p,n,npt,r,s,phi,dpdr,dpds)
        ELSE IF (mod(et,2) == 0) THEN
          CALL quad_basis(p,n,npt,r,s,phi,dpdr,dpds)
        ENDIF
        
        DO pt = 1,npt
          DO m = 1,n  
            i = (m-1)*npt + pt
            
            l(m,pt) = phi(i)
            l(m,npt+pt) = dpdr(i)
            l(m,2*npt+pt) = dpds(i)          
          ENDDO
        ENDDO
        
        CALL DGETRS("N",n,3*npt,coarse%V(1,1,et),mnnds,coarse%ipiv(1,et),l,mnnds,info)
!         IF (info /= 0 ) PRINT*, "LAPACK ERROR"      
        
        dfdr = 0d0
        dfds = 0d0
        dgdr = 0d0
        dgds = 0d0
        f = 0d0
        g = 0d0        
        
        DO pt = 1,npt
          DO i = 1,n
            
            dfdr(pt) = dfdr(pt) + l(i,npt+pt)*coarse%elxy(i,eln,1)
            dfds(pt) = dfds(pt) + l(i,2*npt+pt)*coarse%elxy(i,eln,1)
            dgdr(pt) = dgdr(pt) + l(i,npt+pt)*coarse%elxy(i,eln,2)
            dgds(pt) = dgds(pt) + l(i,2*npt+pt)*coarse%elxy(i,eln,2)
          
            f(pt) = f(pt) + l(i,pt)*coarse%elxy(i,eln,1)
            g(pt) = g(pt) + l(i,pt)*coarse%elxy(i,eln,2)
          ENDDO
        ENDDO
        
        DO pt = 1,npt
          jac(pt) = dfdr(pt)*dgds(pt) - dgdr(pt)*dfds(pt)
        
          f(pt) = f(pt) - x(pt)
          g(pt) = g(pt) - y(pt)
        
          r(pt) = r(pt) - (1d0/jac(pt))*( dgds(pt)*f(pt) - dfds(pt)*g(pt))
          s(pt) = s(pt) - (1d0/jac(pt))*(-dgdr(pt)*f(pt) + dfdr(pt)*g(pt))         
        ENDDO
        
        error = MAX(ABS(MAXVAL(f)) , ABS(MAXVAL(g)))
        IF ( error < tol) THEN
          EXIT
        ENDIF        
               
      ENDDO
      

      IF (it >= maxit) THEN
        PRINT("(A,E22.15)"), "   MAX ITERATIONS EXCEEDED, error = ",error
        PRINT("(2(A,F20.15))"), "   r = ",r(1), "   s = ", s(1)
      ELSE       
!         PRINT("(A,I7,A,E22.15)"), "   iterations: ",it, "  error = ",error
!         PRINT("(2(A,F20.15))"), "   r = ",r(1), "   s = ", s(1)
      ENDIF
      
      hb = 0d0
      DO pt = 1,npt
        DO i = 1,n
          hb(pt) = hb(pt) + l(i,pt)*coarse%elhb(i,eln)
        ENDDO
      ENDDO

      RETURN
      END SUBROUTINE newton
!       

      
      
      END MODULE evaluate