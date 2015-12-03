      MODULE evaluate
      
      USE lapack_interfaces
      
      IMPLICIT NONE
      
      
      
      CONTAINS
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
            
      SUBROUTINE vandermonde(mesh)
      
      USE globals, ONLY: rp,nel_type,grid
      USE basis, ONLY: tri_nodes,tri_basis,quad_nodes,quad_basis
      
      IMPLICIT NONE
      TYPE(grid) :: mesh
      
      INTEGER :: et,n
      REAL(rp) :: r(mesh%mnnds),s(mesh%mnnds)
      INTEGER :: info    
      

      
      ALLOCATE(mesh%V(mesh%mnnds,mesh%mnnds,nel_type+2))
      ALLOCATE(mesh%ipiv(mesh%mnnds,nel_type+2))
      
      ! Evaluate basis functions at reference element nodes
      DO et = 1,nel_type+2
        n = mesh%nnds(et)

        IF (mod(et,2) == 1) THEN
          CALL tri_nodes(1,mesh%np(et),n,r,s)
          CALL tri_basis(mesh%np(et),n,n,r,s,mesh%V(:,:,et))       
        ELSE 
          CALL quad_nodes(1,mesh%np(et),n,r,s)
          CALL quad_basis(mesh%np(et),n,n,r,s,mesh%V(:,:,et))
        ENDIF
        
        
        ! Do LU decomposition 
        CALL DGETRF(n,n,mesh%V(1,1,et),mesh%mnnds,mesh%ipiv(1,et),info)
        
!         DO pt = 1,n
!             PRINT("(100(e15.5))"), (mesh%V(dof,pt,et), dof = 1,n)
!         ENDDO        
!         PRINT*, " "       
        
      ENDDO
      
      
      
      RETURN      
      END SUBROUTINE vandermonde
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
      
      SUBROUTINE re_vert(mesh)
      
      USE globals, ONLY: rp,grid,nel_type
      USE basis, ONLY: tri_nodes, quad_nodes
      
      IMPLICIT NONE
      
      INTEGER :: et,pt
      INTEGER :: n,p
      REAL(rp) :: r(4),s(4)
      TYPE(grid) :: mesh
      
      DO et = 1,nel_type
        IF (mod(et,2) == 1) THEN
          n = mesh%nnds(1)
          p = mesh%np(1)
          CALL tri_nodes(1,p,n,r,s)
        ELSE IF (mod(et,2) == 0) THEN
          n = mesh%nnds(2)
          p = mesh%np(2)
          CALL quad_nodes(1,p,n,r,s)
        ENDIF
        
        DO pt = 1,n
          mesh%rsre(1,pt,et) = r(pt)
          mesh%rsre(2,pt,et) = s(pt)
        ENDDO
        
!         PRINT("(4(F15.5))"), (mesh%rsre(1,pt,et), pt = 1,n)
!         PRINT("(4(F15.5))"), (mesh%rsre(2,pt,et), pt = 1,n)        
!         PRINT*, " "
        
      ENDDO      
      
      RETURN
      END SUBROUTINE re_vert

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     

      SUBROUTINE function_eval()
      
      USE globals, ONLY: rp,nel_type,mnept,nept,ept,eval,base
      USE basis, ONLY: tri_basis,quad_basis,adcirc_basis
      
      IMPLICIT NONE
           
      INTEGER :: p,n,et,m,i,pt,mnnds,mndof,npt
      INTEGER :: info
      
      REAL(rp) :: r(mnept),s(mnept)
      
      mnnds = eval%mnnds
          
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! evaluate shape functions at evaluation points (to compute r,s -> x,y transformtion)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      ALLOCATE(eval%l(mnnds,mnept,nel_type))   
      ALLOCATE(eval%dldr(mnnds,mnept,nel_type))
      ALLOCATE(eval%dlds(mnnds,mnept,nel_type))
      
      DO et = 1,nel_type
        p = eval%np(et)     ! transformation order
        n = eval%nnds(et)   ! transformation nodes
        npt = nept(et)
        
        DO pt = 1,npt
          r(pt) = ept(pt,1,et) 
          s(pt) = ept(pt,2,et)
        ENDDO
        
        IF (mod(et,2) == 1) THEN
          CALL tri_basis(p,n,npt,r,s,eval%l(:,:,et),eval%dldr(:,:,et),eval%dlds(:,:,et))
        ELSE IF (mod(et,2) == 0) THEN
          CALL quad_basis(p,n,npt,r,s,eval%l(:,:,et),eval%dldr(:,:,et),eval%dlds(:,:,et))
        ENDIF       
        
        CALL DGETRS("N",n,npt,eval%V(1,1,et),mnnds,eval%ipiv(1,et),eval%l(1,1,et),mnnds,info)  
        CALL DGETRS("N",n,npt,eval%V(1,1,et),mnnds,eval%ipiv(1,et),eval%dldr(1,1,et),mnnds,info)      
        CALL DGETRS("N",n,npt,eval%V(1,1,et),mnnds,eval%ipiv(1,et),eval%dlds(1,1,et),mnnds,info)              
      
      ENDDO
      
      


      
      RETURN
      END SUBROUTINE function_eval
      
        

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
 
 
      SUBROUTINE newton(mesh,x,y,npt,eln,r,s)

      USE globals, ONLY: rp,grid
      USE basis, ONLY: tri_basis,quad_basis

      IMPLICIT NONE
      INTEGER :: eln,npt
      INTEGER :: it,i,m,pt
      INTEGER et,p,n,mnnds     
      INTEGER :: info
      INTEGER :: maxit
      TYPE(grid), INTENT(IN) :: mesh
      REAL(rp) :: tol,error
      REAL(rp), DIMENSION(npt), INTENT(IN) :: x,y
      REAL(rp), DIMENSION(npt), INTENT(OUT) :: r,s
      REAL(rp), DIMENSION(npt) :: f,g
      REAL(rp), DIMENSION(npt) :: dfdr,dfds,dgdr,dgds,jac
      REAL(rp), DIMENSION(npt*mesh%mnnds,1) :: phi,dpdr,dpds
      REAL(rp) :: l(mesh%mnnds,3*npt)
        
      tol = 1d-10
      maxit = 1000
      info = 0
      
      et = mesh%el_type(eln)
      p = mesh%np(et)  
      n = mesh%nnds(et)
      mnnds = mesh%mnnds
        
      ! Initial guesses  
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
      
        ! Evaluate basis functions at r,s coordinates   
        IF (mod(et,2) == 1) THEN
          CALL tri_basis(p,n,npt,r,s,phi,dpdr,dpds)
        ELSE IF (mod(et,2) == 0) THEN
          CALL quad_basis(p,n,npt,r,s,phi,dpdr,dpds)
        ENDIF
        
        DO pt = 1,npt
          DO m = 1,n  
            i = (m-1)*npt + pt
            
            l(m,pt) = phi(i,1)
            l(m,npt+pt) = dpdr(i,1)
            l(m,2*npt+pt) = dpds(i,1)          
          ENDDO
        ENDDO
        
        ! meshve linear systems to get nodal shape functions/derivatives (l,dldr,dlds) at r,s coordinates
        ! V l(r,s) = phi(r,s), V dldr(r,s) = dpdr(r,s), V dlds(r,s) = dpds(r,s)
        CALL DGETRS("N",n,3*npt,mesh%V(1,1,et),mnnds,mesh%ipiv(1,et),l,mnnds,info)
!         IF (info /= 0 ) PRINT*, "LAPACK ERROR"      
        
        ! Evaluate transformation function/derivatives at r,s coordinates
        dfdr = 0d0
        dfds = 0d0
        dgdr = 0d0
        dgds = 0d0
        f = 0d0
        g = 0d0        
        
        DO pt = 1,npt
          DO i = 1,n
            
            dfdr(pt) = dfdr(pt) + l(i,npt+pt)*mesh%elxy(i,eln,1)
            dfds(pt) = dfds(pt) + l(i,2*npt+pt)*mesh%elxy(i,eln,1)
            dgdr(pt) = dgdr(pt) + l(i,npt+pt)*mesh%elxy(i,eln,2)
            dgds(pt) = dgds(pt) + l(i,2*npt+pt)*mesh%elxy(i,eln,2)
          
            f(pt) = f(pt) + l(i,pt)*mesh%elxy(i,eln,1)
            g(pt) = g(pt) + l(i,pt)*mesh%elxy(i,eln,2)
          ENDDO
        ENDDO
               
        ! Newton iteration               
        DO pt = 1,npt
          jac(pt) = dfdr(pt)*dgds(pt) - dgdr(pt)*dfds(pt)
        
          f(pt) = f(pt) - x(pt) 
          g(pt) = g(pt) - y(pt)
        
          r(pt) = r(pt) - ( dgds(pt)*f(pt) - dfds(pt)*g(pt))/jac(pt)
          s(pt) = s(pt) - (-dgdr(pt)*f(pt) + dfdr(pt)*g(pt))/jac(pt)         
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
      


      RETURN
      END SUBROUTINE newton


      
      
      END MODULE evaluate