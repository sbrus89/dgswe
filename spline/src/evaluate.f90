      MODULE evaluate
      
      USE globals, ONLY: rp
      USE lapack_interfaces
      
      IMPLICIT NONE
        
      
      CONTAINS
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!             
      
      SUBROUTINE vandermonde()
      
      USE globals, ONLY: ctp,np,nnds,mnnds,nel_type,V,ipiv
      USE basis, ONLY: tri_nodes,tri_basis,quad_nodes,quad_basis
      
      IMPLICIT NONE
      INTEGER :: et,pt,dof,i,n
      REAL(rp) :: r(mnnds),s(mnnds)
      INTEGER :: info  
      
      PRINT "(A)", "Computing Vandermode matrix..."
      
      ALLOCATE(V(mnnds,mnnds,nel_type))
      ALLOCATE(ipiv(mnnds,nel_type))
      
      DO et = 1,nel_type

        IF (mod(et,2) == 1) THEN
          CALL tri_nodes(1,np(et),n,r,s)
          CALL tri_basis(np(et),n,n,r,s,V(:,:,et))       
        ELSE IF (mod(et,2) == 0) THEN
          CALL quad_nodes(1,np(et),n,r,s)
          CALL quad_basis(np(et),n,n,r,s,V(:,:,et))
        ENDIF
        
        
        CALL DGETRF(n,n,V(1,1,et),mnnds,ipiv(1,et),info)        
!         DO pt = 1,n
!             PRINT("(100(e15.5))"), (V(dof,pt,et), dof = 1,n)
!         ENDDO        
!         PRINT*, " "

        
      ENDDO
      
      
      
      RETURN      
      END SUBROUTINE vandermonde
       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       

      SUBROUTINE transformation()
      
      USE globals, ONLY: mnnds,nel_type,nnds,np,l,dldr,dlds,V,ipiv
      USE basis, ONLY: tri_nodes,quad_nodes,tri_basis,quad_basis
      
      IMPLICIT NONE
      
      INTEGER :: et,pt,i,j
      INTEGER :: n,p,np1,pn,nnd
      INTEGER :: info
      REAL(rp), DIMENSION(mnnds+1) :: r,s
      
      PRINT "(A)", "Computing interpolating polynomials..."        
      
      ALLOCATE(l(mnnds,mnnds+1,nel_type))
      ALLOCATE(dldr(mnnds,mnnds+1,nel_type))
      ALLOCATE(dlds(mnnds,mnnds+1,nel_type))
      
      DO et = 1,nel_type
        n = nnds(et)
        p = np(et)
      
        IF (mod(et,2) == 1) THEN    
          pn = np(3)
          CALL tri_nodes(1,pn,nnd,r,s)
        ELSE IF (mod(et,2) == 0) THEN
          pn = np(4)
          CALL quad_nodes(1,pn,nnd,r,s)
        ENDIF     
        
        np1 = nnd+1
        
        IF (mod(et,2) == 1) THEN
          r(np1) = -1d0/3d0
          s(np1) = -1d0/3d0
        ELSE IF (mod(et,2) == 0) THEN
          r(np1) = 0d0
          s(np1) = 0d0
        ENDIF        
        
        IF (mod(et,2) == 1) THEN
          CALL tri_basis(p,n,np1,r,s,l(:,:,et),dldr(:,:,et),dlds(:,:,et))
        ELSE IF (mod(et,2) == 0) THEN
          CALL quad_basis(p,n,np1,r,s,l(:,:,et),dldr(:,:,et),dlds(:,:,et))
        ENDIF     
        
!         DO i = 1,np1
!           PRINT*, r(i,et),s(i,et)
!         ENDDO
!         PRINT*, ""        

        CALL DGETRS("N",n,np1,V(1,1,et),mnnds,ipiv(1,et),l(1,1,et),mnnds,info)
        CALL DGETRS("N",n,np1,V(1,1,et),mnnds,ipiv(1,et),dldr(1,1,et),mnnds,info)
        CALL DGETRS("N",n,np1,V(1,1,et),mnnds,ipiv(1,et),dlds(1,1,et),mnnds,info)
        
      ENDDO
      
      
      RETURN
      END SUBROUTINE transformation      


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE grid_size(mesh)
      
      USE globals, ONLY: grid
      
      IMPLICIT NONE
      
      INTEGER :: el
      REAL(rp) :: x1,x2,x3
      REAL(rp) :: y1,y2,y3
      REAL(rp) :: a,b,c,s,r
      
      TYPE(grid) :: mesh
      
      ALLOCATE(mesh%h(mesh%ne))
                 
      DO el = 1,mesh%ne
        x1 = mesh%xy(1,mesh%ect(1,el))
        x2 = mesh%xy(1,mesh%ect(2,el))
        x3 = mesh%xy(1,mesh%ect(3,el))
        
        y1 = mesh%xy(2,mesh%ect(1,el))
        y2 = mesh%xy(2,mesh%ect(2,el))
        y3 = mesh%xy(2,mesh%ect(3,el))
        
        a = sqrt((x1-x2)**2+(y1-y2)**2)
        b = sqrt((x2-x3)**2+(y2-y3)**2)
        c = sqrt((x3-x1)**2+(y3-y1)**2)
        
        s = .5d0*(a+b+c)
        r = sqrt((s-a)*(s-b)*(s-c)/s)
        
        mesh%h(el) = 2d0*r
      ENDDO
      
      RETURN
      END SUBROUTINE grid_size
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE invcpp(n,np,mnp,xyh,xyh2)
      
      USE globals, ONLY: Erad,phi0,lambda0
      
      IMPLICIT NONE
      
      INTEGER :: i,pt
      INTEGER :: n,np,mnp
      REAL(rp) :: xyh(mnp,n,3)
      REAL(rp), OPTIONAL :: xyh2(mnp,n,3)      
      
      IF(PRESENT(xyh2)) THEN
      
        DO i = 1,n
          DO pt = 1,np
            xyh2(pt,i,1) = xyh(pt,i,1)/(Erad*cos(phi0))+lambda0
            xyh2(pt,i,2) = xyh(pt,i,2)/Erad
            xyh2(pt,i,3) = xyh(pt,i,3)           
          ENDDO
        ENDDO      
      
      ELSE
      
        DO i = 1,n
          DO pt = 1,np
            xyh(pt,i,1) = xyh(pt,i,1)/(Erad*cos(phi0))+lambda0
            xyh(pt,i,2) = xyh(pt,i,2)/Erad
            xyh(pt,i,3) = xyh(pt,i,3)           
          ENDDO
        ENDDO
        
      ENDIF
      
      RETURN
      END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      
      END MODULE evaluate