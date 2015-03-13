      SUBROUTINE vandermonde()
      
      USE globals, ONLY: pres,ctp,np,nnds,mnnds,nel_type,Vand,ipiv
      USE basis, ONLY: tri_nodes,tri_basis,quad_nodes,quad_basis
      
      IMPLICIT NONE
      INTEGER :: et,pt,dof,i,n
      REAL(pres) :: r(mnnds),s(mnnds)
      REAL(pres) :: phi(mnnds*mnnds)
      INTEGER :: info        
      
      ALLOCATE(Vand(mnnds,mnnds,nel_type))
      ALLOCATE(ipiv(mnnds,nel_type))
      
      Vand = 0d0
      ipiv = 0
      
      DO et = 1,nel_type
        n = nnds(et)
        IF (mod(et,2) == 1) THEN
          CALL tri_nodes(1,np(et),n,r,s)
          CALL tri_basis(np(et),n,n,r,s,phi)       
        ELSE IF (mod(et,2) == 0) THEN
          CALL quad_nodes(1,np(et),n,r,s)
          CALL quad_basis(np(et),n,n,r,s,phi)
        ENDIF
        
        DO pt = 1,n
          DO dof = 1,n
            i = (dof-1)*n + pt
            Vand(dof,pt,et) = phi(i)
          ENDDO
        ENDDO
        

        
        CALL DGETRF(n,n,Vand(1,1,et),mnnds,ipiv(1,et),info)  
        
!         DO dof = 1,n
!             PRINT("(100(e15.5))"), (Vand(dof,pt,et), pt = 1,n)
!         ENDDO        
!         PRINT*, " "                
        

        
      ENDDO
                  
      RETURN      
      END SUBROUTINE vandermonde