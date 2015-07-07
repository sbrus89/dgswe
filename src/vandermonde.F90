      SUBROUTINE vandermonde()
      
      USE globals, ONLY: pres,np,mnp,nnds,mnnds,norder,nel_type,Va,ipiva,Ve,ipive
      USE basis, ONLY: tri_nodes,tri_basis,quad_nodes,quad_basis,lglpts,jacobi
      USE read_dginp, ONLY: ctp
      
      IMPLICIT NONE
      INTEGER :: et,pt,dof,i,n,p
      REAL(pres) :: r(mnnds),s(mnnds),xi(mnp)
      REAL(pres) :: phi(mnnds*mnnds)
      INTEGER :: info        
      
      
      ALLOCATE(Va(mnnds,mnnds,norder))
      ALLOCATE(ipiva(mnnds,norder))
      
      ALLOCATE(Ve(mnp,mnp,norder))
      ALLOCATE(ipive(mnp,norder))
      
      Va = 0d0
      ipiva = 0
      Ve = 0d0
      ipive = 0
      
      DO et = 1,norder
        n = nnds(et)
        p = np(et)
        IF (mod(et,2) == 1) THEN
          CALL tri_nodes(1,p,n,r,s)
          CALL tri_basis(p,n,n,r,s,phi)       
        ELSE IF (mod(et,2) == 0) THEN
          CALL quad_nodes(1,p,n,r,s)
          CALL quad_basis(p,n,n,r,s,phi)
        ENDIF
        
        DO pt = 1,n
          DO dof = 1,n
            i = (dof-1)*n + pt
            Va(dof,pt,et) = phi(i)
          ENDDO
        ENDDO
                
        CALL DGETRF(n,n,Va(1,1,et),mnnds,ipiva(1,et),info)  
        
!         DO dof = 1,n
!             PRINT("(100(e15.5))"), (Va(dof,pt,et), pt = 1,n)
!         ENDDO        
!         PRINT*, " "    

      ENDDO
      
      
      DO et = 1,nel_type
      
        p = np(et)        
        n = p+1
      
        CALL lglpts(p,xi)       
      
        DO i = 0,p      
      
          CALL jacobi(0,0,i,xi,n,phi)
        
          DO pt = 1,n
            Ve(i+1,pt,et) = phi(pt)
          ENDDO       
        
        ENDDO     
      
        CALL DGETRF(n,n,Ve(1,1,et),mnp,ipive(1,et),info)       
        
      ENDDO
                  
      RETURN      
      END SUBROUTINE vandermonde