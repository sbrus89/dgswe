     MODULE vandermonde
     
     CONTAINS
          
     
     SUBROUTINE vandermonde_matrix(et,p,ndf,V)
     
     USE globals, ONLY: pres
     USE basis, ONLY: tri_nodes,tri_basis,quad_nodes,quad_basis
     
     IMPLICIT NONE
       
     INTEGER :: et,p
     REAL(pres), DIMENSION(:,:) :: V
     
     INTEGER :: npt,ndf
     INTEGER :: info
     REAL(pres), DIMENSION((p+1)**2) :: r,s
   
   
     CALL element_nodes(et,1,p,npt,r,s)
     CALL element_basis(et,p,ndf,npt,r,s,V)
     
     
     END SUBROUTINE vandermonde_matrix
     
     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
     
     
     SUBROUTINE area_vandermonde()
      
      USE globals, ONLY: pres,np,mnp,nnds,mnnds,norder,Va,ipiva
      USE basis, ONLY: tri_nodes,tri_basis,quad_nodes,quad_basis
      
      IMPLICIT NONE
      INTEGER :: et,pt,dof,i,n,p
      REAL(pres) :: r(mnnds),s(mnnds)
      INTEGER :: info        
      
      
      ALLOCATE(Va(mnnds,mnnds,norder))
      ALLOCATE(ipiva(mnnds,norder))
      

      Va = 0d0
      ipiva = 0

      
      DO et = 1,norder
        p = np(et)

        CALL vandermonde_matrix(et,p,n,Va(:,:,et))
                
        CALL DGETRF(n,n,Va(1,1,et),mnnds,ipiva(1,et),info)  
        
!         DO dof = 1,n
!             PRINT("(100(e15.5))"), (Va(dof,pt,et), pt = 1,n)
!         ENDDO        
!         PRINT*, " "    

      ENDDO
      
      END SUBROUTINE area_vandermonde      
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!               
      
      
      SUBROUTINE edge_vandermonde()
      
      USE globals, ONLY: pres,np,mnp,norder,nel_type,Ve,ipive
      USE basis, ONLY: lglpts,jacobi
      
      IMPLICIT NONE
      INTEGER :: et,pt,i,n,p
      REAL(pres) :: xi(mnp)
      REAL(pres) :: phi(mnp*mnp)
      INTEGER :: info          
      
      ALLOCATE(Ve(mnp,mnp,norder))
      ALLOCATE(ipive(mnp,norder)) 
      
      Ve = 0d0
      ipive = 0      
      
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
      END SUBROUTINE edge_vandermonde
      
      
      
      END MODULE vandermonde