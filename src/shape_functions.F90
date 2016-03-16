      MODULE shape_functions
      
      CONTAINS
      
      
      SUBROUTINE shape_functions_eval(nv,p,nnd,npts,r,s,psi,dpdr,dpds)
      
      USE globals, ONLY: rp
      USE basis, ONLY: element_basis
      USE lapack_interfaces
      USE vandermonde, ONLY: vandermonde_matrix
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: nv,p,npts
      INTEGER, INTENT(OUT) :: nnd
      REAL(rp), DIMENSION(:), INTENT(IN) :: r,s
      REAL(rp), DIMENSION(:,:), INTENT(OUT) :: psi
      REAL(rp), DIMENSION(:,:), OPTIONAL, INTENT(OUT) :: dpdr,dpds      
            
      INTEGER :: info
      INTEGER :: ldv,ldp
      INTEGER :: calc_derv
      INTEGER, DIMENSION(:), ALLOCATABLE :: ipiv
      REAL(rp), DIMENSION(:,:), ALLOCATABLE :: V
      
      IF (PRESENT(dpdr) .AND. PRESENT(dpds)) THEN
        calc_derv = 1
      ELSE
        calc_derv = 0
      ENDIF
            
      ldv = (p+1)**2

      ALLOCATE(V(ldv,ldv))
      ALLOCATE(ipiv(ldv))

      ldp = SIZE(psi,1)

      CALL vandermonde_matrix(nv,p,nnd,V)
      
      IF (calc_derv == 0) THEN
      
        CALL element_basis(nv,p,nnd,npts,r,s,psi) 
        CALL DGESV(nnd,npts,V,ldv,ipiv,psi,ldp,info)      
        
      ELSE 
      
        CALL element_basis(nv,p,nnd,npts,r,s,psi,dpdr,dpds) 
        CALL DGESV(nnd,npts,V,ldv,ipiv,psi,ldp,info)    
        CALL DGETRS("N",nnd,npts,V,ldv,ipiv,dpdr,ldp,info) 
        CALL DGETRS("N",nnd,npts,V,ldv,ipiv,dpds,ldp,info)         
      
      ENDIF
      
      END SUBROUTINE shape_functions_eval
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
      
      SUBROUTINE shape_functions_qpts()

      USE globals, ONLY: rp,nel_type,order,nverts,nnds,mnnds, &
                         mnqpta,nqpta,mnqpte,nqpte,np, &
                         qpta,qpte,psia,dpsidr,dpsids
      USE basis, ONLY: element_basis
      USE messenger2, ONLY: myrank
      
      IMPLICIT NONE      
      
      INTEGER :: pt,i,j,dof
      INTEGER :: typ,et,eo,tpts
      INTEGER :: nv,nnd,nqa,nqe,p
      INTEGER :: info      
      REAL(rp) :: r(mnqpta+4*mnqpte),s(mnqpta+4*mnqpte)

      psia = 0d0
      dpsidr = 0d0
      dpsids = 0d0
      
      
      DO typ = 1,2*nel_type
      
        IF (typ <= 4) THEN
          et = typ
        ELSE
          et = typ - 4
        ENDIF
        
        eo = order(typ)
        
        nv = nverts(et)
        nqa = nqpta(et)
        nqe = nv*nqpte(et)
        
        nnd = nnds(eo)
        p = np(eo)      
      
        tpts = nqa+nqe     
      
        DO pt = 1,nqa
          r(pt) = qpta(pt,1,et)
          s(pt) = qpta(pt,2,et)
        ENDDO  
      
        DO i = 1,nqe
          pt = nqa+i
          r(pt) = qpte(i,1,et)
          s(pt) = qpte(i,2,et)
        ENDDO


        CALL shape_functions_eval(nv,p,nnd,tpts,r,s,psia(:,:,typ),dpsidr(:,:,typ),dpsids(:,:,typ))                 
!       
!       PRINT*, "psi : "      
!       DO i = 1,nnd
!         PRINT("(300(F27.17))"), (psia(i,j,et), j = 1,tpts)
!       ENDDO         

!       PRINT*, "dpsidr : "      
!       DO i = 1,nnd
!         PRINT("(300(F27.17))"), (dpsidr(i,j,et), j = 1,tpts)
!       ENDDO     
! 
!       PRINT*, "dpsids : "      
!       DO i = 1,nnd
!         PRINT("(300(F27.17))"), (dpsids(i,j,et), j = 1,tpts)
!       ENDDO            

       ENDDO
       
      END SUBROUTINE shape_functions_qpts
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
       
      SUBROUTINE shape_functions_vertex()

      USE globals, ONLY: rp,norder,nnds,mnnds,np, &
                         psiv
      USE basis, ONLY: element_nodes
      USE messenger2, ONLY: myrank
      
      IMPLICIT NONE      
      
      INTEGER :: pt,i,j,dof
      INTEGER :: et,eo,npts,p,nnd
      INTEGER :: info      
      REAL(rp) :: r(mnnds),s(mnnds)
       

      ! Evaluates linear shape functions at straight/curved element nodal sets
      !
      ! Used to create additional nodes which can then be adjusted to make 
      ! straight elements curved.  See curvilinear.F90
      
      psiv = 0d0
      
      DO eo = 1,norder

        p = np(eo)           
        
        IF (mod(eo,2) == 1) THEN
          et = 1    
        ELSE IF (mod(eo,2) == 0) THEN
          et = 2
        ENDIF  

        CALL element_nodes(eo,1,p,npts,r,s)
        CALL shape_functions_eval(eo,np(et),nnd,npts,r,s,psiv(:,:,eo))   
              
!         IF (myrank == 0) THEN
!           PRINT*, "psiv"
!           DO i = 1,nnd
!             PRINT "(300(f27.17))", (psiv(i,j,eo), j = 1,npts)
!           ENDDO
!           PRINT*," "
!         ENDIF

      ENDDO
      
      END SUBROUTINE shape_functions_vertex
!       
!       
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!       
!       
      SUBROUTINE shape_functions_curve()

      USE globals, ONLY: rp,norder,nnds,mnnds,np, &
                         psic
      USE basis, ONLY: element_nodes
      USE messenger2, ONLY: myrank
      
      IMPLICIT NONE      
      
      INTEGER :: pt,i,j,dof
      INTEGER :: et,eo,npts,p,nnd
      INTEGER :: info      
      REAL(rp) :: r(mnnds),s(mnnds)
      
      ! Evaluates curvilinear shape functions at high-order batymetry nodal sets
      !
      ! Used to compute function-specified bathymetry at high-order batymetry nodes
      ! for curved elements. See curvilinear.F90      
       
      psic = 0d0 
       
      DO eo = 1,norder

        p = np(eo)       

        IF (mod(eo,2) == 1) THEN
          et = 3    
        ELSE IF (mod(eo,2) == 0) THEN
          et = 4
        ENDIF     

        CALL element_nodes(eo,1,p,npts,r,s)
        CALL shape_functions_eval(eo,np(et),nnd,npts,r,s,psic(:,:,eo))             
                
!         IF (myrank == 0) THEN
!           PRINT*, "psic"        
!           DO i = 1,nnd
!             PRINT "(20(f27.17))", (psic(i,j,eo), j = 1,npts)
!           ENDDO
!           PRINT*," "
!         ENDIF

      ENDDO     
      
      END SUBROUTINE shape_functions_curve
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
       
       
      SUBROUTINE shape_functions_edge()

      USE globals, ONLY: rp,nel_type, &
                         mnqpte,nqpte,np,mnp, &
                         qpte,psie,dpsidxi, &
                         Ve,ipive
      USE basis, ONLY: jacobi,djacobi
      USE messenger2, ONLY: myrank
      
      IMPLICIT NONE      
      
      INTEGER :: pt,i
      INTEGER :: et
      INTEGER :: nqe,p,n
      INTEGER :: info      
      REAL(rp) :: xi(mnqpte)
      REAL(rp) :: phi(mnqpte),dphi(mnqpte)
      
      ! Evaulates 1-D edge shape functions at 1-D guass points
      !
      ! Used in computing the Jacobians for the edge integrals.
      ! See edge_transformation.F90
       
      DO et = 1,nel_type
       
        p = np(et)    
        n = p+1
        nqe = nqpte(et)        
       
        DO pt = 1,nqe
          xi(pt) = qpte(pt,2,et) ! These are the 1-D guass points
        ENDDO
      
        DO i = 0,p
      
          CALL jacobi(0,0,i,xi,nqe,phi)
          CALL djacobi(0,0,i,xi,nqe,dphi)
        
          DO pt = 1,nqe          
            psie(i+1,pt,et) = phi(pt)
            dpsidxi(i+1,pt,et) = dphi(pt)
          ENDDO
      
        ENDDO

        CALL DGETRS("N",n,nqe,Ve(1,1,et),mnp,ipive(1,et),psie(1,1,et),mnp,info)
        CALL DGETRS("N",n,nqe,Ve(1,1,et),mnp,ipive(1,et),dpsidxi(1,1,et),mnp,info)
 
      ENDDO

      RETURN
      END SUBROUTINE shape_functions_edge
      
      END MODULE shape_functions