      MODULE shape_functions
      
      CONTAINS
      
      
      SUBROUTINE shape_functions_eval(nv,p,nnd,npts,r,s,V,psi,dpdr,dpds)
      
      USE globals, ONLY, pres
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: nv,p,npts
      INTEGER, INTENT(OUT) :: nnd
      REAL(pres), DIMENSION(:), INTENT(IN) :: r,s
      REAL(pres), DIMENSION(:,:), INTENT(IN) :: V
      REAL(pres), DIMENSION(:,:), INTENT(OUT) :: psi
      REAL(pres), DIMENSION(:,:), INTENT(OUT) :: dpdr,dpds      
            
      INTEGER :: info
      INTEGER :: ld1,ld2
      INTEGER :: calc_derv
      INTEGER, DIMENSION(:), ALLOCATABLE :: ipiv
      
      IF (PRESENT(dpdr) .AND. PRESENT(dpds)) THEN
        calc_derv = 1
      ELSE
        calc_derv = 0
      ENDIF
            
      ld1 = SIZE(V,1)
      ld2 = SIZE(psi,1)
      
      ALLOCATE(ipiv(nnd))
      
      IF (calc_derv == 0) THEN
      
        CALL element_basis(nv,p,nnd,npts,r,s,psi) 
        CALL DGESV(nnd,npts,V,ld1,ipiv,psi,ld2,info)      
        
      ELSE 
      
        CALL element_basis(nv,p,nnd,npts,r,s,psi,dpdr,dpds) 
        CALL DGESV(nnd,npts,V,ld1,ipiv,psi,ld2,info)    
        CALL DGESV(nnd,npts,V,ld1,ipiv,dpdr,ld2,info) 
        CALL DGESV(nnd,npts,V,ld1,ipiv,dpds,ld2,info)         
      
      ENDIF
      
      END SUBROUTINE shape_functions_eval
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
      
      SUBROUTINE shape_functions_qpts()

      USE globals, ONLY: pres,nel_type,order,nverts,nnds,mnnds, &
                         mnqpta,nqpta,mnqpte,nqpte,np, &
                         qpta,qpte,psia,dpsidr,dpsids, &
                         Va,ipiva
      USE basis, ONLY: element_basis
      USE messenger2, ONLY: myrank
      
      IMPLICIT NONE      
      
      INTEGER :: pt,i,j,dof
      INTEGER :: typ,et,eo,tpts
      INTEGER :: nv,nnd,nqa,nqe,p
      INTEGER :: info      
      REAL(pres) :: r(mnqpta+4*mnqpte),s(mnqpta+4*mnqpte)

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


        CALL shape_functions_eval(nv,p,nnd,tpts,r,s,Va(:,:,eo),psi(:,:,typ),dpdr(:,:,typ),dpds(:,:,typ))        
        
      PRINT*, "Psi : "      
      DO i = 1,nnd
        PRINT("(300(F27.17))"), (psia(i,j,et), j = 1,tpts)
      ENDDO            

        CALL element_basis(nv,p,nnd,tpts,r,s,psia(:,:,typ),dpsidr(:,:,typ),dpsids(:,:,typ))             
 
!       PRINT*, "RHS matrix: "      
!       DO i = 1,nnds
!         PRINT("(20(F15.5))"), (psia(i,j,typ), j = 1,tpts)
!       ENDDO      
      
        CALL DGETRS("N",nnd,tpts,Va(1,1,eo),mnnds,ipiva(1,eo),psia(1,1,typ),mnnds,info)
        CALL DGETRS("N",nnd,tpts,Va(1,1,eo),mnnds,ipiva(1,eo),dpsidr(1,1,typ),mnnds,info)      
        CALL DGETRS("N",nnd,tpts,Va(1,1,eo),mnnds,ipiva(1,eo),dpsids(1,1,typ),mnnds,info)  
      
      PRINT*, "Psi : "      
      DO i = 1,nnd
        PRINT("(300(F27.17))"), (psia(i,j,et), j = 1,tpts)
      ENDDO         

       ENDDO
       
      END SUBROUTINE shape_functions_qpts
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
       
      SUBROUTINE shape_functions_vertex()

      USE globals, ONLY: pres,norder,nnds,mnnds,np, &
                         psiv,Va,ipiva
      USE basis, ONLY: tri_nodes,quad_nodes,tri_basis,quad_basis
      USE messenger2, ONLY: myrank
      
      IMPLICIT NONE      
      
      INTEGER :: pt,i,j,dof
      INTEGER :: et,eo,npts,p
      INTEGER :: info      
      REAL(pres) :: r(mnnds),s(mnnds)
       

      ! Evaluates linear shape functions at straight/curved element nodal sets
      !
      ! Used to create additional nodes which can then be adjusted to make 
      ! straight elements curved.  See curvilinear.F90
      
      psiv = 0d0
      
      DO eo = 1,norder
        npts = nnds(eo)
        p = np(eo)
        
        IF (mod(eo,2) == 1) THEN
          et = 1
          CALL tri_nodes(1,p,r,s)
          CALL tri_basis(np(et),npts,r,s,psiv(:,:,eo))       
        ELSE IF (mod(eo,2) == 0) THEN
          et = 2
          CALL quad_nodes(1,p,r,s)
          CALL quad_basis(np(et),npts,r,s,psiv(:,:,eo))
        ENDIF        
        
        CALL DGETRS("N",nnds(et),npts,Va(1,1,et),mnnds,ipiva(1,et),psiv(1,1,eo),mnnds,info)   
        
!         IF (myrank == 0) THEN
!           PRINT*, "psiv"
!           DO i = 1,nnds(et)
!             PRINT "(20(e20.6))", (psiv(i,j,eo), j = 1,npts)
!           ENDDO
!           PRINT*," "
!         ENDIF

      ENDDO
      
      END SUBROUTINE shape_functions_vertex
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
      
      SUBROUTINE shape_functions_curve()

      USE globals, ONLY: pres,norder,nnds,mnnds,np, &
                         psic,Va,ipiva
      USE basis, ONLY: tri_nodes,quad_nodes,tri_basis,quad_basis
      USE messenger2, ONLY: myrank
      
      IMPLICIT NONE      
      
      INTEGER :: pt,i,j,dof
      INTEGER :: et,eo,npts,p
      INTEGER :: info      
      REAL(pres) :: r(mnnds),s(mnnds)
      
      ! Evaluates curvilinear shape functions at high-order batymetry nodal sets
      !
      ! Used to compute function-specified bathymetry at high-order batymetry nodes
      ! for curved elements. See curvilinear.F90      
       
      psic = 0d0 
       
      DO eo = 1,norder
        npts = nnds(eo)
        p = np(eo)
        
        IF (mod(eo,2) == 1) THEN
          et = 3
          CALL tri_nodes(1,p,r,s)
          CALL tri_basis(np(et),npts,r,s,psic(:,:,eo))       
        ELSE IF (mod(eo,2) == 0) THEN
          et = 4
          CALL quad_nodes(1,p,r,s)
          CALL quad_basis(np(et),npts,r,s,psic(:,:,eo))
        ENDIF               
        
        CALL DGETRS("N",nnds(et),npts,Va(1,1,et),mnnds,ipiva(1,et),psic(1,1,eo),mnnds,info)   
        
!         IF (myrank == 0) THEN
!           PRINT*, "psic"        
!           DO i = 1,nnds(et)
!             PRINT "(20(e20.6))", (psic(i,j,eo), j = 1,npts)
!           ENDDO
!           PRINT*," "
!         ENDIF

      ENDDO     
      
      END SUBROUTINE shape_functions_curve
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
       
       
      SUBROUTINE shape_functions_edge()

      USE globals, ONLY: pres,nel_type, &
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
      REAL(pres) :: xi(mnqpte)
      REAL(pres) :: phi(mnqpte),dphi(mnqpte)
      
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