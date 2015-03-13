      SUBROUTINE shape_functions()

      USE globals, ONLY: pres,nel_type,nverts,nnds,mnnds,mnqpta,nqpta,mnqpte,nqpte,np, &
                         qpta,qpte,psia,dpsidr,dpsids, &
                         Vand,ipiv
      USE basis, ONLY: tri_basis,quad_basis
      
      IMPLICIT NONE      
      
      INTEGER :: pt,i,j,dof
      INTEGER :: et,tpts
      INTEGER :: nv,nnd,nqa,nqe,p
      INTEGER :: info      
      REAL(pres) :: r(mnqpta+4*mnqpte),s(mnqpta+4*mnqpte)
      REAL(pres) :: phi(mnnds*(mnqpta+4*mnqpte)),dphidr(mnnds*(mnqpta+4*mnqpte)),dphids(mnnds*(mnqpta+4*mnqpte))

      psia = 0d0
      dpsidr = 0d0
      dpsids = 0d0
      
      DO et = 1,nel_type
      
        nv = nverts(et)
        nnd = nnds(et)
        nqa = nqpta(et)
        nqe = nv*nqpte(et)
        p = np(et)      
      
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


        IF(mod(et,2) == 1) THEN
          CALL tri_basis(p,nnd,tpts,r,s,phi,dphidr,dphids)             
        ELSE IF (mod(et,2) == 0) THEN
          CALL quad_basis(p,nnd,tpts,r,s,phi,dphidr,dphids)  
        ENDIF
      
        DO pt = 1,tpts
          DO dof = 1,nnd
            i = (dof-1)*tpts + pt  
          
            psia(dof,pt,et) = phi(i)
            dpsidr(dof,pt,et) = dphidr(i)
            dpsids(dof,pt,et) = dphids(i)
          ENDDO
        ENDDO      
      

!       PRINT*, "RHS matrix: "      
!       DO i = 1,nnds
!         PRINT("(20(F15.5))"), (l(i,j), j = 1,tpts)
!       ENDDO      
      
        CALL DGETRS("N",nnd,tpts,Vand(1,1,et),mnnds,ipiv(1,et),psia(1,1,et),mnnds,info)
        CALL DGETRS("N",nnd,tpts,Vand(1,1,et),mnnds,ipiv(1,et),dpsidr(1,1,et),mnnds,info)      
        CALL DGETRS("N",nnd,tpts,Vand(1,1,et),mnnds,ipiv(1,et),dpsids(1,1,et),mnnds,info)  
      
!       PRINT*, "Psi : "      
!       DO i = 1,nnd
!         PRINT("(300(F27.17))"), (psia(i,j,et), j = 1,tpts)
!       ENDDO         
      
 
      ENDDO

      RETURN
      END SUBROUTINE shape_functions