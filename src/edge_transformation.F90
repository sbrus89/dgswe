      SUBROUTINE edge_transformation(et,np,nqpts)
      
      USE globals, ONLY: pres,qpte,psie,dpsidxi,ned,ged2el,ged2led,el_type,ect,xy,elxy,ctp,detJe
      USE basis, ONLY: lglpts,jacobi,djacobi
      
      IMPLICIT NONE
      
      INTEGER :: i,j,k,n,m,pt,nd,ed
      INTEGER :: np,nnds,nnp,et,nqpts
      INTEGER :: ipiv(np+1),info     
      INTEGER :: el_in,el_ex,led_in
      INTEGER :: nvert,ind
      REAL(pres) :: xi(np+1+nqpts)      
      REAL(pres) :: Pi(np+1+nqpts)
      REAL(pres) :: dPi(np+1+nqpts)   
      REAL(pres) :: V(np+1,np+1)
      REAL(pres) :: phi(np+1,nqpts),dphidxi(np+1,nqpts)
      REAL(pres) :: l(np+1,nqpts),dldxi(np+1,nqpts)   
      REAL(pres) :: x,y,dxdxi,dydxi
      
      
      nnds = np+1
      
      CALL lglpts(np,xi) 
      
      DO i = 1,nqpts
        pt = nnds + i
        xi(pt) = qpte(i,et)
      ENDDO
      
      m = 0
      DO i = 0,np
      
        m = m + 1
      
        CALL jacobi(0,0,i,xi,nnds+nqpts,Pi)
        CALL djacobi(0,0,i,xi,nnds+nqpts,dPi)
        
        DO pt = 1,nnds
          V(m,pt) = Pi(pt)
        ENDDO
        
        DO k = 1,nqpts
          pt = nnds + k
          
          phi(m,k) = Pi(pt)
          dphidxi(m,k) = dPi(pt)
        ENDDO
      
      ENDDO
      
      
      
      DO pt = 1,nqpts
        DO m = 1,nnds
          l(m,pt) = phi(m,pt)
          dldxi(m,pt) = dphidxi(m,pt)
        ENDDO
      ENDDO
      
      CALL DGETRF(nnds,nnds,V,nnds,ipiv,info)
      CALL DGETRS("N",nnds,nqpts,V,nnds,ipiv,l,nnds,info)
      CALL DGETRS("N",nnds,nqpts,V,nnds,ipiv,dldxi,nnds,info)
      
      DO pt = 1,nqpts
        DO m = 1,nnds
          psie(m,pt,et) = l(m,pt)
          dpsidxi(m,pt,et) = dldxi(m,pt)
        ENDDO
      ENDDO
      
      
      
      DO ed = 1,ned
      
        el_in = ged2el(1,ed)        
        led_in = ged2led(1,ed)
        
        IF(mod(el_type(el_in),2) == 1) THEN
          nvert = 3
        ELSE IF (mod(el_type(el_in),2) == 0) THEN
          nvert = 4
        ENDIF
        
        IF(el_type(el_in) == et) THEN
          
          DO pt = 1,nqpts
          
            dxdxi = 0d0
            dydxi = 0d0
            
            DO nd = 1,np+1
              IF (led_in == nvert-1 .and. nd == np+1) THEN
                ind = 1
              ELSE 
                ind = mod(led_in,nvert)*np+nd
              ENDIF
              
              x = elxy(ind,el_in,1)
              y = elxy(ind,el_in,2)
              
              dxdxi = dxdxi + dldxi(nd,pt)*x
              dydxi = dydxi + dldxi(nd,pt)*y
              
            ENDDO
            
            detJe(ed,pt) = sqrt(dxdxi*dxdxi + dydxi*dydxi)
          
          ENDDO
          
        ENDIF
      
      ENDDO
          
      
      
      RETURN
      END SUBROUTINE edge_transformation