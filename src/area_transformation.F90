      SUBROUTINE area_transformation(et,np,nnds,nqpta,nqpte)

      USE globals, ONLY: pres,ne,ned,el_type,nelnds,xy,ect,ged2el,ged2led, &
                         ndof,wpta,dpdr,dpds,phia,mmi, &
                         psia,dpsidr,dpsids, &
                         detJa,dpdx_init,dpdy_init,phia,phia_int_init, &
                         depth,dhbdx_init,dhbdy_init, &
                         nx_pt,ny_pt

      IMPLICIT NONE
      
      INTEGER :: i,j,nd,el,pt,m,dof,ind,ed,led
      INTEGER :: np,nnds,et,nnp,nqpta,nqpte,ndf
      INTEGER :: ipiv(nnds),ipiv2(ndof(et)),work(ndof(et)),info
      INTEGER :: nvert
      REAL(pres) :: V(nnds,nnds)
      REAL(pres) :: l(nnds,nqpta+4*nqpte),dldr(nnds,nqpta+4*nqpte),dlds(nnds,nqpta+4*nqpte)
      REAL(pres) :: phi(nnds,nqpta+4*nqpte),dphidr(nnds,nqpta+4*nqpte),dphids(nnds,nqpta+4*nqpte) 
      REAL(pres) :: x,y
      REAL(pres) :: dxdr,dxds,dydr,dyds
      REAL(pres) :: drdx,drdy,dsdx,dsdy,detJ
      REAL(pres) :: nx,ny
      REAL(pres) :: hb 
      REAL(pres) :: mm(ndof(et),ndof(et))
      REAL(pres) :: area,x1,x2,x3,y1,y2,y3
      


      IF (mod(et,2) == 1) THEN
         
        CALL tri_trans(et,np,nnds,nqpta,nqpte,V,phi,dphidr,dphids)
        nvert = 3
  
      ELSE IF (mod(et,2) == 0) THEN
      
        CALL quad_trans(et,np,nnds,nqpta,nqpte,V,phi,dphidr,dphids)
        nvert = 4
            
      ENDIF
      
      PRINT*, " "
      PRINT*, "Vandermonde matrix: "
      DO i = 1,nnds
          PRINT("(20(F10.4))"), (V(i,j), j = 1,nnds)
      ENDDO
      
      
      
      
      DO pt = 1,nqpta+nvert*nqpte
        DO m = 1,nnds
          l(m,pt) = phi(m,pt)
          dldr(m,pt) = dphidr(m,pt)
          dlds(m,pt) = dphids(m,pt)
        ENDDO
      ENDDO
      
      CALL DGETRF(nnds,nnds,V,nnds,ipiv,info)
      CALL DGETRS("N",nnds,nqpta+nvert*nqpte,V,nnds,ipiv,l,nnds,info)
      CALL DGETRS("N",nnds,nqpta+nvert*nqpte,V,nnds,ipiv,dldr,nnds,info)      
      CALL DGETRS("N",nnds,nqpta+nvert*nqpte,V,nnds,ipiv,dlds,nnds,info)      
      
 
      DO pt = 1,nqpta+nvert*nqpte
        DO m = 1,nnds
          psia(m,pt,et) = l(m,pt)
          dpsidr(m,pt,et) = dldr(m,pt)
          dpsids(m,pt,et) = dlds(m,pt)
        ENDDO
      ENDDO 
     

      ndf = ndof(et)
      
      DO el = 1,ne
        IF (el_type(el) == et) THEN
        
          mm = 0d0
        
          DO pt = 1,nqpta        
            dxdr = 0d0
            dxds = 0d0
            dydr = 0d0
            dyds = 0d0
          
            DO nd = 1,nelnds(el)
              x = xy(1,ect(nd,el))
              y = xy(2,ect(nd,el))
          
              dxdr = dxdr + dldr(nd,pt)*x
              dxds = dxds + dlds(nd,pt)*x
              dydr = dydr + dldr(nd,pt)*y
              dyds = dyds + dlds(nd,pt)*y
                        
            ENDDO
            
            detJa(el,pt) = dxdr*dyds - dxds*dydr
            
            drdx =  dyds/detJa(el,pt)
            drdy = -dxds/detJa(el,pt)
            dsdx = -dydr/detJa(el,pt)
            dsdy =  dxdr/detJa(el,pt)
            
            DO dof = 1,ndf
              ind = (dof-1)*nqpta + pt        
              
              dpdx_init(el,ind) = wpta(pt,et)*(dpdr(dof,pt,et)*drdx + dpds(dof,pt,et)*dsdx)*detJa(el,pt)
              dpdy_init(el,ind) = wpta(pt,et)*(dpdr(dof,pt,et)*drdy + dpds(dof,pt,et)*dsdy)*detJa(el,pt)  
              
              phia_int_init(el,ind) = wpta(pt,et)*phia(dof,pt,et)*detJa(el,pt)
            
            ENDDO
            
            DO nd = 1,nelnds(el)
              hb = depth(ect(nd,el))
              
              dhbdx_init(el,pt) = dhbdx_init(el,pt) + (dldr(nd,pt)*drdx + dlds(nd,pt)*dsdx)*hb
              dhbdy_init(el,pt) = dhbdy_init(el,pt) + (dldr(nd,pt)*drdy + dlds(nd,pt)*dsdy)*hb              
            ENDDO
            
            DO i = 1,ndf
              DO j = 1,ndf
                mm(j,i) = mm(j,i) + wpta(pt,et)*phia(i,pt,et)*phia(j,pt,et)*detJa(el,pt)
              ENDDO
            ENDDO
                     
          ENDDO
          
          CALL DGETRF(ndf,ndf,mm,ndf,ipiv2,info)          
          CALL DGETRI(ndf,mm,ndf,ipiv2,work,ndf,info)
          
          m = 1
          DO i = 1,ndf
            DO j = 1,ndf
              mmi(el,m) = mm(i,j)
              m = m + 1
            ENDDO
          ENDDO
          
!           print*, ' ' 
!           DO i = 1,ndf
!             print "(I5,16(e23.14))", el, (mm(i,j), j = 1,ndf)
!           ENDDO
!           IF (et == 1) THEN
!             x1 = xy(1,ect(1,el))
!             y1 = xy(2,ect(1,el))
! 
!             x2 = xy(1,ect(2,el))
!             y2 = xy(2,ect(2,el))
! 
!              x3 = xy(1,ect(3,el))
!              y3 = xy(2,ect(3,el))
! 
!             area = .5d0*((x2*y3-x3*y2) + (x3*y1-x1*y3) + (x1*y2-x2*y1))          
!             print*, 1d0/area
!           ENDIF
          

        ENDIF
      ENDDO
      
      

      
      DO ed = 1,ned
      
        el = ged2el(1,ed)
        led = ged2led(1,ed)
        
        IF (el_type(el) == et ) THEN
        
        DO i = 1,nqpte
          pt = nqpta + (led-1)*nqpte+i
          
          dxdr = 0d0
          dxds = 0d0
          dydr = 0d0
          dyds = 0d0
          
          DO nd = 1,nelnds(el)
            x = xy(1,ect(nd,el))
            y = xy(2,ect(nd,el))
          
            dxdr = dxdr + dldr(nd,pt)*x
            dxds = dxds + dlds(nd,pt)*x
            dydr = dydr + dldr(nd,pt)*y
            dyds = dyds + dlds(nd,pt)*y
                        
          ENDDO
            
          detJ = dxdr*dyds - dxds*dydr
            
          drdx =  dyds/detJ
          drdy = -dxds/detJ
          dsdx = -dydr/detJ
          dsdy =  dxdr/detJ
          
          IF (nvert == 3) THEN
              SELECT CASE(led)
                CASE(1)
                  nx = drdx + dsdx
                  ny = drdy + dsdy
                CASE(2)
                  nx = -drdx
                  ny = -drdy                
                CASE(3)
                  nx = -dsdx
                  ny = -dsdy               
              END SELECT  
              nx_pt(ed,i) = nx/sqrt(nx*nx+ny*ny)
              ny_pt(ed,i) = ny/sqrt(nx*nx+ny*ny)                     
          ELSE IF(nvert == 4) THEN        
              SELECT CASE(led)
                CASE(1)
                  nx = drdx
                  ny = drdy           
                CASE(2)
                  nx = dsdx
                  ny = dsdy            
                CASE(3)
                  nx = -drdx
                  ny = -drdy             
                CASE(4)
                  nx = -dsdx
                  ny = -dsdy      
              END SELECT       
              nx_pt(ed,i) = nx/sqrt(nx*nx+ny*ny)
              ny_pt(ed,i) = ny/sqrt(nx*nx+ny*ny)       
          ENDIF
        
        ENDDO
        
        ENDIF
        
      ENDDO
  
      RETURN
      END SUBROUTINE area_transformation
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
      SUBROUTINE tri_trans(et,np,nnds,nqpta,nqpte,V,phi,dphidr,dphids)
      
      USE globals, ONLY: pres,qpta,qpte
      USE basis, ONLY: nodes,jacobi,djacobi      
      
      IMPLICIT NONE

      INTEGER :: i,j,k,pt,m,led
      INTEGER :: np,nnds,et,nqpta,nqpte
      REAL(pres) :: r(nnds+nqpta+3*nqpte),s(nnds+nqpta+3*nqpte)
      REAL(pres) :: a(nnds+nqpta+3*nqpte),b(nnds+nqpta+3*nqpte)
      REAL(pres) :: Pi(nnds+nqpta+3*nqpte),Pj(nnds+nqpta+3*nqpte)
      REAL(pres) :: dPi(nnds+nqpta+3*nqpte),dPj(nnds+nqpta+3*nqpte)      
      REAL(pres) :: V(nnds,nnds)  
      REAL(pres) :: phi(nnds,nqpta+3*nqpte),dphidr(nnds,nqpta+3*nqpte),dphids(nnds,nqpta+3*nqpte)
      REAL(pres) :: dpda,dpdb,dadr,dads,ii
      
      r = 0d0
      s = 0d0 
      V = 0d0      
      
      ! Get triangular reference element nodes
      CALL nodes(1,np,nnds,r,s)   
      
      DO i = 1,nqpta
        pt = nnds + i
        r(pt) = qpta(i,1,et)
        s(pt) = qpta(i,2,et)
      ENDDO  
      
      DO led = 1,3
        DO i = 1,nqpte
          pt = nnds+nqpta+(led-1)*nqpte + i
          SELECT CASE(led)
            CASE(1)
              r(pt) = -qpte(i,et)
              s(pt) =  qpte(i,et)
            CASE(2)
              r(pt) = -1d0
              s(pt) = -qpte(i,et)
            CASE(3)
              r(pt) =  qpte(i,et)
              s(pt) = -1d0
          END SELECT
        ENDDO
      ENDDO

      ! Change quadrature points from r,s (master element) to a,b extended coordinates
      DO pt = 1,nnds+nqpta+nqpte
        IF(s(pt) /= 1d0) THEN
          a(pt) = 2d0*(1d0+r(pt))/(1d0-s(pt))-1d0 
        ELSE 
          a(pt) = -1d0
        ENDIF
        b(pt) = s(pt)      
      ENDDO
        
      m = 0
      DO i = 0,np
        DO j = 0,np-i

          m = m+1
          
          Pi = 0d0
          Pj = 0d0
          
          dPi = 0d0          
          dPj = 0d0          

          CALL jacobi(0    ,0,i,a,nnds+nqpta+3*nqpte,Pi)
          CALL jacobi(2*i+1,0,j,b,nnds+nqpta+3*nqpte,Pj)
          
          
          CALL djacobi(0    ,0,i,a,nnds+nqpta+3*nqpte,dPi)          
          CALL djacobi(2*i+1,0,j,b,nnds+nqpta+3*nqpte,dPj)          
          
          DO pt = 1,nnds 
            V(m,pt) = 2d0*Pi(pt)*Pj(pt)*(1d0-b(pt))**i
          ENDDO
          
          ! Calculate function values
          DO k = 1,nqpta+3*nqpte 
            pt = nnds + k
            
            phi(m,k) = 2d0*Pi(pt)*Pj(pt)*(1d0-b(pt))**i
          ENDDO

          ii = real(i,pres)                    
          
          ! Calculate derivative values
          DO k = 1,nqpta+3*nqpte
            pt = nnds + k
            
            dadr = 2d0/(1d0-s(pt))
            dads = 2d0*(1d0+r(pt))/(1d0-s(pt))**2d0            

            dpda = 2d0*dPi(pt)*Pj(pt)*(1d0-b(pt))**ii
            dpdb = 2d0*Pi(pt)*(dPj(pt)*(1d0-b(pt))**ii - ii*(1d0-b(pt))**(ii-1d0)*Pj(pt))
            
            dphidr(m,k) = dpda*dadr
            dphids(m,k) = dpda*dads + dpdb
          ENDDO          

        ENDDO
      ENDDO  
      
      
      RETURN
      END SUBROUTINE tri_trans
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
      SUBROUTINE quad_trans(et,np,nnds,nqpta,nqpte,V,phi,dphidr,dphids)
      
      USE globals, ONLY: pres,qpta,qpte
      USE basis, ONLY: lglpts,jacobi,djacobi      
      
      IMPLICIT NONE
      
      INTEGER :: i,j,k,n,m,pt,led
      INTEGER :: np,nnds,nnp,et,nqpta,nqpte
      REAL(pres) :: xi(np+1)      
      REAL(pres) :: r(nnds+nqpta+4*nqpte),s(nnds+nqpta+4*nqpte)
      REAL(pres) :: Pi(nnds+nqpta+4*nqpte),Pj(nnds+nqpta+4*nqpte)
      REAL(pres) :: dPi(nnds+nqpta+4*nqpte),dPj(nnds+nqpta+4*nqpte)      
      REAL(pres) :: V(nnds,nnds)
      REAL(pres) :: phi(nnds,nqpta+4*nqpte),dphidr(nnds,nqpta+4*nqpte),dphids(nnds,nqpta+4*nqpte)
      
      r = 0d0
      s = 0d0      
      
      ! Get 1-D LGL points 
      CALL lglpts(np,xi)
        
      ! Do tensor product, ordering the nodes counter-clockwise

      ! Find number of loops around refence quad element, excluding middle points
      IF (np <= 2)THEN
        nnp = 1
      ELSE
        IF(mod(np,2) == 1) THEN
          nnp = np-1
        ELSE IF (mod(np,2) == 0) THEN
          nnp = np-2
        ENDIF
      ENDIF
        
      n = 1
      DO k = 1,nnp ! loop over number of loops
         
        ! Edge 4
        DO i = k,np+1 - (k-1)
          j = k
          r(n) = xi(i)
          s(n) = xi(j)
          
          n = n+1
        ENDDO
        
        ! Edge 1
        DO j = 2 + (k-1),np+1 - (k-1)
          i = np+1 - (k-1)
          r(n) = xi(i)
          s(n) = xi(j)
          
          n = n+1
        ENDDO

        ! Edge 2
        DO i = np - (k-1),1 + (k-1),-1
          j = np+1 - (k-1)
          r(n) = xi(i)
          s(n) = xi(j)
          
          n = n+1      
        ENDDO

        ! Edge 3
        DO j = np - (k-1),2 + (k-1),-1
          i = 1 + (k-1)
          r(n) = xi(i)
          s(n) = xi(j)

          n = n+1
        ENDDO    
          
      ENDDO
        
      ! middle point
      IF (mod(np+1,2) == 1) THEN
        i = np/2 + 1
        r(n) = xi(i)
        s(n) = xi(i)
          
      ENDIF
        
        
      PRINT*, ' '
      PRINT*, 'Quadrilateral Points'
      DO n = 1,nnds
        PRINT("(2(f10.4))"), r(n),s(n)
      ENDDO   
      
      DO i = 1,nqpta
        pt = nnds + i
        
        r(pt) = qpta(i,1,et)
        s(pt) = qpta(i,2,et)
      ENDDO
      
      DO led = 1,4
        DO i = 1,nqpte
          pt = nnds+nqpta+(led-1)*nqpte + i
          SELECT CASE(led)
            CASE(1)
              r(pt) = 1d0
              s(pt) = qpte(i,et)
            CASE(2)
              r(pt) = -qpte(i,et)
              s(pt) =  1d0
            CASE(3)
              r(pt) = -1d0
              s(pt) = -qpte(i,et)
            CASE(4)
              r(pt) =  qpte(i,et)
              s(pt) = -1d0
          END SELECT
        ENDDO
      ENDDO
      
      
      m = 0
      DO i = 0,np
        DO j = 0,np

          m = m+1

          Pi = 0d0
          Pj = 0d0
          
          dPi = 0d0
          dPj = 0d0

          CALL jacobi(0,0,i,r,nnds+nqpta+4*nqpte,Pi)
          CALL jacobi(0,0,j,s,nnds+nqpta+4*nqpte,Pj)
          
          CALL djacobi(0,0,i,r,nnds+nqpta+4*nqpte,dPi)
          CALL djacobi(0,0,j,s,nnds+nqpta+4*nqpte,dPj)          

          DO pt = 1,nnds 
            V(m,pt) = 2d0*Pi(pt)*Pj(pt)
          ENDDO
          
          DO k = 1,nqpta+4*nqpte
            pt = nnds + k
          
            phi(m,k) = 2d0*Pi(pt)*Pj(pt)
            dphidr(m,k) = 2d0*dPi(pt)*Pj(pt)
            dphids(m,k) = 2d0*Pi(pt)*dPj(pt)
          ENDDO

        ENDDO
      ENDDO      
      
      
      
      RETURN
      END SUBROUTINE quad_trans
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    

