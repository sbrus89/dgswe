      SUBROUTINE area_transformation(et,np,nnds,nqpta,nqpte)

      USE globals, ONLY: pres,ne,ned,el_type,nelnds,xy,ect,ged2el,ged2led,elxy,elhb, &
                         ndof,wpta,dpdr,dpds,phia,mmi_init, &
                         psia,dpsidr,dpsids, &
                         detJa,dpdx_init,dpdy_init,phia,phia_int_init, &
                         depth,dhbdx_init,dhbdy_init, &
                         nx_pt,ny_pt,Spe,cfac, &
                         coord_sys,r_earth,sphi0

      IMPLICIT NONE
      
      INTEGER :: i,j,nd,el,pt,m,dof,ind,ed,led
      INTEGER :: np,nnds,et,nnp,nqpta,nqpte,ndf
      INTEGER :: ipiv(nnds),ipiv2(ndof(et)),work(ndof(et)),info
      INTEGER :: nvert
      REAL(pres) :: V(nnds,nnds)
      REAL(pres) :: l(nnds,nqpta+4*nqpte),dldr(nnds,nqpta+4*nqpte),dlds(nnds,nqpta+4*nqpte)
      REAL(pres) :: phi(nnds,nqpta+4*nqpte),dphidr(nnds,nqpta+4*nqpte),dphids(nnds,nqpta+4*nqpte) 
      REAL(pres) :: x,y,xpt,ypt
      REAL(pres) :: dxdr,dxds,dydr,dyds
      REAL(pres) :: drdx,drdy,dsdx,dsdy,detJ
      REAL(pres) :: Sp
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
      
!       PRINT*, " "
!       PRINT*, "Vandermonde matrix: "
!       DO i = 1,nnds
!           PRINT("(20(F20.15))"), (V(i,j), j = 1,nnds)
!       ENDDO
!       PRINT*, " "
!       PRINT*, "RHS matrix: "      
!       DO i = 1,nnds
!         PRINT("(20(F15.5))"), (phi(i,j), j = 1,nqpta+nvert*nqpte)
!       ENDDO
      
      
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
          psia(m,pt,et)   = l(m,pt)
          dpsidr(m,pt,et) = dldr(m,pt)
          dpsids(m,pt,et) = dlds(m,pt)
        ENDDO
      ENDDO 
     

      ndf = ndof(et)
      
      DO el = 1,ne
        IF (el_type(el) == et) THEN
        
          mm = 0d0
        
     pts: DO pt = 1,nqpta        
            dxdr = 0d0
            dxds = 0d0
            dydr = 0d0
            dyds = 0d0
            
            xpt = 0d0
            ypt = 0d0
          
            DO nd = 1,nelnds(el)
              x = elxy(nd,el,1)
              y = elxy(nd,el,2)
          
              dxdr = dxdr + dldr(nd,pt)*x
              dxds = dxds + dlds(nd,pt)*x
              dydr = dydr + dldr(nd,pt)*y
              dyds = dyds + dlds(nd,pt)*y
              
              xpt = xpt + l(nd,pt)*x
              ypt = ypt + l(nd,pt)*y
                        
            ENDDO
            
            detJa(el,pt) = dxdr*dyds - dxds*dydr
            
            drdx =  dyds/detJa(el,pt)
            drdy = -dxds/detJa(el,pt)
            dsdx = -dydr/detJa(el,pt)
            dsdy =  dxdr/detJa(el,pt)
            
            IF (coord_sys == 1) THEN
              Sp = 1d0
            ELSE
              Sp = cos(sphi0)/cos(ypt/r_earth)        
            ENDIF
                     
            
            DO dof = 1,ndf
              ind = (dof-1)*nqpta + pt        
              
              dpdx_init(el,ind) = wpta(pt,et)*(dpdr(dof,pt,et)*drdx + dpds(dof,pt,et)*dsdx)*detJa(el,pt)*Sp
              dpdy_init(el,ind) = wpta(pt,et)*(dpdr(dof,pt,et)*drdy + dpds(dof,pt,et)*dsdy)*detJa(el,pt)  
              
              phia_int_init(el,ind) = wpta(pt,et)*phia(dof,pt,et)*detJa(el,pt)
            
            ENDDO
            
            dhbdx_init(el,pt) = 0d0
            dhbdy_init(el,pt) = 0d0
            DO nd = 1,nelnds(el)  ! This assumes there is an equal order representation between the bathymetry and the coordinate transformation
              hb = elhb(nd,el)
              
              dhbdx_init(el,pt) = dhbdx_init(el,pt) + (dldr(nd,pt)*drdx + dlds(nd,pt)*dsdx)*hb*Sp
              dhbdy_init(el,pt) = dhbdy_init(el,pt) + (dldr(nd,pt)*drdy + dlds(nd,pt)*dsdy)*hb              
            ENDDO
            
            DO i = 1,ndf
              DO j = 1,ndf
                mm(j,i) = mm(j,i) + wpta(pt,et)*phia(i,pt,et)*phia(j,pt,et)*detJa(el,pt)
              ENDDO
            ENDDO
                     
          ENDDO pts
          
          CALL DGETRF(ndf,ndf,mm,ndf,ipiv2,info)          
          CALL DGETRI(ndf,mm,ndf,ipiv2,work,ndf,info)
          
          m = 1
          DO i = 1,ndf
            DO j = 1,ndf
              mmi_init(el,m) = mm(i,j)
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
          
          xpt = 0d0
          ypt = 0d0
          
          DO nd = 1,nelnds(el)
            x = elxy(nd,el,1)
            y = elxy(nd,el,2)
          
            dxdr = dxdr + dldr(nd,pt)*x
            dxds = dxds + dlds(nd,pt)*x
            dydr = dydr + dldr(nd,pt)*y
            dyds = dyds + dlds(nd,pt)*y
            
            xpt = xpt + l(nd,pt)*x
            ypt = ypt + l(nd,pt)*y
                        
          ENDDO
          
          IF (coord_sys == 1) THEN
            Sp = 1d0
          ELSE
            Sp = cos(sphi0)/cos(ypt/r_earth)        
          ENDIF          
            
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
          ENDIF
          
          nx_pt(ed,i) = nx/sqrt(nx*nx+ny*ny)
          ny_pt(ed,i) = ny/sqrt(nx*nx+ny*ny)           
          
          Spe(ed,i) = Sp
          cfac(ed,i) = ny_pt(ed,i)**2+(nx_pt(ed,i)*Sp)**2
        
        ENDDO
        
        ENDIF
        
      ENDDO
  
      RETURN
      END SUBROUTINE area_transformation
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
      SUBROUTINE tri_trans(et,np,nnds,nqpta,nqpte,V,phi,dphidr,dphids)
      
      USE globals, ONLY: pres,qpta,qpte
      USE basis, ONLY: tri_nodes,tri_basis     
      
      IMPLICIT NONE

      INTEGER :: i,j,k,pt,m,led,dof
      INTEGER :: np,nnds,et,nqpta,nqpte
      INTEGER :: tnds
      REAL(pres) :: r(nnds+nqpta+3*nqpte),s(nnds+nqpta+3*nqpte)     
      REAL(pres) :: p(nnds*(nnds+nqpta+3*nqpte)),dpdr(nnds*(nnds+nqpta+3*nqpte)),dpds(nnds*(nnds+nqpta+3*nqpte))
      REAL(pres) :: V(nnds,nnds)  
      REAL(pres) :: phi(nnds,nqpta+3*nqpte),dphidr(nnds,nqpta+3*nqpte),dphids(nnds,nqpta+3*nqpte)       
      
      ! Get triangular reference element nodes
      CALL tri_nodes(1,np,nnds,r,s)   
      
      DO i = 1,nqpta
        pt = nnds + i
        r(pt) = qpta(i,1,et)
        s(pt) = qpta(i,2,et)
      ENDDO  
      
      DO i = 1,3*nqpte
        pt = nnds+nqpta+i
        r(pt) = qpte(i,1,et)
        s(pt) = qpte(i,2,et)
      ENDDO

      tnds = nnds+nqpta+3*nqpte
      CALL tri_basis(np,nnds,tnds,r,s,p,dpdr,dpds)       
      
      DO pt = 1,nnds
        DO dof = 1,nnds
          i = (dof-1)*tnds + pt
          V(dof,pt) = p(i)
        ENDDO
      ENDDO
      
      DO pt = 1,nqpta+3*nqpte
        DO dof = 1,nnds
          i = nnds + (dof-1)*tnds + pt        
          phi(dof,pt) = p(i)
          dphidr(dof,pt) = dpdr(i)
          dphids(dof,pt) = dpds(i)
        ENDDO
      ENDDO
      
      
      RETURN
      END SUBROUTINE tri_trans
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
      SUBROUTINE quad_trans(et,np,nnds,nqpta,nqpte,V,phi,dphidr,dphids)
      
      USE globals, ONLY: pres,qpta,qpte
      USE basis, ONLY: quad_nodes,quad_basis     
      
      IMPLICIT NONE
      
      INTEGER :: i,dof,pt
      INTEGER :: np,nnds,et,nqpta,nqpte
      INTEGER:: tnds
      REAL(pres) :: r(nnds+nqpta+4*nqpte),s(nnds+nqpta+4*nqpte)    
      REAL(pres) :: p(nnds*(nnds+nqpta+4*nqpte)),dpdr(nnds*(nnds+nqpta+4*nqpte)),dpds(nnds*(nnds+nqpta+4*nqpte))
      REAL(pres) :: V(nnds,nnds)
      REAL(pres) :: phi(nnds,nqpta+4*nqpte),dphidr(nnds,nqpta+4*nqpte),dphids(nnds,nqpta+4*nqpte)  
         
      CALL quad_nodes(1,np,nnds,r,s)   
              
      DO i = 1,nqpta
        pt = nnds + i        
        r(pt) = qpta(i,1,et)
        s(pt) = qpta(i,2,et)
      ENDDO
      
      DO i = 1,4*nqpte
        pt = nnds+nqpta+i
        r(pt) = qpte(i,1,et)
        s(pt) = qpte(i,2,et)
      ENDDO
      
      tnds = nnds+nqpta+4*nqpte
      CALL quad_basis(np,nnds,tnds,r,s,p,dpdr,dpds)
      
      DO pt = 1,nnds 
        DO dof = 1,nnds
          i = (dof-1)*tnds + pt
          V(dof,pt) = p(i)
        ENDDO
      ENDDO
          
      DO pt = 1,nqpta+4*nqpte        
        DO dof = 1,nnds
          i = nnds + (dof-1)*tnds + pt
          phi(dof,pt) = p(i)
          dphidr(dof,pt) = dpdr(i)
          dphids(dof,pt) = dpds(i)
        ENDDO        
      ENDDO   
      
      
      
      RETURN
      END SUBROUTINE quad_trans
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    

