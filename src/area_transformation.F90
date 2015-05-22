      SUBROUTINE area_transformation()

      USE globals, ONLY: pres,ne,ned,el_type,nelnds,mnnds,mndof,xy,ect,ged2el,ged2led,elxy, &
                         ndof,nverts,nnds,nqpta,nqpte,wpta,dpdr,dpds,phia,mmi_init, &
                         psia,dpsidr,dpsids, &
                         detJa,dpdx_init,dpdy_init,phia,phia_int_init, &
                         nx_pt,ny_pt,Spe,cfac, &
                         coord_sys,r_earth,sphi0

      IMPLICIT NONE
      
      INTEGER :: i,j,nd,el,pt,m,dof,ed,led,edpt
      INTEGER :: et,nnd,nqa,nqe,ndf,nv
      INTEGER :: ipiv(mndof),work(mndof*mndof)
      INTEGER ::info
      REAL(pres) :: mm(mndof,mndof)      
      REAL(pres) :: x,y,xpt,ypt
      REAL(pres) :: dxdr,dxds,dydr,dyds
      REAL(pres) :: drdx,drdy,dsdx,dsdy,detJ
      REAL(pres) :: Sp
      REAL(pres) :: nx,ny

    
          

          
      DO el = 1,ne
        et = el_type(el)
        
        nnd = nnds(et)
        nqa = nqpta(et)
        ndf = ndof(et)
        
        mm = 0d0
        
   pts: DO pt = 1,nqa        
          dxdr = 0d0
          dxds = 0d0
          dydr = 0d0
          dyds = 0d0
            
          xpt = 0d0
          ypt = 0d0
          
          DO nd = 1,nnd
            x = elxy(nd,el,1)
            y = elxy(nd,el,2)
          
            dxdr = dxdr + dpsidr(nd,pt,et)*x
            dxds = dxds + dpsids(nd,pt,et)*x
            dydr = dydr + dpsidr(nd,pt,et)*y
            dyds = dyds + dpsids(nd,pt,et)*y
              
            xpt = xpt + psia(nd,pt,et)*x
            ypt = ypt + psia(nd,pt,et)*y
                      
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
            i = (dof-1)*nqa + pt        
              
            dpdx_init(el,i) = wpta(pt,et)*(dpdr(dof,pt,et)*drdx + dpds(dof,pt,et)*dsdx)*detJa(el,pt)*Sp
            dpdy_init(el,i) = wpta(pt,et)*(dpdr(dof,pt,et)*drdy + dpds(dof,pt,et)*dsdy)*detJa(el,pt)  
              
            phia_int_init(el,i) = wpta(pt,et)*phia(dof,pt,et)*detJa(el,pt)
           
          ENDDO
            
            
          DO i = 1,ndf
            DO j = 1,ndf
              mm(j,i) = mm(j,i) + wpta(pt,et)*phia(i,pt,et)*phia(j,pt,et)*detJa(el,pt)
            ENDDO
          ENDDO

        ENDDO pts
          
          
        CALL DGETRF(ndf,ndf,mm,mndof,ipiv,info)       
        CALL DGETRI(ndf,mm,mndof,ipiv,work,ndf*ndf,info)

        m = 1
        DO i = 1,ndf
          DO j = 1,ndf
            mmi_init(el,m) = mm(i,j)
            m = m + 1
          ENDDO
        ENDDO          
          
      ENDDO
      

      
      nx_pt = 0d0
      ny_pt = 0d0          
      
      DO ed = 1,ned
      
        el = ged2el(1,ed)
        led = ged2led(1,ed)
        
        et = el_type(el) 
        
        nnd = nnds(et)
        nqa = nqpta(et)
        nqe = nqpte(et)
        nv = nverts(et)        
        
        DO i = 1,nqe
          pt = nqa + (led-1)*nqe+i
          edpt = (led-1)*nqe+i
          
          dxdr = 0d0
          dxds = 0d0
          dydr = 0d0
          dyds = 0d0
          
          xpt = 0d0
          ypt = 0d0
          

          DO nd = 1,nnd
            x = elxy(nd,el,1)
            y = elxy(nd,el,2)
          
            dxdr = dxdr + dpsidr(nd,pt,et)*x
            dxds = dxds + dpsids(nd,pt,et)*x
            dydr = dydr + dpsidr(nd,pt,et)*y
            dyds = dyds + dpsids(nd,pt,et)*y
            
            xpt = xpt + psia(nd,pt,et)*x
            ypt = ypt + psia(nd,pt,et)*y
                        
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
          
          IF (nv == 3) THEN
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
              
          ELSE IF(nv == 4) THEN        
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
        
      ENDDO

!        
  
      RETURN
      END SUBROUTINE area_transformation
      