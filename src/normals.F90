      SUBROUTINE normals()

      USE globals, ONLY: pres,ned,el_type,ged2el,ged2led,elxy, &
                         nverts,nnds,nqpta,nqpte, &
                         psia,dpsidr,dpsids, &
                         nx_pt,ny_pt,Spe,cfac, &
                         coord_sys,r_earth,sphi0

      IMPLICIT NONE
      
      INTEGER :: i,nd,el,pt,ed,led,edpt
      INTEGER :: et,nnd,nqa,nqe,nv    
      REAL(pres) :: x,y,xpt,ypt
      REAL(pres) :: dxdr,dxds,dydr,dyds
      REAL(pres) :: drdx,drdy,dsdx,dsdy,detJ
      REAL(pres) :: Sp
      REAL(pres) :: nx,ny
    

      
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
      END SUBROUTINE normals