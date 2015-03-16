      SUBROUTINE bathymetry_interp(et)

      USE globals, ONLY: pres,el_type,nelnds,ndof,nnds,nqpta,nqpte,nverts, &
                         ne,ned,elxy,elhb,ged2el,ged2led, &
                         hbqpta_init,dhbdx_init,dhbdy_init,hbqpte_init, &
                         psia,dpsidr,dpsids, &
                         coord_sys,sphi0,r_earth
                         

      IMPLICIT NONE
      
      INTEGER :: i,el,pt,nd,ed
      INTEGER :: et,edpt,el1,el2,led1,led2
      INTEGER :: ndf,nnd,nqa,nqe,nv
      REAL(pres) :: x,y,hb
      REAL(pres) :: xpt,ypt,detJa,Sp
      REAL(pres) :: dxdr,dxds,dydr,dyds
      REAL(pres) :: drdx,drdy,dsdx,dsdy


      OPEN(unit=45, file='dhb.d')
      OPEN(unit=46, file='hb.d')
      IF (et == 3) THEN
        WRITE(45,*) ne,16
        WRITE(46,*) ne, 16
      ENDIF      
      
      ndf = ndof(et)
      nnd = nnds(et)
      nqa = nqpta(et)
      nqe = nqpte(et)
      nv = nverts(et)      
      
      
      DO el = 1,ne
        IF (el_type(el) == et) THEN

        
     pts: DO pt = 1,nqa        
            dxdr = 0d0
            dxds = 0d0
            dydr = 0d0
            dyds = 0d0
            
            xpt = 0d0
            ypt = 0d0
          
            DO nd = 1,nelnds(el)
              x = elxy(nd,el,1)
              y = elxy(nd,el,2)
          
              dxdr = dxdr + dpsidr(nd,pt,et)*x
              dxds = dxds + dpsids(nd,pt,et)*x
              dydr = dydr + dpsidr(nd,pt,et)*y
              dyds = dyds + dpsids(nd,pt,et)*y
              
              xpt = xpt + psia(nd,pt,et)*x
              ypt = ypt + psia(nd,pt,et)*y
                        
            ENDDO
            
            detJa = dxdr*dyds - dxds*dydr
            
            drdx =  dyds/detJa
            drdy = -dxds/detJa
            dsdx = -dydr/detJa
            dsdy =  dxdr/detJa    
            
            IF (coord_sys == 1) THEN
              Sp = 1d0
            ELSE
              Sp = cos(sphi0)/cos(ypt/r_earth)        
            ENDIF            
                                             
            dhbdx_init(el,pt) = 0d0
            dhbdy_init(el,pt) = 0d0
            hbqpta_init(el,pt) = 0d0
            DO nd = 1,nelnds(el)  
              hb = elhb(nd,el)
              
              hbqpta_init(el,pt) =  hbqpta_init(el,pt) + psia(nd,pt,et)*hb
              
              dhbdx_init(el,pt) = dhbdx_init(el,pt) + (dpsidr(nd,pt,et)*drdx + dpsids(nd,pt,et)*dsdx)*hb*Sp
              dhbdy_init(el,pt) = dhbdy_init(el,pt) + (dpsidr(nd,pt,et)*drdy + dpsids(nd,pt,et)*dsdy)*hb              
            ENDDO
                     
                     
             WRITE(45,"(4(e24.17,1x))") xpt,ypt,dhbdx_init(el,pt),dhbdy_init(el,pt)      
             WRITE(46,"(4(e24.17,1x))") xpt,ypt,hbqpta_init(el,pt)
          ENDDO pts          
          

        ENDIF
      ENDDO
      
      CLOSE(45)
      CLOSE(46)      
      
      
      
      DO ed = 1,ned
      
        el1 = ged2el(1,ed)
        led1 = ged2led(1,ed)
        el2 = ged2el(2,ed)
        led2 = ged2led(2,ed)        
        
        IF (el_type(el1) == et ) THEN
        
        DO i = 1,nqe
          pt = nqa + (led1-1)*nqe+i
          edpt = (led1-1)*nqe+i          
          
          hbqpte_init(el1,edpt) = 0d0          
          DO nd = 1,nelnds(el1)                                 
            hbqpte_init(el1,edpt) = hbqpte_init(el1,edpt) + psia(nd,pt,et)*elhb(nd,el1)     
          ENDDO   
        ENDDO          
          
        ENDIF          
          
        IF (el2 /= 0) THEN    
        IF (el_type(el2) == et) THEN
        DO i = 1,nqe        
          pt = nqa + (led2-1)*nqe+i
          edpt = (led2-1)*nqe+i          
          
          hbqpte_init(el2,edpt) = 0d0          
          DO nd = 1,nelnds(el2)                                 
            hbqpte_init(el2,edpt) = hbqpte_init(el2,edpt) + psia(nd,pt,et)*elhb(nd,el2)     
          ENDDO           
        ENDDO
        ENDIF
        ENDIF
        
      ENDDO          

      RETURN
      END SUBROUTINE bathymetry_interp