      SUBROUTINE bathymetry_interp()

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
      REAL(pres) :: xpt,ypt,detJ,Sp
      REAL(pres) :: dxdr,dxds,dydr,dyds
      REAL(pres) :: drdx,drdy,dsdx,dsdy


      OPEN(unit=45, file='dhb.d')
      OPEN(unit=46, file='hb.d')
      WRITE(45,*) ne,16
      WRITE(46,*) ne, 16
    
      
      
      DO el = 1,ne
        et = el_type(el)     
        
     pts: DO pt = 1,nqpta(et)        
            dxdr = 0d0
            dxds = 0d0
            dydr = 0d0
            dyds = 0d0
            
            xpt = 0d0
            ypt = 0d0
          
            DO nd = 1,nnds(et)
              x = elxy(nd,el,1)
              y = elxy(nd,el,2)
          
              dxdr = dxdr + dpsidr(nd,pt,et)*x
              dxds = dxds + dpsids(nd,pt,et)*x
              dydr = dydr + dpsidr(nd,pt,et)*y
              dyds = dyds + dpsids(nd,pt,et)*y
              
              xpt = xpt + psia(nd,pt,et)*x
              ypt = ypt + psia(nd,pt,et)*y
                        
            ENDDO
            
            detJ = dxdr*dyds - dxds*dydr
            
            drdx =  dyds/detJ
            drdy = -dxds/detJ
            dsdx = -dydr/detJ
            dsdy =  dxdr/detJ 
            
            IF (coord_sys == 1) THEN
              Sp = 1d0
            ELSE
              Sp = cos(sphi0)/cos(ypt/r_earth)        
            ENDIF            
                                             
            dhbdx_init(el,pt) = 0d0
            dhbdy_init(el,pt) = 0d0
            hbqpta_init(el,pt) = 0d0
            DO nd = 1,nnds(et)  
              hb = elhb(nd,el)
              
              hbqpta_init(el,pt) =  hbqpta_init(el,pt) + psia(nd,pt,et)*hb
              
              dhbdx_init(el,pt) = dhbdx_init(el,pt) + (dpsidr(nd,pt,et)*drdx + dpsids(nd,pt,et)*dsdx)*hb*Sp
              dhbdy_init(el,pt) = dhbdy_init(el,pt) + (dpsidr(nd,pt,et)*drdy + dpsids(nd,pt,et)*dsdy)*hb              
            ENDDO
                     
             IF (et == 3) THEN        
               WRITE(45,"(4(e24.17,1x))") xpt,ypt,dhbdx_init(el,pt),dhbdy_init(el,pt)      
               WRITE(46,"(4(e24.17,1x))") xpt,ypt,hbqpta_init(el,pt)
             ENDIF
             
          ENDDO pts          
                  
      ENDDO
      
      CLOSE(45)
      CLOSE(46)      
      
      
      
      DO ed = 1,ned
      
        el1 = ged2el(1,ed)
        led1 = ged2led(1,ed)
        el2 = ged2el(2,ed)
        led2 = ged2led(2,ed)        
        
        et = el_type(el1) 
        nqa = nqpta(et)
        nqe = nqpte(et)        
        
        DO i = 1,nqe
          pt = nqa + (led1-1)*nqe+i
          edpt = (led1-1)*nqe+i          
          
          hbqpte_init(el1,edpt) = 0d0          
          DO nd = 1,nnds(et)                                 
            hbqpte_init(el1,edpt) = hbqpte_init(el1,edpt) + psia(nd,pt,et)*elhb(nd,el1)     
          ENDDO   
        ENDDO          
                            
          
        IF (el2 /= 0) THEN    

          et = el_type(el2) 
          nqa = nqpta(et)
          nqe = nqpte(et)          
        
          DO i = 1,nqe        
            pt = nqa + (led2-1)*nqe+i
            edpt = (led2-1)*nqe+i          
          
            hbqpte_init(el2,edpt) = 0d0          
            DO nd = 1,nnds(et)                                 
              hbqpte_init(el2,edpt) = hbqpte_init(el2,edpt) + psia(nd,pt,et)*elhb(nd,el2)     
            ENDDO           
          ENDDO

        ENDIF
        
      ENDDO         
      

      RETURN
      END SUBROUTINE bathymetry_interp