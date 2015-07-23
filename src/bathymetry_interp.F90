      SUBROUTINE bathymetry_interp()

      USE globals, ONLY: pres,el_type,ed_type,nelnds,mndof,ndof,mnnds,nnds,nqpta,nqpte,nverts,order, &
                         ne,ned,elxy,elhb,ged2el,ged2led, &
                         hbqpta_init,dhbdx_init,dhbdy_init,hbqpte_init,hbm,hbqpted, &
                         Va,ipiva,psia,dpsidr,dpsids, &
                         coord_sys,sphi0,r_earth,recv_edge
                         
      USE read_dginp, ONLY: out_direc                         
                         

      IMPLICIT NONE
      
      INTEGER :: i,j,el,pt,nd,ed,dof,led
      INTEGER :: typ,et,eo,edpt,el1,el2,led1,led2,gp1,gp2
      INTEGER :: ndf,nnd,nqa,nqe,nv
      INTEGER :: info
      REAL(pres) :: x,y,hb
      REAL(pres) :: xpt,ypt,detJ,Sp
      REAL(pres) :: dxdr,dxds,dydr,dyds
      REAL(pres) :: drdx,drdy,dsdx,dsdy
      REAL(pres) :: hbn(mnnds)


      OPEN(unit=45, file='dhb.d')
      OPEN(unit=46, file='hb.d')
      WRITE(45,*) ne,16
      WRITE(46,*) ne, 16
    
      
      
      DO el = 1,ne
      
        et = el_type(el)
        typ = et + 4
        eo = order(typ)        
        
!         PRINT "(I5,200(f10.4))", et,(elhb(nd,el),   nd = 1,nnds(eo))          
        
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
            DO nd = 1,nnds(eo)  
              hb = elhb(nd,el)
              
              hbqpta_init(el,pt) =  hbqpta_init(el,pt) + psia(nd,pt,typ)*hb
              
              dhbdx_init(el,pt) = dhbdx_init(el,pt) + (dpsidr(nd,pt,typ)*drdx + dpsids(nd,pt,typ)*dsdx)*hb*Sp
              dhbdy_init(el,pt) = dhbdy_init(el,pt) + (dpsidr(nd,pt,typ)*drdy + dpsids(nd,pt,typ)*dsdy)*hb              
            ENDDO
                     
!              IF (et == 1) THEN        
               WRITE(45,"(I7,4(e24.17,1x))") el,xpt,ypt,dhbdx_init(el,pt),dhbdy_init(el,pt)      
               WRITE(46,"(I7,4(e24.17,1x))") el,xpt,ypt,hbqpta_init(el,pt)
!              ENDIF
             
          ENDDO pts          
          
                 
                  
      ENDDO
      
      CLOSE(45)
      CLOSE(46)      
      

      
      

      
!       OPEN(unit=47, file='hbe.d')     
      
      DO ed = 1,ned
      
        el1 = ged2el(1,ed)
        led1 = ged2led(1,ed)          
        el2 = ged2el(2,ed)
        led2 = ged2led(2,ed)           
          
        IF (ed_type(ed) == 0) THEN    
        
          DO i = 1,2
          
            el = ged2el(i,ed)
            led = ged2led(i,ed)          
            
            et = el_type(el)
          
            ! Interior edge => assumed to be straight
            IF (mod(et,2) == 1) THEN
              et = 1
            ELSE IF (mod(et,2) == 0) THEN
              et = 2
            ENDIF
          
            nqa = nqpta(et)
            nqe = nqpte(et)                    
        
            typ = et + 4
            eo = order(typ)
           
            DO j = 1,nqe
              pt = nqa + (led-1)*nqe+j
              edpt = (led-1)*nqe+j       
            
              hbqpte_init(el,edpt) = 0d0          
              DO nd = 1,nnds(eo)                                 
                hbqpte_init(el,edpt) = hbqpte_init(el,edpt) + psia(nd,pt,typ)*elhb(nd,el)     
              ENDDO   
              hbqpted(ed,j) = hbqpte_init(el,edpt)
              
            ENDDO
            
!             PRINT*, (hbqpte_init(el,edpt), edpt = (led-1)*nqe+1,(led-1)*nqe+nqe)
          
          ENDDO
!           PRINT*, "  "


        ELSE IF (ed_type(ed) == 10) THEN
        
          et = el_type(el1) 
          nqa = nqpta(et)
          nqe = nqpte(et)   

          typ = et + 4          
          eo = order(typ)      
        
          DO i = 1,nqe        
            pt = nqa + (led1-1)*nqe+i         
          
            hbqpted(ed,i) = 0d0          
            DO nd = 1,nnds(eo)                                 
              hbqpted(ed,i) = hbqpted(ed,i) + psia(nd,pt,typ)*elhb(nd,el1)     
            ENDDO           
          ENDDO                            

!             PRINT("(4(F24.17))"), (hbqpted(ed,i), i = 1,nqe)   
            
        ELSE
        
          et = el_type(el1) 
          
          ! Specified flow/open edge assumed to be straight
          IF (mod(et,2) == 1) THEN
            et = 1
          ELSE IF (mod(et,2) == 0) THEN
            et = 2
          ENDIF      
          
          nqa = nqpta(et)
          nqe = nqpte(et)           

          typ = et + 4          
          eo = order(typ)      
        
          DO i = 1,nqe        
            pt = nqa + (led1-1)*nqe+i            
          
            hbqpted(ed,i) = 0d0          
            DO nd = 1,nnds(eo)                                 
              hbqpted(ed,i) = hbqpted(ed,i) + psia(nd,pt,typ)*elhb(nd,el1)     
            ENDDO           
          ENDDO                            
        
          
        ENDIF
        
      ENDDO         
      
!       edpt = 0
!       DO ed = 1,ned
!       
!         el1 = ged2el(1,ed)
!         led1 = ged2led(1,ed)          
!         el2 = ged2el(2,ed)
!         led2 = ged2led(2,ed)  
!         
!         DO pt = 1,nqpte(1)
!         
!           gp1 = (led1-1)*nqpte(1) + pt
!           gp2 = (led2-1)*nqpte(1) + nqpte(1) - pt + 1
!           
!           IF (abs(hbqpte_init(el1,gp1)- hbqpte_init(el2,gp2)) > 1d-13 .and. recv_edge(ed) == 0 ) THEN
!             PRINT*, "edge batymetry values don't match"
!             PRINT*, "recv_edge =", recv_edge(ed)
!             PRINT*, hbqpte_init(el1,gp1),hbqpte_init(el2,gp2)   
!             edpt = edpt + 1
!           ENDIF
!           
!         
!         ENDDO
! 
!       ENDDO      
!       PRINT*, "# of unmatched points", edpt
      
!       CLOSE(47)
     
      ALLOCATE(hbm(mnnds,ne))     
      hbm = 0d0 
      DO el = 1,ne
      
        et = el_type(el)
        typ = et + 4
        eo = order(typ)     
                
        DO nd = 1,nnds(eo)
          hbm(nd,el) = elhb(nd,el)
        ENDDO
        
        CALL DGETRS("T",nnds(eo),1,Va(1,1,eo),mnnds,ipiva(1,eo),hbm(1,el),mnnds,info)                         
        
      ENDDO
      
      OPEN(unit=65,FILE=trim(out_direc) //'hb_modal.d')
      DO dof = 1,mnnds
        WRITE(65,"(16000(e24.17,1x))") (hbm(dof,el), el = 1,ne)
      ENDDO
      CLOSE(65)
      


      RETURN
      END SUBROUTINE bathymetry_interp