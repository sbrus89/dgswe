      SUBROUTINE area_transformation()

      USE globals, ONLY: pres,ne,el_type,mndof,elxy, &
                         ndof,nnds,nqpta,wpta,dpdr,dpds,phia,mmi_init, &
                         psia,dpsidr,dpsids, &
                         detJa,dpdx_init,dpdy_init,phia,phia_int_init, &
                         coord_sys,r_earth,sphi0

      IMPLICIT NONE
      
      INTEGER :: i,j,nd,el,pt,m,dof
      INTEGER :: et,nnd,nqa,ndf
      INTEGER :: ipiv(mndof),work(mndof*mndof)
      INTEGER ::info
      REAL(pres) :: mm(mndof,mndof)      
      REAL(pres) :: x,y,xpt,ypt
      REAL(pres) :: dxdr,dxds,dydr,dyds
      REAL(pres) :: drdx,drdy,dsdx,dsdy
      REAL(pres) :: Sp

    
          

          
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
      
  
      RETURN
      END SUBROUTINE area_transformation
      