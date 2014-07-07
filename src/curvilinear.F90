      SUBROUTINE curvilinear()

      USE globals, ONLY: pres,nnfbed,nfbedn,ged2el,ged2led,nelnds,el_type,ctp,elxy,xy,ect,elhb
      USE basis, ONLY: quad_nodes,tri_nodes

      IMPLICIT NONE
      
      INTEGER :: ed,led,nd,pt
      INTEGER :: ged,el,ind
      INTEGER :: nvert
      REAL(pres) :: x,y,ytest,hb
      REAL(pres) :: rq((ctp+1)**2),sq((ctp+1)**2)
      REAL(pres) :: rt((ctp+1)*(ctp+2)/2),st((ctp+1)*(ctp+2)/2)      
      
      CALL tri_nodes(1,ctp,(ctp+1)*(ctp+2)/2,rt,st)   
      
      CALL quad_nodes(1,ctp,(ctp+1)*(ctp+1),rq,sq)      
      
      
      DO ed = 1,nnfbed
        ged = nfbedn(ed)
        
        el = ged2el(1,ged)
        led = ged2led(1,ged)
        
        IF (nelnds(el) == 3) THEN
          el_type(el) = 3 
          nelnds(el) = (ctp+1)*(ctp+2)/2
          nvert = 3
                   
          DO nd = 1,nelnds(el)            
              
            x = .5d0*(-(rt(nd)+st(nd))*xy(1,ect(1,el)) + (rt(nd)+1d0)*xy(1,ect(2,el)) + (st(nd)+1d0)*xy(1,ect(3,el)))
            y = .5d0*(-(rt(nd)+st(nd))*xy(2,ect(1,el)) + (rt(nd)+1d0)*xy(2,ect(2,el)) + (st(nd)+1d0)*xy(2,ect(3,el)))            
            
            elxy(nd,el,1) = x
            elxy(nd,el,2) = y
          ENDDO
        
        ELSE IF (nelnds(el) == 4) THEN
          el_type(el) = 4
          nelnds(el) = (ctp+1)*(ctp+1)
          nvert = 4
                    
          DO nd = 1,nelnds(el)
          
            x = .25d0*((1d0-rq(nd))*(1d0-sq(nd))*xy(1,ect(1,el)) + (1d0+rq(nd))*(1d0-sq(nd))*xy(1,ect(2,el)) + (1d0+rq(nd))*(1d0+sq(nd))*xy(1,ect(3,el)) + (1d0-rq(nd))*(1d0+sq(nd))*xy(1,ect(4,el)))
            y = .25d0*((1d0-rq(nd))*(1d0-sq(nd))*xy(2,ect(1,el)) + (1d0+rq(nd))*(1d0-sq(nd))*xy(2,ect(2,el)) + (1d0+rq(nd))*(1d0+sq(nd))*xy(2,ect(3,el)) + (1d0-rq(nd))*(1d0+sq(nd))*xy(2,ect(4,el)))
            
            elxy(nd,el,1) = x
            elxy(nd,el,2) = y
          ENDDO          
        ENDIF
        
        DO nd = 1,ctp-1
          pt = mod(led,nvert)*ctp + 1 + nd
          
          ytest = elxy(pt,el,2)
          x = elxy(pt,el,1) 
          
          IF (ytest < 250d0) THEN
            y = 0d0 + 100d0*(1d0/(COSH(4d0*(x-2000d0)/500d0)))
          ELSE IF (ytest > 250d0) THEN
            y = 500d0 - 100d0*(1d0/(COSH(4d0*(x-2000d0)/500d0)))
          ENDIF
          
          elxy(pt,el,2) = y
        ENDDO
        
        DO nd = 1,nelnds(el)
          elhb(nd,el) = 10d0
        ENDDO
        
      ENDDO
      
      
      RETURN
      END SUBROUTINE curvilinear