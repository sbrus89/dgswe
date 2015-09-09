      SUBROUTINE curvilinear()

      USE globals, ONLY: pres,nnfbed,nfbedn,ged2el,ged2led,pi, &
                         nelnds,el_type,order,elxy,xy,ect,elhb,ne,mnnds,nnds, &
                         psiv,psic
      USE basis, ONLY: quad_nodes,tri_nodes
      USE read_dginp, ONLY: ctp
      USE transformation, ONLY: element_transformation

      IMPLICIT NONE
      
      INTEGER :: ed,led,nd,pt
      INTEGER :: ged,el,ind
      INTEGER :: et,typ,eo
      INTEGER :: nvert,nnd,npts,nt,nq
      REAL(pres) :: xpt,ypt,ytest,hb
      REAL(pres) :: x(mnnds),y(mnnds)      
      REAL(pres) :: rq(mnnds),sq(mnnds)
      REAL(pres) :: rt(mnnds),st(mnnds)  
      REAL(pres) :: xyhb(mnnds,ne,2)
      
      CALL tri_nodes(1,ctp,nt,rt,st)   
      
      CALL quad_nodes(1,ctp,nq,rq,sq)   
      
      elhb = 50d0
      
      DO el = 1,ne
      
        et = el_type(el)
        nnd = nnds(et)        
        IF (mod(et,2) == 1) THEN
          npts = nnds(5)
          et = 5
        ELSE IF (mod(et,2) == 0) THEN
          npts = nnds(6)
          et = 6
        ENDIF
        
        DO pt = 1,npts    
          
          DO nd = 1,nnd
            x(nd) = xy(1,ect(nd,el))
            y(nd) = xy(2,ect(nd,el))
          ENDDO

          CALL element_transformation(nnd,x,y,psiv(:,pt,et),xpt,ypt)
          
          elhb(pt,el) = 10d0
!           elhb(pt,el) = 10d0 - 5d0*cos(2d0*pi/500d0*ypt)        
          xyhb(pt,el,1) = xpt
          xyhb(pt,el,2) = ypt
        ENDDO  
        
      ENDDO      
      
      
      
      
      
      DO ed = 1,nnfbed
        ged = nfbedn(ed)
        
        el = ged2el(1,ged)
        led = ged2led(1,ged)
        
        IF (nelnds(el) == 3) THEN
          el_type(el) = 3 
          nnd = nnds(3)
          nvert = 3

                   
!           DO nd = 1,nnd                           
!             x = .5d0*(-(rt(nd)+st(nd))*xy(1,ect(1,el))  &
!                      + (rt(nd)+1d0)*xy(1,ect(2,el))  & 
!                      + (st(nd)+1d0)*xy(1,ect(3,el)))
!             y = .5d0*(-(rt(nd)+st(nd))*xy(2,ect(1,el))  &
!                      + (rt(nd)+1d0)*xy(2,ect(2,el))  &
!                      + (st(nd)+1d0)*xy(2,ect(3,el)))            
!             
!             elxy(nd,el,1) = x
!             elxy(nd,el,2) = y
!           ENDDO
        
        ELSE IF (nelnds(el) == 4) THEN
          el_type(el) = 4
          nnd = nnds(4)
          nvert = 4
                    
!           DO nd = 1,nelnds(el)
!           
!             x = .25d0*((1d0-rq(nd))*(1d0-sq(nd))*xy(1,ect(1,el)) &
!                      + (1d0+rq(nd))*(1d0-sq(nd))*xy(1,ect(2,el)) &
!                      + (1d0+rq(nd))*(1d0+sq(nd))*xy(1,ect(3,el)) &
!                      + (1d0-rq(nd))*(1d0+sq(nd))*xy(1,ect(4,el)))
!             y = .25d0*((1d0-rq(nd))*(1d0-sq(nd))*xy(2,ect(1,el)) &
!                      + (1d0+rq(nd))*(1d0-sq(nd))*xy(2,ect(2,el)) &
!                      + (1d0+rq(nd))*(1d0+sq(nd))*xy(2,ect(3,el)) &
!                      + (1d0-rq(nd))*(1d0+sq(nd))*xy(2,ect(4,el)))
!             
!             elxy(nd,el,1) = x
!             elxy(nd,el,2) = y
!           ENDDO  

        ENDIF
        
        et = el_type(el)        
        nelnds(el) = nnd  
        
        
        DO pt = 1,nnd     
          
          DO nd = 1,nvert
            x(nd) = xy(1,ect(nd,el))
            y(nd) = xy(2,ect(nd,el))
          ENDDO

          CALL element_transformation(nvert,x,y,psiv(:,pt,et),xpt,ypt)
          
          elxy(pt,el,1) = xpt
          elxy(pt,el,2) = ypt
        ENDDO      
        
        
        
        DO nd = 1,ctp-1
          pt = mod(led,nvert)*ctp + 1 + nd
          
          ytest = elxy(pt,el,2)
          xpt = elxy(pt,el,1) 
          
          IF (ytest < 250d0) THEN
            ypt = 0d0 + 100d0*(1d0/(COSH(4d0*(xpt-2000d0)/500d0)))
          ELSE IF (ytest > 250d0) THEN
            ypt = 500d0 - 100d0*(1d0/(COSH(4d0*(xpt-2000d0)/500d0)))
          ENDIF
          
          elxy(pt,el,2) = ypt
        ENDDO       
        
        
        
        
        et = el_type(el)
        nnd = nnds(et)        
        IF (mod(et,2) == 1) THEN
          npts = nnds(5)
          et = 5
        ELSE IF (mod(et,2) == 0) THEN
          npts = nnds(6)
          et = 6
        ENDIF
        
        DO pt = 1,npts              

          CALL element_transformation(nnd,elxy(:,el,1),elxy(:,el,2),psic(:,pt,et),xpt,ypt)
          
          elhb(pt,el) = 10d0
!           elhb(pt,el) = 10d0 - 5d0*cos(2d0*pi/500d0*ypt)     
          xyhb(pt,el,1) = xpt
          xyhb(pt,el,2) = ypt          
        ENDDO  
        
        
      ENDDO
      
      OPEN(unit=242,file='bathy.d')
      DO el = 1,ne
      
        et = el_type(el)     
        IF (mod(et,2) == 1) THEN
          npts = nnds(5)
        ELSE IF (mod(et,2) == 0) THEN
          npts = nnds(6)
        ENDIF
        
        DO pt = 1,npts
          WRITE(242,"(3(e24.17,1x))") xyhb(pt,el,1),xyhb(pt,el,2),elhb(pt,el)
        ENDDO
      ENDDO
      CLOSE(242)

      
      
      RETURN
      END SUBROUTINE curvilinear