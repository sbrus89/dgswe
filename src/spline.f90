      PROGRAM spline

      USE globals, ONLY: rp,base,eval,ctp,nverts, &
                         ax,bx,cx,dx,ay,by,cy,dy,dt
      USE allocation, ONLY: sizes
                         
      USE calc_spline, ONLY: calc_cubic_spline,eval_cubic_spline, &
                             newton,spline_init
      USE check, ONLY: check_angle,check_deformation
      USE find_element, ONLY: in_element
      USE evaluate, ONLY: vandermonde,transformation

      IMPLICIT NONE
      INTEGER :: i,j,k,n,seg,sind,eind,num,qpts,btype
      INTEGER :: el,eln,nd,ndn,led,n1ed1,n2ed1,n1bed,n2bed,nvert
      INTEGER :: el_in,found
      INTEGER :: segtype
      INTEGER :: n1,n2
      INTEGER :: base_bed
      INTEGER :: neval,nbase
      REAL(rp) :: htest,t,tpt,x,y,xs,ys,r,sig
      REAL(rp) :: d1,d2,d3,t1,t2,xr(2),xa(2)
      REAL(rp) :: n1x,n1y,n2x,n2y,n3x,n3y,n4x,n4y,edlen
      REAL(rp) :: theta1,theta2
      



      OPEN(unit=30,file='spline.out')
      OPEN(unit=60,file='eval_nodes.out')        
      
      PRINT "(A)", " "
      
      CALL read_input()
      
      CALL sizes()

      CALL read_grid(base)
      CALL read_grid(eval)
      
      CALL connect(base)
      CALL connect(eval)
      
      sig = 0d0
      
      
      CALL spline_init(num)
      
      CALL vandermonde()  
      
      CALL transformation()        
    

      WRITE(30,*) num
      WRITE(60,*) num

      DO seg = 1,base%nbou
      
        segtype = base%fbseg(2,seg)
        
        IF( segtype == 0 .OR. segtype == 10 .OR. segtype == 20  .OR. &   ! land boundaries
            segtype == 1 .OR. segtype == 11 .OR. segtype == 21 ) THEN    ! island boundaries
        
          n = base%fbseg(1,seg)    ! n nodes, n-1 subintervals

          PRINT "(A)", " "
          PRINT "(A,I5)", "Normal flow boundary ",seg
          PRINT "(A,I5)", "Normal flow boundary nodes ",n


!           PAUSE
 


          !!!!!!!!!!!!!!!!!!!
          ! x value spline 
          !!!!!!!!!!!!!!!!!!!


          CALL calc_cubic_spline(1,seg,n,sig,ax,bx,cx,dx,dt)


          !!!!!!!!!!!!!!!!!!!
          ! y value spline 
          !!!!!!!!!!!!!!!!!!!


          CALL calc_cubic_spline(2,seg,n,sig,ay,by,cy,dy,dt)
          
                  
          
!           t = 0d0
!           DO nd = 1,n-1
! 
!            n1bed = base%fbnds(nd,seg)
!            n2bed = base%fbnds(nd+1,seg)   
!            
!            n1x = base%xy(1,n1bed)
!            n1y = base%xy(2,n1bed)
!           
!            n2x = base%xy(1,n2bed)
!            n2y = base%xy(2,n2bed)      
!            
!           
!            
!            CALL check_angle(seg,n,nd,theta1,theta2,edlen)
!            
!            
!             
!             IF ( theta1 > 20d0 .AND. theta2 > 20d0) THEN         
! 
!       elem: DO el = 1,base%nepn(n1bed)
!               eln = base%epn(el,n1bed)
!               
!               nvert = nverts(base%el_type(eln))
!               
!               DO led = 1,nvert
!  
!                 n1ed1 = base%vct(mod(led+0,nvert)+1,eln)
!                 n2ed1 = base%vct(mod(led+1,nvert)+1,eln) 
! 
!                 IF(((n1ed1 == n1bed).AND.(n2ed1 == n2bed)).OR. &
!                    ((n1ed1 == n2bed).AND.(n2ed1 == n1bed))) THEN
!                    
!                    print*, "  ", eln
!                    
!                    DO i = 1,ctp-1
!                      r = -1d0 + real(i,rp)*2d0/real(ctp,rp)
!                      tpt = .5d0*dt*(r + 1d0) + t               
!                      
!                      xm = .5d0*(1d0-r)*n1x + .5d0*(1d0+r)*n2x
!                      ym = .5d0*(1d0-r)*n1y + .5d0*(1d0+r)*n2y
!                      
!                      CALL newton(tpt,t,xm,ym,ax(nd),bx(nd),cx(nd),dx(nd),ay(nd),by(nd),cy(nd),dy(nd),x,y)
! 
!                      CALL check_deformation(base%minedlen(eln),xm,ym,x,y)
!                      
!                      ndn = mod(led,nvert)*ctp+i+1                        
!                      
!                      base%xy(1,base%ect(ndn,eln)) = x
!                      base%xy(2,base%ect(ndn,eln)) = y                                                         
!                        
!                    ENDDO
!                    
!                    t = t + dt
!                    EXIT elem
! 
!                 ENDIF
!                 
!               ENDDO
!             ENDDO elem
!             ENDIF                      
!             
!             
!           ENDDO          


          CALL write_spline(n)
          
          
  
          neval = eval%fbseg(1,seg)
  
          WRITE(60,*) ctp*(neval-1) + 1

          DO i = 1,neval-1  
             
            n1 = eval%fbnds(i,seg)
            n2 = eval%fbnds(i+1,seg)
                   
            n1x = eval%xy(1,n1)
            n1y = eval%xy(2,n1)
          
            n2x = eval%xy(1,n2)
            n2y = eval%xy(2,n2) 
            
            xa(1) = .5d0*(n1x + n2x) ! average coordinates to avoid ambiguity with verticies
            xa(2) = .5d0*(n1y + n2y)
            
            PRINT*, "FINDING ELEMENT FOR POINT: ", i, " NODE: ",n1
            CALL in_element(seg,n1,xa,base_bed)                        
            nd = base_bed    
            
            t = 0d0        ! find starting parameter value for found edge
            DO j = 1,nd-1
              t = t + dt(j)
            ENDDO
            
            IF (i == neval-1) THEN
              n = ctp
            ELSE 
              n = ctp-1
            ENDIF
                   
!             PRINT*, "EVALUATING SPLINE COORDINATES"        
            DO j = 0,n                                     
            
              r = -1d0 + real(j,rp)*2d0/real(ctp,rp)   
                     
              xr(1) = .5d0*(1d0-r)*n1x + .5d0*(1d0+r)*n2x
              xr(2) = .5d0*(1d0-r)*n1y + .5d0*(1d0+r)*n2y
              
     
              tpt = .5d0*dt(nd)*(r + 1d0) + t          ! initial guess for iteration                            
                   
              CALL newton(tpt,t,xr,ax(nd),bx(nd),cx(nd),dx(nd),ay(nd),by(nd),cy(nd),dy(nd),x,y)
!               PRINT*, 2d0/dt*(tpt-t)-1d0
              
!               CALL eval_cubic_spline(tpt,t,ax(nd),bx(nd),cx(nd),dx(nd),x)
!               CALL eval_cubic_spline(tpt,t,ay(nd),by(nd),cy(nd),dy(nd),y)              
              
              WRITE(60,*) x,y
                       
            ENDDO  

            PRINT*, "------------------------------------------------------------"            
            PRINT*, " "
            
          ENDDO

              
        ENDIF 
        

      ENDDO
      CLOSE(30)
      CLOSE(60)
      
      
      
      
      CALL write_grid(base)

      END PROGRAM spline
      
      

      
      
