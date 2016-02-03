      PROGRAM spline

      USE globals, ONLY: rp,base,eval,ctp,nverts, &
                         ax,bx,cx,dx,ay,by,cy,dy
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
      REAL(rp) :: htest,dt,t,tpt,x,y,xs,ys,r,sig
      REAL(rp) :: d1,d2,d3,t1,t2,xm,ym
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
          PRINT "(A)", " "

!           PAUSE
 


          !!!!!!!!!!!!!!!!!!!
          ! x value spline 
          !!!!!!!!!!!!!!!!!!!

          ! Load nodal boundary x coordinates
          k = 1   
          DO i = 1,n
            ax(k) = base%xy(1,base%fbnds(i,seg))
            k = k+1
          ENDDO

          CALL calc_cubic_spline(sig,n,ax,bx,cx,dx,dt)


          !!!!!!!!!!!!!!!!!!!
          ! y value spline 
          !!!!!!!!!!!!!!!!!!!

          ! Load nodal boundary y coordinates
          k = 1
          DO i = 1,n
            ay(k) = base%xy(2,base%fbnds(i,seg))
            k = k+1
          ENDDO

          CALL calc_cubic_spline(sig,n,ay,by,cy,dy,dt)
          
          
          PRINT*, dt          
          
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


          CALL write_spline(n,dt)
          
          
  
          neval = eval%fbseg(1,seg)
          nbase = base%fbseg(1,seg)
  
          WRITE(60,*) ctp*(neval-1) + 1

          DO i = 1,neval-1  
            nd = eval%fbnds(i,seg)
             
            PRINT*, "FINDING ELEMENT FOR POINT: ",i, " NODE: ",nd
            CALL in_element(nd,eval%xy(1:2,nd),el_in,found)          

            
            nvert = nverts(base%el_type(el_in))            
            
            found = 0
     ledge: DO led = 1,nvert
 
              n1ed1 = base%vct(mod(led+0,nvert)+1,el_in)
              n2ed1 = base%vct(mod(led+1,nvert)+1,el_in) 
              
              DO j = 1,nbase-1
              
                n1bed = base%fbnds(j,seg)
                n2bed = base%fbnds(j+1,seg)   
                                                        
                IF(((n1ed1 == n1bed).AND.(n2ed1 == n2bed)).OR. &
                   ((n1ed1 == n2bed).AND.(n2ed1 == n1bed))) THEN
!                    PRINT*, "n1bed = ",n1bed, "n2bed = ",n2bed                   
                   found = 1
                   
                   t = real(j-1,rp)*dt
                   base_bed = j
                   
                   EXIT ledge
                ENDIF        
              
              ENDDO
              
            ENDDO ledge
            
            IF (found == 0) THEN
              PRINT*, "Boundary edge not found"
            ENDIF           
                   
                   
            n1 = eval%fbnds(i,seg)
            n2 = eval%fbnds(i+1,seg)
                   
            n1x = eval%xy(1,n1)
            n1y = eval%xy(2,n1)
          
            n2x = eval%xy(1,n2)
            n2y = eval%xy(2,n2)     
            
            IF (i == neval-1) THEN
              n = ctp
            ELSE 
              n = ctp-1
            ENDIF
                   
            DO j = 0,n
              r = -1d0 + real(j,rp)*2d0/real(ctp,rp)
              tpt = .5d0*dt*(r + 1d0) + t               
                     
              xm = .5d0*(1d0-r)*n1x + .5d0*(1d0+r)*n2x
              ym = .5d0*(1d0-r)*n1y + .5d0*(1d0+r)*n2y
                     
              nd = base_bed                     
              CALL newton(tpt,t,xm,ym,ax(nd),bx(nd),cx(nd),dx(nd),ay(nd),by(nd),cy(nd),dy(nd),x,y)
              
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
      
      

      
      
