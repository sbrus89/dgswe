      PROGRAM spline

      USE globals, ONLY: rp,base,eval,ctp,nverts, &
                         ax,bx,cx,dx,ay,by,cy,dy,dt, &
                         rpts,theta_tol,sig
      USE allocation, ONLY: sizes
                         
      USE calc_spline, ONLY: calc_cubic_spline,eval_cubic_spline, &
                             newton,spline_init
      USE check, ONLY: check_angle,check_deformation,l2_project,quad_interp
      USE find_element, ONLY: in_element,check_elem
      USE evaluate, ONLY: vandermonde,transformation  

      IMPLICIT NONE
      INTEGER :: i,j,k,n,seg,sind,eind,num,qpts,btype
      INTEGER :: el,eln,nd,ndn,led,n1ed1,n2ed1,n1bed,n2bed,nvert
      INTEGER :: el_in,found
      INTEGER :: segtype
      INTEGER :: n1,n2
      INTEGER :: base_bed
      INTEGER :: neval,nbase
      REAL(rp) :: htest,ti,tpt,x,y,r,ra,xs,ys
      REAL(rp) :: d1,d2,d3,t1,t2,xr(2),xa(2)
      REAL(rp) :: n1x,n1y,n2x,n2y,n3x,n3y,n4x,n4y,edlen
      REAL(rp) :: theta1,theta2

      



      OPEN(unit=30,file='spline.out')
      OPEN(unit=60,file='eval_nodes.out') 
      OPEN(unit=90,file='max_deform.out')
      
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
      WRITE(90,*) num

      DO seg = 1,base%nbou
      
        segtype = base%fbseg(2,seg)
        
        IF( segtype == 0 .OR. segtype == 10 .OR. segtype == 20  .OR. &   ! land boundaries
            segtype == 1 .OR. segtype == 11 .OR. segtype == 21 ) THEN    ! island boundaries
        
          n = base%fbseg(1,seg)    ! n nodes, n-1 subintervals

          PRINT "(A)", " "          
          PRINT "(A)", " "          
          PRINT "(A)", " "          
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
          
                  
           
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! modify spline if necessary
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          WRITE(90,*) n-1
          
          DO nd = 1,n-1
          
            ti = 0d0        ! find starting parameter value for found edge
            DO j = 1,nd-1
              ti = ti + dt(j)
            ENDDO            
                                            
            CALL check_angle(seg,nd,dt(nd),ti,ax(nd),bx(nd),cx(nd),dx(nd),ay(nd),by(nd),cy(nd),dy(nd))            
                                              
            CALL check_deformation(seg,nd,dt(nd),ti,ax(nd),bx(nd),cx(nd),dx(nd),ay(nd),by(nd),cy(nd),dy(nd))                           
                        
          ENDDO          



          CALL write_spline(n)
          
          
          
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! evaluate spline at eval grid points
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         
  
          neval = eval%fbseg(1,seg)
  
          WRITE(60,*) ctp*(neval-1) + 1

          DO i = 1,neval-1  
             
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
                   
!             PRINT*, "EVALUATING SPLINE COORDINATES"        
            DO j = 0,n                                     
            
              r = rpts(j+1)   
              
              IF (j == 0) THEN
                ra = r + 1d-2          ! add an offset avoid ambiguity with verticies
              ELSE IF (j == ctp) THEN
                ra = r - 1d-2
              ELSE 
                ra = r
              ENDIF
                                                            
              xa(1) = .5d0*(1d0-ra)*n1x + .5d0*(1d0+ra)*n2x 
              xa(2) = .5d0*(1d0-ra)*n1y + .5d0*(1d0+ra)*n2y                                  
              
              PRINT*, "FINDING ELEMENT FOR POINT: ", i, " NODE: ",n1
              CALL in_element(seg,n1,n2,xa,el_in,base_bed)                        
              nd = base_bed   
              
              xr(1) = .5d0*(1d0-r)*n1x + .5d0*(1d0+r)*n2x
              xr(2) = .5d0*(1d0-r)*n1y + .5d0*(1d0+r)*n2y                
              
              ! check to make sure vertex offset used to find element isn't 
              ! too large that the found element isn't connected to the vertex point
              IF (j == 0 .OR. j == ctp) THEN           
                CALL check_elem(xr,el_in)         
              ENDIF               
            
              ti = 0d0        ! find starting parameter value for found edge
              DO k = 1,nd-1
                ti = ti + dt(k)
              ENDDO              
                      
                   
              CALL newton(r,dt(nd),ti,xr,ax(nd),bx(nd),cx(nd),dx(nd),ay(nd),by(nd),cy(nd),dy(nd),x,y)
              

!               ! Try new initial guess if minimum was not found in (-1,1) interval              
!               IF (r < -1d0) THEN
!                 PRINT "(A,F24.17)", "R = ", r
!                 r = -1d0
!                 CALL newton(r,dt(nd),ti,xr,ax(nd),bx(nd),cx(nd),dx(nd),ay(nd),by(nd),cy(nd),dy(nd),x,y)
!               ENDIF
              
!               ! Evaluate spline at specified parameter value (no distance minimixation)              
!               tpt = .5d0*dt(nd)*(r + 1d0) + ti               
!               CALL eval_cubic_spline(tpt,ti,ax(nd),bx(nd),cx(nd),dx(nd),x)
!               CALL eval_cubic_spline(tpt,ti,ay(nd),by(nd),cy(nd),dy(nd),y)              
              
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
      
      

      
      
