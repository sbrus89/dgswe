      PROGRAM spline

      USE globals, ONLY: rp,base,eval,ctp,nverts
      USE allocation, ONLY: sizes
                         
      USE calc_spline, ONLY: calc_cubic_spline,eval_cubic_spline,newton
      USE check, ONLY: check_angle,check_deformation

      IMPLICIT NONE
      INTEGER :: i,j,k,n,seg,sind,eind,num,qpts,btype,nmax
      INTEGER :: el,eln,nd,ndn,led,n1ed1,n2ed1,n1bed,n2bed,nvert
      REAL(rp) :: htest,dt,t,tpt,x,y,xs,ys,r,sig
      REAL(rp) :: d1,d2,d3,t1,t2,xm,ym
      REAL(rp) :: n1x,n1y,n2x,n2y,n3x,n3y,n4x,n4y,edlen
      REAL(rp) :: theta1,theta2
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: ax,bx,cx,dx
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: ay,by,cy,dy


      OPEN(unit=30,file='spline.out')
      
      PRINT "(A)", " "
      
      CALL read_input()
      
      CALL sizes()

      CALL read_grid(base)
      CALL read_grid(eval)
      
      CALL connect(base)
      CALL connect(eval)
      
      sig = 1d0
      
      num = 0
      nmax = 0 
      DO seg = 1,base%nbou
        IF(base%fbseg(2,seg) == 10 .OR. base%fbseg(2,seg) == 11 .OR. base%fbseg(2,seg) == 101)THEN          
          num = num + 1
          IF (base%fbseg(1,seg) > nmax) THEN
            nmax = base%fbseg(1,seg)
          ENDIF
        ENDIF
      ENDDO

      PRINT "(A)", " "
      PRINT "(A,I5)", "Total number of type 0 normal flow boundaries ",num
      PRINT "(A,I5)", "Max number of nodes in a flow boundary segment ",nmax
      PRINT "(A)", " "
      
      ALLOCATE(ax(nmax),cx(nmax),bx(nmax-1),dx(nmax-1))
      ALLOCATE(ay(nmax),cy(nmax),by(nmax-1),dy(nmax-1))
    

      WRITE(30,*) num

      DO seg = 1,base%nbou

        IF(base%fbseg(2,seg) == 10 .OR. base%fbseg(2,seg) == 11 .OR. base%fbseg(2,seg) == 101)THEN
        
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
          
          t = 0d0
          DO nd = 1,n-1

           n1bed = base%fbnds(nd,seg)
           n2bed = base%fbnds(nd+1,seg)   
           
           n1x = base%xy(1,n1bed)
           n1y = base%xy(2,n1bed)
          
           n2x = base%xy(1,n2bed)
           n2y = base%xy(2,n2bed)      
           
          
           
           CALL check_angle(seg,n,nd,theta1,theta2,edlen)
           
           
            
            IF ( theta1 > 20d0 .AND. theta2 > 20d0) THEN         

      elem: DO el = 1,base%nepn(n1bed)
              eln = base%epn(el,n1bed)
              
              nvert = nverts(base%el_type(eln))
              
              DO led = 1,nvert
 
                n1ed1 = base%vct(mod(led+0,nvert)+1,eln)
                n2ed1 = base%vct(mod(led+1,nvert)+1,eln) 

                IF(((n1ed1 == n1bed).AND.(n2ed1 == n2bed)).OR. &
                   ((n1ed1 == n2bed).AND.(n2ed1 == n1bed))) THEN
                   
                   print*, "  ", eln
                   
                   DO i = 1,ctp-1
                     r = -1d0 + real(i,rp)*2d0/real(ctp,rp)
                     tpt = .5d0*dt*(r + 1d0) + t               
                     
                     xm = .5d0*(1d0-r)*n1x + .5d0*(1d0+r)*n2x
                     ym = .5d0*(1d0-r)*n1y + .5d0*(1d0+r)*n2y
                     
                     CALL newton(tpt,t,xm,ym,ax(nd),bx(nd),cx(nd),dx(nd),ay(nd),by(nd),cy(nd),dy(nd),x,y)

                     CALL check_deformation(base%minedlen(eln),xm,ym,x,y)
                     
                     ndn = mod(led,nvert)*ctp+i+1                        
                     
                     base%xy(1,base%ect(ndn,eln)) = x
                     base%xy(2,base%ect(ndn,eln)) = y                                                         
                       
                   ENDDO
                   
                   t = t + dt
                   EXIT elem

                ENDIF
                
              ENDDO
            ENDDO elem
            ENDIF
            
            
          ENDDO          

  
          WRITE(30,*) n,dt

          DO i = 1,n-1
            WRITE(30,"(4(E25.12))"), ax(i),bx(i),cx(i),dx(i)
          ENDDO
          WRITE(30,"(4(E25.12))"), ax(n),0.0,0.0,0.0

          DO i = 1,n-1
            WRITE(30,"(4(E25.12))"), ay(i),by(i),cy(i),dy(i)
          ENDDO
          WRITE(30,"(4(E25.12))"), ay(n),0.0,0.0,0.0


        ENDIF

      ENDDO
      CLOSE(30)
      
      
      
      
      CALL write_grid(base)

      END PROGRAM spline
      
      

      
      
