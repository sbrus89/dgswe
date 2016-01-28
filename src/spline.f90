      PROGRAM spline

      USE globals, ONLY: rp,ne,nn,ctp,xy,ect,vct,depth,nelnds, &
                         nope,neta,obnds,nbou,nvel,fbseg,fbnds, &
                         nepn,epn,nverts,el_type, &
                         mnelnds,minedlen
                         
      USE calc_spline, ONLY: cubic_spline                         

      IMPLICIT NONE
      INTEGER :: i,j,k,n,seg,sind,eind,num,qpts,btype,nmax
      INTEGER :: el,eln,nd,ndn,led,n1ed1,n2ed1,n1bed,n2bed,n3bed,n4bed,nvert
      REAL(rp) :: htest,mult,dt,t,tpt,x,y,xs,ys,r,sig
      REAL(rp) :: d1,d2,d3,t1,t2,xm,ym
      REAL(rp) :: n1x,n1y,n2x,n2y,n3x,n3y,n4x,n4y,edlen
      REAL(rp) :: newton,angle,theta1,theta2
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: ax,bx,cx,dx
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: ay,by,cy,dy
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: Ml,Md,Mu,v

      OPEN(unit=30,file='spline.out')
      
      PRINT "(A)", " "

      CALL read_grid()
      
      CALL connect()
      
      sig = 1d0
      
      num = 0
      nmax = 0 
      DO seg = 1,nbou
        IF(fbseg(2,seg) == 10 .OR. fbseg(2,seg) == 11 .OR. fbseg(2,seg) == 101)THEN          
          num = num + 1
          IF (fbseg(1,seg) > nmax) THEN
            nmax = fbseg(1,seg)
          ENDIF
        ENDIF
      ENDDO

      PRINT "(A)", " "
      PRINT "(A,I5)", "Total number of type 0 normal flow boundaries ",num
      PRINT "(A,I5)", "Max number of nodes in a flow boundary segment ",nmax
      PRINT "(A)", " "
      
      ALLOCATE(ax(nmax),cx(nmax),bx(nmax-1),dx(nmax-1))
      ALLOCATE(ay(nmax),cy(nmax),by(nmax-1),dy(nmax-1))
      ALLOCATE(Ml(nmax),Md(nmax),Mu(nmax),v(nmax))      

      WRITE(30,*) num

      DO seg = 1,nbou

        IF(fbseg(2,seg) == 10 .OR. fbseg(2,seg) == 11 .OR. fbseg(2,seg) == 101)THEN
          n = fbseg(1,seg)    ! n nodes, n-1 subintervals

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
            ax(k) = xy(1,fbnds(i,seg))
            k = k+1
          ENDDO

          CALL cubic_spline(sig,n,ax,bx,cx,dx,dt)


          !!!!!!!!!!!!!!!!!!!
          ! y value spline 
          !!!!!!!!!!!!!!!!!!!

          ! Load nodal boundary y coordinates
          k = 1
          DO i = 1,n
            ay(k) = xy(2,fbnds(i,seg))
            k = k+1
          ENDDO

          CALL cubic_spline(sig,n,ay,by,cy,dy,dt)
          
          
          PRINT*, dt          
          
          t = 0d0
          DO nd = 1,n-1
            n1bed = fbnds(nd,seg)
            n2bed = fbnds(nd+1,seg) 
            IF (nd == n-1) THEN
              n3bed = fbnds(2,seg)
            ELSE
              n3bed = fbnds(nd+2,seg)            
            ENDIF
            
            IF (nd == 1) THEN
              n4bed = fbnds(n-1,seg)
            ELSE
              n4bed = fbnds(nd-1,seg)            
            ENDIF           
            
            n1x = xy(1,n1bed)
            n1y = xy(2,n1bed)
            
            n2x = xy(1,n2bed)
            n2y = xy(2,n2bed)
            
            n3x = xy(1,n3bed)
            n3y = xy(2,n3bed)            
            
            n4x = xy(1,n4bed)
            n4y = xy(2,n4bed)            
            
            theta1 = angle(n1x,n1y,n2x,n2y,n3x,n3y)
            theta2 = angle(n4x,n4y,n1x,n1y,n2x,n2y) 
            
            PRINT*, theta1,theta2
            

            edlen = sqrt((n1x-n2x)**2+(n1y-n2y)**2)       
            
            IF ( theta1 > 20d0 .AND. theta2 > 20d0) THEN         

      elem: DO el = 1,nepn(n1bed)
              eln = epn(el,n1bed)
              
              nvert = nverts(el_type(eln))
              
              DO led = 1,nvert
 
                n1ed1 = vct(mod(led+0,nvert)+1,eln)
                n2ed1 = vct(mod(led+1,nvert)+1,eln) 

                IF(((n1ed1 == n1bed).AND.(n2ed1 == n2bed)).OR. &
                   ((n1ed1 == n2bed).AND.(n2ed1 == n1bed))) THEN
                   
                   print*, "  ", eln
                   
                   DO i = 1,ctp-1
                     r = -1d0 + real(i,rp)*2d0/real(ctp,rp)
                     tpt = .5d0*dt*(r + 1d0) + t               
                     
                     xm = .5d0*(1d0-r)*n1x + .5d0*(1d0+r)*n2x
                     ym = .5d0*(1d0-r)*n1y + .5d0*(1d0+r)*n2y
                     
                     tpt = newton(tpt,t,xm,ym,ax(nd),bx(nd),cx(nd),dx(nd),ay(nd),by(nd),cy(nd),dy(nd))

                     x = ax(nd) + bx(nd)*(tpt-t) + cx(nd)*(tpt-t)**2 + dx(nd)*(tpt-t)**3
                     y = ay(nd) + by(nd)*(tpt-t) + cy(nd)*(tpt-t)**2 + dy(nd)*(tpt-t)**3
                     
                     xs = x
                     ys = y
                     
                     d1 = sqrt((xm-x)**2 + (ym-y)**2)
                                          
                     
                     ndn = mod(led,nvert)*ctp+i+1   
                     
                     r = -1d0
                     DO WHILE (d1 > .1d0*minedlen(eln) .AND. r < 1d0)
                       r = r + .1d0
                       
                       x = .5d0*(1d0-r)*xs + .5d0*(1d0+r)*xm
                       y = .5d0*(1d0-r)*ys + .5d0*(1d0+r)*ym
                       
                       d1 = sqrt((xm-x)**2 + (ym-y)**2)
                     ENDDO

                     IF (d1 <= .1d0*minedlen(eln)) THEN
                       xy(1,ect(ndn,eln)) = x
                       xy(2,ect(ndn,eln)) = y
                     ELSE
                       PRINT*, "DEFORMATION TOO LARGE"
                     ENDIF                       
                       
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
      
      CALL write_grid()

      END PROGRAM spline
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      REAL(kind(0d0)) FUNCTION newton(told,t,xm,ym,ax,bx,cx,dx,ay,by,cy,dy)
      
      IMPLICIT NONE      
      
      INTEGER, PARAMETER :: rp = kind(0d0)
      INTEGER :: it,maxit
      REAL(rp) :: told,tnew,t
      REAL(rp) :: xm,ym
      REAL(rp) :: ax,bx,cx,dx,ay,by,cy,dy
      REAL(rp) :: f,fp,fpp,g,gp,gpp
      REAL(rp) :: d,dp,tol
      
      tol = 1d-4
      maxit = 10000
      
iter: DO it = 1,maxit

        f = ax + bx*(told-t) + cx*(told-t)**2 + dx*(told-t)**3
        fp = bx + 2d0*cx*(told-t) + 3d0*dx*(told-t)**2
        fpp = 2d0*cx + 6d0*dx*(told-t)
        
        g = ay + by*(told-t) + cy*(told-t)**2 + dy*(told-t)**3
        gp = by + 2d0*cy*(told-t) + 3d0*dy*(told-t)**2
        gpp = 2d0*cy + 6d0*dy*(told-t)
        
        d = 2d0*(f-xm)*fp + 2d0*(g-ym)*gp
        dp = 2d0*fp**2 + 2d0*(f-xm)*fpp + 2d0*gp**2 + 2d0*(g-ym)*gpp
        
        tnew = told - d/dp
        
        IF (ABS(d) < tol) THEN
!           PRINT*, "iterations", it
          EXIT iter
        ENDIF
        
        told = tnew
        
      ENDDO iter
      
      IF (it >= maxit) THEN
        PRINT*, "MAX ITERATIONS EXCEEDED"
      ENDIF     
      
      newton = tnew
      
      RETURN 
      END FUNCTION
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      REAL(kind(0d0)) FUNCTION angle(x1,y1,x2,y2,x3,y3)
      
      IMPLICIT NONE  
      
      INTEGER, PARAMETER :: rp = kind(0d0)
      REAL(rp), PARAMETER  ::  pi=3.141592653589793D0
      REAL(rp) :: x1,x2,x3,xy,y1,y2,y3,y4
      REAL(rp) :: xv1,xv2,yv1,yv2
      REAL(rp) :: adotb,la,lb
      
      
      xv1 = x1-x2
      yv1 = y1-y2
      
      xv2 = x3-x2
      yv2 = y3-y2      
      
      adotb = xv1*xv2 + yv1*yv2
      
      la = sqrt(xv1**2 + yv1**2)
      lb = sqrt(xv2**2 + yv2**2)
      
      angle = acos(adotb/(la*lb))*(180d0/pi)
      
      RETURN
      END FUNCTION