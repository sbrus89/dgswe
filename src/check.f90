      MODULE check

      USE globals, ONLY:rp,pi

      IMPLICIT NONE

      CONTAINS
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE check_angle(seg,nd,dt,t,ax,bx,cx,dx,ay,by,cy,dy)
      
      USE globals, ONLY: base,theta_tol
      
      IMPLICIT NONE     
      
      INTEGER, INTENT(IN) :: nd,seg
      REAL(rp), INTENT(IN) :: dt,t
      REAL(rp), INTENT(INOUT) :: ax,bx,cx,dx,ay,by,cy,dy      
      
      INTEGER :: n
      INTEGER :: n1bed,n2bed,n3bed,n4bed      
      REAL(rp) :: n1x,n1y,n2x,n2y,n3x,n3y,n4x,n4y
      REAL(rp) :: theta1,theta2   
      REAL(rp) :: edlen

      n = base%fbseg(1,seg)    ! n nodes, n-1 subintervals
      
      n1bed = base%fbnds(nd,seg)
      n2bed = base%fbnds(nd+1,seg) 
      IF (nd == n-1) THEN
        n3bed = base%fbnds(2,seg)
      ELSE
        n3bed = base%fbnds(nd+2,seg)            
      ENDIF
            
      IF (nd == 1) THEN
        n4bed = base%fbnds(n-1,seg)
      ELSE
        n4bed = base%fbnds(nd-1,seg)            
      ENDIF           
            
      n1x = base%xy(1,n1bed)
      n1y = base%xy(2,n1bed)
          
      n2x = base%xy(1,n2bed)
      n2y = base%xy(2,n2bed)
         
      n3x = base%xy(1,n3bed)
      n3y = base%xy(2,n3bed)            
           
      n4x = base%xy(1,n4bed)
      n4y = base%xy(2,n4bed)            
          
      theta1 = angle(n1x,n1y,n2x,n2y,n3x,n3y)
      theta2 = angle(n4x,n4y,n1x,n1y,n2x,n2y) 
            
      PRINT*, theta1,theta2
            

      edlen = sqrt((n1x-n2x)**2+(n1y-n2y)**2) 
      
      IF ( theta1 < theta_tol .OR. theta2 < theta_tol) THEN               

!       CALL l2_project(dt,ax,bx,cx,dx,ay,by,cy,dy)              
        CALL quad_interp(nd,seg,dt,t,ax,bx,cx,dx,ay,by,cy,dy)
              
      ENDIF                
            
      END SUBROUTINE check_angle
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE check_deformation(seg,nd,dt,ti,ax,bx,cx,dx,ay,by,cy,dy)
      
      USE globals, ONLY: base
      USE find_element, ONLY: in_element
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: seg,nd
      REAL(rp), INTENT(IN) :: dt,ti
      REAL(rp), INTENT(INOUT) :: ax,bx,cx,dx,ay,by,cy,dy
      REAL(rp):: xd
      INTEGER :: n1,n2
      REAL(rp) :: xn1(2),xn2(2)
      REAL(rp) :: xa(2)
      REAL(rp) :: xs,ys
      REAL(rp) :: r,d1
      INTEGER :: el_in,bed
        

        n1 = base%fbnds(nd,seg)
        n2 = base%fbnds(nd+1,seg) 

        xn1(1) = base%xy(1,n1)
        xn1(2) = base%xy(2,n1)
        
        xn2(1) = base%xy(1,n2)
        xn2(2) = base%xy(2,n2)
        
        xa(1) = .5d0*(xn1(1)+xn2(1))
        xa(2) = .5d0*(xn1(2)+xn2(2))
      
        CALL in_element(seg,n1,xa,el_in,bed)  
        CALL max_diff(dt,ti,xn1,xn2,ax,bx,cx,dx,ay,by,cy,dy,xd)
        
        IF (xd > .1d0*base%minedlen(el_in)) THEN         
          PRINT*, "DEFORMATION TOO LARGE"
          PAUSE
        ENDIF
      
      RETURN
      END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      REAL(rp) FUNCTION angle(x1,y1,x2,y2,x3,y3)
      
      IMPLICIT NONE  
      
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      SUBROUTINE l2_project(dt,ax,bx,cx,dx,ay,by,cy,dy)
      
      IMPLICIT NONE
      
      REAL(rp), INTENT(IN) :: dt
      REAL(rp), INTENT(INOUT) :: ax,bx,cx,dx
      REAL(rp), INTENT(INOUT) :: ay,by,cy,dy
      
      REAL(rp) :: Al2(3,3),Bl2(3,2)
      INTEGER :: ipiv(3),info      
      
      Al2(1,1) = dt
      Al2(1,2) = dt**2/2d0
      Al2(1,3) = dt**3/3d0
              
      Al2(2,1) = dt**2/2d0
      Al2(2,2) = dt**3/3d0
      Al2(2,3) = dt**4/4d0
              
      Al2(3,1) = dt**3/3d0              
      Al2(3,2) = dt**4/4d0
      Al2(3,3) = dt**5/5d0
              
      Bl2(1,1) = ax*dt        + bx*dt**2/2d0 + cx*dt**3/3d0 + dx*dt**4/4d0
      Bl2(2,1) = ax*dt**2/2d0 + bx*dt**3/3d0 + cx*dt**4/4d0 + dx*dt**5/5d0
      Bl2(3,1) = ax*dt**3/3d0 + bx*dt**4/4d0 + cx*dt**5/5d0 + dx*dt**6/6d0
              
      Bl2(1,2) = ay*dt        + by*dt**2/2d0 + cy*dt**3/3d0 + dy*dt**4/4d0
      Bl2(2,2) = ay*dt**2/2d0 + by*dt**3/3d0 + cy*dt**4/4d0 + dy*dt**5/5d0
      Bl2(3,2) = ay*dt**3/3d0 + by*dt**4/4d0 + cy*dt**5/5d0 + dy*dt**6/6d0              
              
              
      CALL DGESV(3, 2, Al2, 3, ipiv, Bl2, 3, info)
      
              
      ax = Bl2(1,1)
      bx = Bl2(2,1)
      cx = Bl2(3,1)
      dx = 0d0
              
      ay = Bl2(1,2)
      by = Bl2(2,2)
      cy = Bl2(3,2)
      dy = 0d0      
      
      RETURN
      END SUBROUTINE l2_project


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      SUBROUTINE quad_interp(i,seg,dt,t,ax,bx,cx,dx,ay,by,cy,dy)
      
      USE globals, ONLY: base
      USE calc_spline, ONLY: newton
      
      IMPLICIT NONE       
      
      INTEGER, INTENT(IN) :: i,seg       
      REAL(rp), INTENT(IN) :: t,dt
      REAL(rp), INTENT(INOUT) :: ax,bx,cx,dx
      REAL(rp), INTENT(INOUT) :: ay,by,cy,dy      
      

      INTEGER :: n1,n2
      REAL(rp) :: tpt
      REAL(rp) :: x,y
      REAL(rp) :: xr(2)
      REAL(rp) :: n1x,n1y,n2x,n2y
      
      REAL(rp) :: Al2(3,3),Bl2(3,2)
      INTEGER :: ipiv(3),info         
      
      
      n1 = base%fbnds(i,seg)
      n2 = base%fbnds(i+1,seg)
                   
      n1x = base%xy(1,n1)
      n1y = base%xy(2,n1)
          
      n2x = base%xy(1,n2)
      n2y = base%xy(2,n2)       
      
      xr(1) = .5d0*n1x + .5d0*n2x
      xr(2) = .5d0*n1y + .5d0*n2y
              
      tpt = .5d0*dt + t          ! initial guess for iteration                                               
      CALL newton(tpt,t,xr,ax,bx,cx,dx,ay,by,cy,dy,x,y) 
      
      Al2(1,1) = 1d0
      Al2(1,2) = 0d0
      Al2(1,3) = 0d0
              
      Al2(2,1) = 1d0
      Al2(2,2) = .5d0*dt
      Al2(2,3) = (.5d0*dt)**2
              
      Al2(3,1) = 1d0             
      Al2(3,2) = dt
      Al2(3,3) = dt**2
              
      Bl2(1,1) = n1x
      Bl2(2,1) = x
      Bl2(3,1) = n2x
                         
      Bl2(1,2) = n1y
      Bl2(2,2) = y
      Bl2(3,2) = n2y              
              
      CALL DGESV(3, 2, Al2, 3, ipiv, Bl2, 3, info)
      
      ax = Bl2(1,1)
      bx = Bl2(2,1)
      cx = Bl2(3,1)
      dx = 0d0
              
      ay = Bl2(1,2)
      by = Bl2(2,2)
      cy = Bl2(3,2)
      dy = 0d0        
      
      
      RETURN
      END SUBROUTINE quad_interp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      SUBROUTINE max_diff(dt,ti,xn1,xn2,ax,bx,cx,dx,ay,by,cy,dy,xd)
      
      USE calc_spline, ONLY: eval_cubic_spline
      
      IMPLICIT NONE      
          
      REAL(rp), INTENT(IN) :: dt,ti
      REAL(rp), INTENT(IN) :: xn1(2),xn2(2)
      REAL(rp), INTENT(IN) :: ax,bx,cx,dx,ay,by,cy,dy
      REAL(rp), INTENT(OUT) :: xd
      
      INTEGER :: it,maxit
      REAL(rp) :: tol,eps        
      REAL(rp) :: t,r,lambda 
      REAL(rp) :: xe,xep,xepp,ye,yep,yepp
      REAL(rp) :: xs,xsp,xspp,ys,ysp,yspp   
      REAL(rp) :: x1,x2,y1,y2
      REAL(rp) :: F,dFdt,dFdr,dFdl
      REAL(rp) :: G,dGdt,dGdr,dGdl
      REAL(rp) :: E,dEdt,dEdr,dEdl
      REAL(rp) :: dDdt,dDdtdt,dDdrdt,dDdr,dDdrdr
      REAL(rp) :: w,dwdt,dwdtdt,dwdrdt,dwdr,dwdrdr
      REAL(rp) :: J(3,3),H(3,1)        
      INTEGER :: ipiv(3),info
      
      tol = 1d-8
      maxit = 1000
      
      t = ti + .5d0*dt
      r = 0d0
      lambda = 0d0
      
      x1 = xn1(1)
      x2 = xn2(1)
      
      y1 = xn1(2)
      y2 = xn2(2)
      
iter: DO it = 1,maxit

        xe = .5d0*(1d0-r)*x1 + .5d0*(1d0+r)*x2
        ye = .5d0*(1d0-r)*y1 + .5d0*(1d0+r)*y2
        
        xep = -.5d0*x1 + .5d0*x2
        yep = -.5d0*y1 + .5d0*y2
        
        xepp = 0d0
        yepp = 0d0

        CALL eval_cubic_spline(t,ti,ax,bx,cx,dx,xs,xsp,xspp)               
        CALL eval_cubic_spline(t,ti,ay,by,cy,dy,ys,ysp,yspp)   
        
        dDdt = 2d0*((xs-xe)*xsp + (ys-ye)*ysp)
        dDdtdt = 2d0*(xsp**2 + xspp*(xs-xe) + ysp**2 + yspp*(ys-ye))
        dDdrdt = -2d0*(xsp*xep+ysp*yep)
        dDdr = -2d0*((xs-xe)*xep + (ys-ye)*yep)
        dDdrdr = -2d0*(xep**2 + yep**2)
        
        w = (x1-xe)*(xs-xe) + (y1-ye)*(ys-ye)
        dwdt = (x1-xe)*xsp + (y1-ye)*ysp
        dwdtdt = (x1-xe)*xspp + (y1-ye)*yspp
        dwdrdt = -xep*xs-yep*ysp
        dwdr = -xep*(x1-2d0*xe+xs) - yep*(y1-2d0*ye+ys)
        dwdrdr = -xepp*(x1-2d0*xe+xs) + 2d0*xep**2 - yepp*(y1-2d0*ye+ys) + 2d0*yep**2
        
        
        
        dFdt = dDdtdt - lambda*dwdtdt
        dFdr = dDdrdt - lambda*dwdrdt
        dFdl = -dwdt
        
        dGdt = dDdrdt - lambda*dwdrdt
        dGdr = dDdrdr - lambda*dwdrdr
        dGdl = -dwdr
        
        dEdt = dwdt
        dEdr = dwdr
        dEdl = 0d0
        
        J(1,1) = dFdt
        J(2,1) = dGdt
        J(3,1) = dEdt
        
        J(1,2) = dFdr
        J(2,2) = dGdr       
        J(3,2) = dEdr
        
        J(1,3) = dFdl
        J(2,3) = dGdl
        J(3,3) = dEdl
        
        F = dDdt - lambda*dwdt
        G = dDdr - lambda*dwdr
        E = w        
        
        H(1,1) = -F
        H(2,1) = -G
        H(3,1) = -E
        
        CALL DGESV(3,1,J,3,ipiv,H,3,info)
        
        t = t + H(1,1)
        r = r + H(2,1)
        lambda = lambda + H(3,1)
        
        eps = max(H(1,1),H(2,1),H(3,1))
        
        IF (ABS(eps) < tol) THEN
!           PRINT*, "iterations", it
!           PRINT*, d
          EXIT iter
        ENDIF
        
        
      ENDDO iter
      
      xe = .5d0*(1d0-r)*x1 + .5d0*(1d0+r)*x2
      ye = .5d0*(1d0-r)*y1 + .5d0*(1d0+r)*y2

      CALL eval_cubic_spline(t,ti,ax,bx,cx,dx,xs)               
      CALL eval_cubic_spline(t,ti,ay,by,cy,dy,ys)         
      
      xd = sqrt((xs-xe)**2 + (ys-ye)**2)
      
      WRITE(90,"(7(E26.17))") xe,ye,xs,ys,t,r,lambda
      
      IF (it >= maxit) THEN
        PRINT "(A,E28.16)", "MAX ITERATIONS EXCEEDED IN FINDING MAX DEFORMATION, ERROR: ", ABS(eps)        
      ENDIF     

      
      RETURN
      END SUBROUTINE max_diff
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      
      END MODULE