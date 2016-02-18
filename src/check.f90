      MODULE check

      USE globals, ONLY:rp,pi

      IMPLICIT NONE

      CONTAINS
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE check_angle(seg,nd,dt,ti,ax,bx,cx,dx,ay,by,cy,dy)
      
      USE globals, ONLY: base,theta_tol
      USE calc_spline, ONLY: newton      
      
      IMPLICIT NONE     
      
      INTEGER, INTENT(IN) :: nd,seg
      REAL(rp), INTENT(IN) :: dt,ti
      REAL(rp), INTENT(INOUT) :: ax,bx,cx,dx,ay,by,cy,dy      
      
      INTEGER :: n
      INTEGER :: n1bed,n2bed,n3bed,n4bed      
      REAL(rp) :: n1x,n1y,n2x,n2y,n3x,n3y,n4x,n4y
      REAL(rp) :: theta1,theta2   
      REAL(rp) :: edlen
      REAL(rp) :: x,y
      REAL(rp) :: xr(2)
      REAL(rp) :: xi(3),yi(3)
      REAL(rp) :: r

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

        xr(1) = .5d0*n1x + .5d0*n2x
        xr(2) = .5d0*n1y + .5d0*n2y
              
        r = 0d0                                        
        CALL newton(r,dt,ti,xr,ax,bx,cx,dx,ay,by,cy,dy,x,y) 
        
        xi(1) = n1x
        xi(2) = x
        xi(3) = n2x
        
        yi(1) = n1y
        yi(2) = y
        yi(3) = n2y
        
        r = 0d0
        CALL quad_interp(r,dt,ti,xi,yi,ax,bx,cx,dx,ay,by,cy,dy)
              
      ENDIF                
            
      END SUBROUTINE check_angle
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE check_deformation(seg,nd,dt,ti,ax,bx,cx,dx,ay,by,cy,dy)
      
      USE globals, ONLY: base,deform_tol
      USE find_element, ONLY: in_element
      USE calc_spline, ONLY: eval_cubic_spline
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: seg,nd
      REAL(rp), INTENT(IN) :: dt,ti
      REAL(rp), INTENT(INOUT) :: ax,bx,cx,dx,ay,by,cy,dy
      INTEGER :: i
      INTEGER :: n1,n2
      REAL(rp) :: xn1(2),xn2(2)
      REAL(rp) :: xa(2)
      REAL(rp) :: xe,ye
      REAL(rp) :: x,y
      REAL(rp) :: xs,ys
      REAL(rp) :: xsp,ysp
      REAL(rp) :: r,t,lam
      REAL(rp) :: xd,max_dist
      REAL(rp) :: ri(2),xi(4),yi(4)
      INTEGER :: el_in,bed 
        

        n1 = base%fbnds(nd,seg)
        n2 = base%fbnds(nd+1,seg) 

        xn1(1) = base%xy(1,n1)
        xn1(2) = base%xy(2,n1)
        
        xn2(1) = base%xy(1,n2)
        xn2(2) = base%xy(2,n2)
        
        xa(1) = .5d0*(xn1(1)+xn2(1))
        xa(2) = .5d0*(xn1(2)+xn2(2))
      
        CALL in_element(seg,n1,n2,xa,el_in,bed)  
        max_dist = deform_tol*base%minedlen(el_in)
        
        r = 0d0 
        t = .5d0*dt + ti
        lam = 0d0
        CALL max_diff(dt,ti,xn1,xn2,ax,bx,cx,dx,ay,by,cy,dy,xd,xs,ys,xe,ye,t,r,lam)
        
        
        IF (r < -1d0 .or. r > 1d0) THEN
          IF (r < -1d0) THEN      ! Try new intial guesses if the max deformation is found outside        
            r = .75d0             ! the (-1, 1) r interval.  This seems to work for now, but something      
          ELSE IF (r > 1d0) THEN  ! more sophisticated might be needed later on.    
            r = -.75d0               
          ENDIF
        
          t = .5d0*dt*(r+1d0) + ti
          lam = 0d0
          CALL max_diff(dt,ti,xn1,xn2,ax,bx,cx,dx,ay,by,cy,dy,xd,xs,ys,xe,ye,t,r,lam)  
        ENDIF
        
        IF (r < -1d0 .or. r > 1d0) THEN
          PRINT*, "R VALUE OUT OF RANGE: ", r 
          PAUSE
        ENDIF
        
        
        
        WRITE(90,"(7(E26.17))") xe,ye,xs,ys,t,r,lam
        
        
        
        IF (xd > max_dist) THEN         
          PRINT*, "DEFORMATION TOO LARGE"
          
          xi(1) = xn1(1)
          xi(4) = xn2(1)
          
          yi(1) = xn1(2)
          yi(4) = xn2(2)          
          
          IF ( r > -1d0/3d0 .and. r < 1d0/3d0) THEN ! If max is in the middle of the interval, match the limited 
            ri(1) = r - .1d0                        ! value a little to the left and right of the detected max
            ri(2) = r + .1d0
            
            DO i = 1,2
          
              xe = .5d0*(1d0-ri(i))*xn1(1) + .5d0*(1d0+ri(i))*xn2(1)
              ye = .5d0*(1d0-ri(i))*xn1(2) + .5d0*(1d0+ri(i))*xn2(2)
            
              x = xs ! initial guess
              y = ys
              CALL set_dist(xn1(1),xn1(2),xe,ye,max_dist,x,y)
            
              xi(i+1) = x
              yi(i+1) = y
              
            ENDDO         
            
            CALL cubic_interp(0,ri,dt,ti,xi,yi,ax,bx,cx,dx,ay,by,cy,dy)            
          ELSE                                      ! If max is to the left or right of the middle of the interval,          
                                                    ! match the limited value at the detected max and the dervative
            ri(1) = r                               ! value at the opposite end
                                           
            xe = .5d0*(1d0-ri(1))*xn1(1) + .5d0*(1d0+ri(1))*xn2(1)
            ye = .5d0*(1d0-ri(1))*xn1(2) + .5d0*(1d0+ri(1))*xn2(2)            
            
            x = xs ! initial guess
            y = ys
            CALL set_dist(xn1(1),xn1(2),xe,ye,max_dist,x,y)            
            
            xi(2) = x
            yi(2) = y
            
            IF (r <= -1d0/3d0) THEN
              ri(2) = 1d0
            ENDIF 
            
            IF (r >= 1d0/3d0) THEN
              ri(2) = -1d0
            ENDIF                        
            
            t = .5d0*dt*(ri(2)+1d0) + ti
            CALL eval_cubic_spline(t,ti,ax,bx,cx,dx,xs,xsp)               
            CALL eval_cubic_spline(t,ti,ay,by,cy,dy,ys,ysp)
            
            xi(3) = xsp
            yi(3) = ysp
          
            CALL cubic_interp(1,ri,dt,ti,xi,yi,ax,bx,cx,dx,ay,by,cy,dy)             
          
          ENDIF          
          
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

      SUBROUTINE quad_interp(r,dt,ti,x,y,ax,bx,cx,dx,ay,by,cy,dy)      
      
      IMPLICIT NONE       
      
      REAL(rp), INTENT(IN) :: r
      REAL(rp), INTENT(IN) :: dt,ti
      REAL(rp), INTENT(IN) :: x(3),y(3)      
      REAL(rp), INTENT(INOUT) :: ax,bx,cx,dx
      REAL(rp), INTENT(INOUT) :: ay,by,cy,dy           
      
      REAL(rp) :: t
      REAL(rp) :: A(3,3),B(3,2)
      INTEGER :: ipiv(3),info         
      
      
      t = .5d0*dt*(r+1d0) + ti

      
      A(1,1) = 1d0
      A(1,2) = 0d0
      A(1,3) = 0d0
              
      A(2,1) = 1d0
      A(2,2) = (t-ti)
      A(2,3) = (t-ti)**2
              
      A(3,1) = 1d0             
      A(3,2) = dt
      A(3,3) = dt**2
              
      B(1,1) = x(1)
      B(2,1) = x(2)
      B(3,1) = x(3)
                         
      B(1,2) = y(1)
      B(2,2) = y(2)
      B(3,2) = y(3)             
              
      CALL DGESV(3, 2, A, 3, ipiv, B, 3, info)
      
      ax = B(1,1)
      bx = B(2,1)
      cx = B(3,1)
      dx = 0d0
              
      ay = B(1,2)
      by = B(2,2)
      cy = B(3,2)
      dy = 0d0        
      
      
      RETURN
      END SUBROUTINE quad_interp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      SUBROUTINE cubic_interp(deriv,r,dt,ti,x,y,ax,bx,cx,dx,ay,by,cy,dy)      
      
      IMPLICIT NONE       
          
      INTEGER, INTENT(IN) :: deriv    
      REAL(rp), INTENT(IN) :: r(2),dt,ti
      REAL(rp), INTENT(IN) :: x(4),y(4)      
      REAL(rp), INTENT(INOUT) :: ax,bx,cx,dx
      REAL(rp), INTENT(INOUT) :: ay,by,cy,dy      
      
      REAL(rp) :: t(2)      
      REAL(rp) :: A(4,4),B(4,2)
      INTEGER :: ipiv(4),info                    
      
      t(1) = .5d0*dt*(r(1)+1d0)+ti    ! intermediate boundary edge interpolation points
      t(2) = .5d0*dt*(r(2)+1d0)+ti
      
      A(1,1) = 1d0      ! iterpolate at beginning boundary node (t = ti)
      A(1,2) = 0d0
      A(1,3) = 0d0
      A(1,4) = 0d0      
              
      A(2,1) = 1d0
      A(2,2) = (t(1)-ti)
      A(2,3) = (t(1)-ti)**2
      A(2,4) = (t(1)-ti)**3     
      
      IF (deriv == 0) THEN
        A(3,1) = 1d0             
        A(3,2) = (t(2)-ti)
        A(3,3) = (t(2)-ti)**2
        A(3,4) = (t(2)-ti)**3         
      ELSE
        A(3,1) = 0d0             
        A(3,2) = 1d0
        A(3,3) = 2d0*(t(2)-ti)
        A(3,4) = 3d0*(t(2)-ti)**2          
      ENDIF
              
      A(4,1) = 1d0      ! iterpolate at end boundary node (t = ti + dt)       
      A(4,2) = dt
      A(4,3) = dt**2
      A(4,4) = dt**3    
              
      B(1,1) = x(1)
      B(2,1) = x(2)
      B(3,1) = x(3)
      B(4,1) = x(4)      
                         
      B(1,2) = y(1)
      B(2,2) = y(2)
      B(3,2) = y(3)    
      B(4,2) = y(4)        
              
      CALL DGESV(4, 2, A, 4, ipiv, B, 4, info)
      
      ax = B(1,1)
      bx = B(2,1)
      cx = B(3,1)
      dx = B(4,1)
              
      ay = B(1,2)
      by = B(2,2)
      cy = B(3,2)
      dy = B(4,2)        
      
      
      RETURN
      END SUBROUTINE cubic_interp      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      SUBROUTINE max_diff(dt,ti,xn1,xn2,ax,bx,cx,dx,ay,by,cy,dy,xd,xs,ys,xe,ye,t,r,lambda)
      
      USE calc_spline, ONLY: eval_cubic_spline
      
      IMPLICIT NONE      
          
      REAL(rp), INTENT(IN) :: dt,ti
      REAL(rp), INTENT(IN) :: xn1(2),xn2(2)
      REAL(rp), INTENT(IN) :: ax,bx,cx,dx,ay,by,cy,dy
      REAL(rp), INTENT(OUT) :: xd,xs,ys,xe,ye
      REAL(rp), INTENT(INOUT) :: t,r,lambda
      
      INTEGER :: it,maxit
      REAL(rp) :: tol,eps        
      REAL(rp) :: xep,xepp,yep,yepp
      REAL(rp) :: xsp,xspp,ysp,yspp   
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
      
!       t = ti + .5d0*dt
!       r = 0d0
!       lambda = 0d0
      
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
        
        eps = max(abs(H(1,1)),abs(H(2,1)),abs(H(3,1)))
        
        IF (eps < tol) THEN
          EXIT iter
        ENDIF
        
        
      ENDDO iter
      
      xe = .5d0*(1d0-r)*x1 + .5d0*(1d0+r)*x2
      ye = .5d0*(1d0-r)*y1 + .5d0*(1d0+r)*y2

      CALL eval_cubic_spline(t,ti,ax,bx,cx,dx,xs)               
      CALL eval_cubic_spline(t,ti,ay,by,cy,dy,ys)         
      
      xd = sqrt((xs-xe)**2 + (ys-ye)**2)     
      
      IF (it >= maxit) THEN
        PRINT "(A,E28.16)", "MAX ITERATIONS EXCEEDED IN FINDING MAX DEFORMATION, ERROR: ", ABS(eps)        
      ENDIF     

      
      RETURN
      END SUBROUTINE max_diff
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE set_dist(x1,y1,xe,ye,dist,x,y)
      
      IMPLICIT NONE
            
      REAL(rp), INTENT(IN) :: x1,y1      
      REAL(rp), INTENT(IN) :: xe,ye
      REAL(rp), INTENT(IN) :: dist
      REAL(rp), INTENT(INOUT) :: x,y
      
            
      INTEGER :: it,maxit
      REAL(rp) :: tol,eps
      
      REAL(rp) :: D,dDdx,dDdy
      REAL(rp) :: w,dwdx,dwdy
      REAL(rp) :: J
      
      
      tol = 1d-8
      maxit = 1000      
      
iter: DO it = 1,maxit
      
        D = (x-xe)**2 + (y-ye)**2 - dist**2
        w = (x1-xe)*(x-xe) + (y1-ye)*(y-ye)
      
        dDdx = 2d0*(x-xe)
        dDdy = 2d0*(y-ye)
        
        dwdx = x1-xe
        dwdy = y1-ye
        
        J = dDdx*dwdy - dDdy*dwdx
      
        x = x - ( dwdy*D - dDdy*w)/J
        y = y - (-dwdx*D + dDdx*w)/J 
        
        eps = max(abs(D),abs(w))
        
        IF (eps < tol) THEN
          EXIT iter
        ENDIF           
      
      ENDDO iter
      
      IF (it >= maxit) THEN
        PRINT "(A,E28.16)", "MAX ITERATIONS EXCEEDED IN FINDING SET DEFORMATION, ERROR: ", ABS(eps)        
      ENDIF          
      
      RETURN
      END SUBROUTINE set_dist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      
      END MODULE