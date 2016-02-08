      MODULE check

      USE globals, ONLY:rp,pi

      IMPLICIT NONE

      CONTAINS
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE check_angle(seg,n,nd,theta1,theta2,edlen)
      
      USE globals, ONLY: base
      
      IMPLICIT NONE     
      
      INTEGER, INTENT(IN) :: nd,n,seg
      
      INTEGER :: n1bed,n2bed,n3bed,n4bed      
      REAL(rp) :: n1x,n1y,n2x,n2y,n3x,n3y,n4x,n4y
!       REAL(rp) :: angle
      REAL(rp), INTENT(OUT) :: theta1,theta2   
      REAL(rp), INTENT(OUT) :: edlen

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
            
      END SUBROUTINE check_angle
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE check_deformation(minedlen,xm,ym,x,y)
      
      IMPLICIT NONE
      
      REAL(rp), INTENT(IN) :: minedlen
      REAL(rp), INTENT(IN) :: xm,ym
      REAL(rp), INTENT(INOUT) :: x,y
      REAL(rp) :: xs,ys
      REAL(rp) :: r,d1
      
      xs = x
      ys = y
      
      d1 = sqrt((xm-x)**2 + (ym-y)**2)
         
      r = -1d0
      DO WHILE (d1 > .1d0*minedlen .AND. r < 1d0)
        r = r + .1d0
                       
        x = .5d0*(1d0-r)*xs + .5d0*(1d0+r)*xm
        y = .5d0*(1d0-r)*ys + .5d0*(1d0+r)*ym
                       
        d1 = sqrt((xm-x)**2 + (ym-y)**2)
      ENDDO

      IF (d1 <= .1d0*minedlen) THEN
        RETURN
      ELSE
        PRINT*, "DEFORMATION TOO LARGE"
        x = xs
        y = ys
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

      SUBROUTINE quad_interp(i,seg,t,dt,ax,bx,cx,dx,ay,by,cy,dy)
      
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
      
      END MODULE