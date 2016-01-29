      MODULE check

      USE globals, ONLY:rp,pi

      IMPLICIT NONE

      CONTAINS
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE check_angle(seg,n,nd,theta1,theta2,edlen)
      
      USE globals, ONLY: fbnds,xy
      
      IMPLICIT NONE     
      
      INTEGER, INTENT(IN) :: nd,n,seg
      
      INTEGER :: n1bed,n2bed,n3bed,n4bed      
      REAL(rp) :: n1x,n1y,n2x,n2y,n3x,n3y,n4x,n4y
!       REAL(rp) :: angle
      REAL(rp), INTENT(OUT) :: theta1,theta2   
      REAL(rp), INTENT(OUT) :: edlen

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
      
      END MODULE