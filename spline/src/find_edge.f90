
      SUBROUTINE find_edge(pt1,pt2,xm,el_in,led,base_bou,base_bed,base_led)
      
      USE globals, ONLY: rp,base,eval,nverts
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: pt1,pt2
      REAL(rp), INTENT(IN) :: xm(2)
      INTEGER, INTENT(IN) :: el_in      
      INTEGER, INTENT(IN) :: led(4)
      
      INTEGER, INTENT(OUT) :: base_bed,base_bou,base_led     
      
      INTEGER :: nv,ged,bou,n1    
      INTEGER :: i,j,ed
      INTEGER :: found
      INTEGER :: n1bed,n2bed,n1ed1,n2ed1
      REAL(rp) :: x1(2),x2(2),x3(2),x4(2)
      REAL(rp) :: ax,ay,bx,by,cx,cy,dx,dy
      REAL(rp) :: r,t
      REAL(rp) :: eps

      ! Try to find boundary edge based on lowest difference between 
      ! reference element area and sub-element area sum and
      ! the edge with lowest sub-element area  
      
       nv = nverts(base%el_type(el_in))          
      
      
       found = 0      
       
 edge: DO i = 1,nv
         ged = base%el2ged(el_in,led(i))
        
        
        DO ed = 1,base%nnfbed
          
          IF (ged == base%nfbedn(ed)) THEN
            bou = base%nfbednn(ed,1)
            n1 = base%nfbednn(ed,3)
            
            found = 1
            base_bou = bou
            base_bed = n1
            base_led = led(i)
            
!             PRINT*, base_bou,base_bed,base%nfbednn(ed,2), base%fbnds(n1,bou), base%fbnds(n1+1,bou),el_in
            
            EXIT edge
          ENDIF
        ENDDO
      
      ENDDO edge
      
      
      
      ! Try to find the base boundary edge based on intersection between 
      ! line perpendicular to eval edge and base edge
            
      eps = 1d-12      
      IF (found == 0) THEN
        
        
        x1(1) = eval%xy(1,pt1)
        x1(2) = eval%xy(2,pt1)
        
        x2(1) = eval%xy(1,pt2)
        x2(2) = eval%xy(2,pt2)        
        

   bnd: DO ed = 1,base%nnfbed
              
          ged = base%nfbedn(ed)
          n1bed = base%ged2nn(1,ged)
          n2bed = base%ged2nn(2,ged)
          
          
          x3(1) = base%xy(1,n1bed)
          x3(2) = base%xy(2,n1bed)
          
          x4(1) = base%xy(1,n2bed)
          x4(2) = base%xy(2,n2bed)
          
          ax = xm(1)
          ay = xm(2)
          
          bx = -(xm(2)-x1(2))/(xm(1)-x1(1))
          by = 1d0
          
          cx = .5d0*(x3(1)+x4(1))
          cy = .5d0*(x3(2)+x4(2))
          
          dx = .5d0*(x4(1)-x3(1))
          dy = .5d0*(x4(2)-x3(2))
          
          t = (-dy*(cx-ax) + dx*(cy-ay))/(-bx*dy+by*dx)
          r = (-by*(cx-ax) + bx*(cy-ay))/(-bx*dy+by*dx)               
          
          IF ((r>=-1d0-eps .and. r<=1d0+eps) ) THEN
          
            bou = base%nfbednn(ed,1)
            n1 = base%nfbednn(ed,3)
            
            found = 1
            base_bou = bou
            base_bed = n1 
          
!             PRINT*, "n1bed = ",n1bed, "n2bed = ",n2bed 
!             PRINT*, "R = ", r
!             PRINT*, "T = ", t
!             PRINT*, "X = ", xm(1) + bx*t
!             PRINT*, "Y = ", xm(2) + t
!             
!             PRINT*, "X1 = ", x1(1), "Y1 = ", x1(2)            
!             PRINT*, "XM = ", xm(1), "YM = ", xm(2)
!             PRINT*, "X2 = ", x2(1), "Y2 = ", x2(2)   

            EXIT bnd
            
          ENDIF          
            
              
              
       ENDDO bnd     
       
       IF (found == 0) THEN
         PRINT*, "BOUNDARY EDGE NOT FOUND"
         STOP
       ENDIF
        
        
      ENDIF          
      
      RETURN
      END SUBROUTINE find_edge

      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 