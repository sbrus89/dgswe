
      SUBROUTINE find_edge(seg,pt1,pt2,xm,el_in,led,base_bed)
      
      USE globals, ONLY: rp,base,eval,nverts
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: pt1,pt2
      REAL(rp), INTENT(IN) :: xm(2)
      INTEGER, INTENT(IN) :: seg    
      INTEGER, INTENT(IN) :: el_in      
      INTEGER, INTENT(IN) :: led(4)
      
      INTEGER, INTENT(OUT) :: base_bed     
      
      INTEGER :: nv,ged    
      INTEGER :: i,j
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

         n1ed1 = base%ged2nn(1,ged)
         n2ed1 = base%ged2nn(2,ged)
                
                         
   bseg: DO j = 1,base%fbseg(1,seg)-1
              
           n1bed = base%fbnds(j,seg)
           n2bed = base%fbnds(j+1,seg)                    
                                                        
           IF(((n1ed1 == n1bed).AND.(n2ed1 == n2bed)).OR. &
              ((n1ed1 == n2bed).AND.(n2ed1 == n1bed))) THEN
              PRINT*, "n1bed = ",n1bed, "n2bed = ",n2bed    
!             PRINT*, n1bed, base%xy(1,n1bed), base%xy(2,n1bed)

              found = 1                   
              base_bed = j
                
              EXIT edge
           ENDIF         
              
        ENDDO bseg
      
      ENDDO edge
      
      
      
      ! Try to find the base boundary edge based on intersection between 
      ! line perpendicular to eval edge and base edge
            
      eps = 1d-12      
      IF (found == 0) THEN
        
        
        x1(1) = eval%xy(1,pt1)
        x1(2) = eval%xy(2,pt1)
        
        x2(1) = eval%xy(1,pt2)
        x2(2) = eval%xy(2,pt2)        
        
 bseg2: DO j = 1,base%fbseg(1,seg)-1
              
          n1bed = base%fbnds(j,seg)
          n2bed = base%fbnds(j+1,seg)             
          
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
          
            found = 1                   
            base_bed = j          
          
!             PRINT*, "n1bed = ",n1bed, "n2bed = ",n2bed 
!             PRINT*, "R = ", r
!             PRINT*, "T = ", t
!             PRINT*, "X = ", xm(1) + bx*t
!             PRINT*, "Y = ", xm(2) + t
!             
!             PRINT*, "X1 = ", x1(1), "Y1 = ", x1(2)            
!             PRINT*, "XM = ", xm(1), "YM = ", xm(2)
!             PRINT*, "X2 = ", x2(1), "Y2 = ", x2(2)   

            EXIT bseg2
            
          ENDIF          
            
              
              
       ENDDO bseg2     
       
       IF (found == 0) THEN
         PRINT*, "BOUNDARY EDGE NOT FOUND"
         STOP
       ENDIF
        
        
      ENDIF          
      
      RETURN
      END SUBROUTINE find_edge

      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 