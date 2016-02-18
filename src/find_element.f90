      MODULE find_element

      USE globals, ONLY: rp
      USE kdtree2_module
      USE lapack_interfaces

      CONTAINS
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      

      SUBROUTINE in_element(seg,pt1,pt2,xt,el_found,bed)

      USE globals, ONLY: base,tree_xy,srchdp,closest,fbnds

      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: seg  
      INTEGER, INTENT(IN) :: pt1,pt2        
      REAL(rp), INTENT(IN) :: xt(2)
      
      INTEGER, INTENT(OUT) :: bed
      INTEGER, INTENT(OUT) :: el_found
            
      INTEGER :: srch,i
      INTEGER :: el,eln,clnd
      INTEGER :: found,ed_found      
      INTEGER :: min_el
      REAL(rp) :: diff,min_diff
      REAL(rp) :: tol
            
      
      tol = 1d-5 
        CALL kdtree2_n_nearest(tp=tree_xy,qv=xt,nn=srchdp,results=closest) ! find what element xt is in               
        
        ! Test elements to see which element point is located in    
        found = 0      
        min_diff = 999d0
search: DO srch = 1,srchdp

          clnd = fbnds(closest(srch)%idx)
          
!           PRINT("(A,I5)"), "   closest node: ", clnd

   elem:  DO el = 1,base%nepn(clnd) 
 
            eln = base%epn(el,clnd)
            
!             PRINT("(A,I5)"), "   testing: ", eln                               

            ! Compute sum of sub-triangle areas
            CALL sub_element(pt1,eln,xt,diff,ed_found)            
            
            IF (diff < min_diff) THEN  ! keep track of element with minimum difference in sum of sub-triangle areas                                  
              min_diff = diff          ! to return if tolerance is not met
              min_el = eln
            ENDIF
          
            ! The station is in the element if the reference element area and sum of sub triangle are the same
            ! Make sure vertexes use elements with boundary edges
            IF (diff < tol .AND. base%bel_flag(eln) == 1) THEN
              PRINT("(A,I5)"), "   element found", eln                            
                      
              el_found = eln        
              found = 1                        
              
              ! find base edge (to get correct spline coefficients)                 
              CALL find_edge(pt1,pt2,xt,seg,el_found,ed_found,bed)
                          
              EXIT search            
            ENDIF    
            
!             PRINT*, " "  
          ENDDO elem
        
        ENDDO search    
        

        
        IF (found == 0) THEN
          el_found = min_el        

          PRINT*, "ELEMENT NOT FOUND FOR POINT: ",pt1  
          PRINT*, "USING ELEMENT ", el_found, "(AREA DIFF = ",min_diff, ")"       

          CALL sub_element(pt1,el_found,xt,diff,ed_found)            
          CALL find_edge(pt1,pt2,xt,seg,el_found,ed_found,bed)
          
        ENDIF     
 

      RETURN
      END SUBROUTINE in_element
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE find_edge(pt1,pt2,xm,seg,el_in,led,base_bed)
      
      USE globals, ONLY: base,eval,nverts
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: pt1,pt2
      REAL(rp), INTENT(IN) :: xm(2)
      INTEGER, INTENT(IN) :: seg    
      INTEGER, INTENT(IN) :: el_in      
      INTEGER, INTENT(IN) :: led
      
      INTEGER, INTENT(OUT) :: base_bed     
      
      INTEGER :: nvert      
      INTEGER :: j
      INTEGER :: found
      INTEGER :: n1bed,n2bed,n1ed1,n2ed1
      REAL(rp) :: x1(2),x2(2),x3(2),x4(2)
      REAL(rp) :: ax,ay,bx,by,cx,cy,dx,dy
      REAL(rp) :: r,t


      
       nvert = nverts(base%el_type(el_in))            
            
       found = 0
              
      ! Try to find boundary edge based on lowest difference between 
      ! reference element area and sub-element area sum and
      ! the edge with lowest sub-element area              
 bseg: DO j = 1,base%fbseg(1,seg)-1
              
         n1bed = base%fbnds(j,seg)
         n2bed = base%fbnds(j+1,seg)             
 
         n1ed1 = base%vct(mod(led+0,nvert)+1,el_in)
         n2ed1 = base%vct(mod(led+1,nvert)+1,el_in)           
                                                        
         IF(((n1ed1 == n1bed).AND.(n2ed1 == n2bed)).OR. &
            ((n1ed1 == n2bed).AND.(n2ed1 == n1bed))) THEN
            PRINT*, "n1bed = ",n1bed, "n2bed = ",n2bed    
!           PRINT*, n1bed, base%xy(1,n1bed), base%xy(2,n1bed)

            found = 1                   
            base_bed = j
              
            EXIT bseg
         ENDIF        
              
              
      ENDDO bseg
      
      ! Try to find the base boundary edge based on intersection between 
      ! line perpendicular to eval edge and base edge
            
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
          
          IF ((r>=-1d0 .and. r<=1d0) ) THEN
          
            found = 1                   
            base_bed = j          
          
            PRINT*, "n1bed = ",n1bed, "n2bed = ",n2bed 
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
        
        
      ENDIF          
      
      RETURN
      END SUBROUTINE find_edge

      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      SUBROUTINE sub_element(pt,eln,xt,diff,closest_ed)
      
      USE globals, ONLY: base,nverts,rsre           
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: pt
      INTEGER, INTENT(IN) :: eln
      REAL(rp), INTENT(IN) :: xt(2)
      
      INTEGER, INTENT(OUT) :: closest_ed
      REAL(rp), INTENT(OUT) :: diff
      
      INTEGER :: i    
      INTEGER :: et,nvert
      INTEGER :: n1,n2
      REAL(rp) :: area,area_sum
      REAL(rp) :: stri_area,stri_min
      REAL(rp) :: dist
      REAL(rp) :: r(2)
      REAL(rp) :: x(3),y(3)
      
      et = base%el_type(eln)
      nvert = nverts(et)                                           
          
      ! Compute the local (r,s) coordinates of the (x,y) station location
      CALL newton(pt,xt,eln,r)
          
      ! Find reference element area
      IF (mod(et,2) == 1) THEN
        area = 2d0
      ELSE IF (mod(et,2) == 0) THEN
        area = 4d0
      ENDIF          
          
      ! Compute sum of sub-triangle areas
      area_sum = 0d0
      stri_min = 999d0
      DO i = 1,nvert
        n1 = mod(i+0,nvert)+1
        n2 = mod(i+1,nvert)+1           

        x(1) = rsre(1,n1,et)
        y(1) = rsre(2,n1,et)
            
        x(2) = rsre(1,n2,et)
        y(2) = rsre(2,n2,et)
            
        x(3) = r(1)
        y(3) = r(2)
              
        stri_area = .5d0*abs((x(2)-x(1))*(y(3)-y(1)) - (x(3)-x(1))*(y(2)-y(1)))
        
!         PRINT "(A,I7,A,F14.7,A,I7,A,I7)", "stri: ",i," area: ",stri_area, " n1: ", base%vct(n1,eln), " n2: ",base%vct(n2,eln)             
              
        IF (stri_area < stri_min) THEN ! keep track of minimum sub-triangle area to determine
          stri_min = stri_area         ! which edge the point lies on, or is closest to
          closest_ed = i
        ENDIF
            
        area_sum = area_sum + stri_area
      ENDDO
          
!     PRINT("(A,F20.15,A,F20.15)"), "   area = ",area, "   area sum = ", area_sum
!     PRINT*, " "
          
      diff = abs(area-area_sum) 
      
      
      DO i = 1,nvert
        
        x(1) = rsre(1,i,et)
        y(1) = rsre(2,i,et)
            
        x(3) = r(1)
        y(3) = r(2)
              
!         dist = sqrt((x(1)-x(3))**2 + (y(1)-y(3))**2)        
!         PRINT "(A,I7,A,F14.7,A,I7,A,I7)", "node: ",i," dist: ",dist, " n1: ", base%vct(i,eln)     
              
        IF (stri_area < stri_min) THEN ! keep track of minimum sub-triangle area to determine
          stri_min = stri_area         ! which edge the point lies on, or is closest to
          closest_ed = i
        ENDIF
            
        area_sum = area_sum + stri_area
      ENDDO      
     
      
      RETURN
      END SUBROUTINE sub_element     

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
      
      SUBROUTINE newton(pt,x,eln,r)

      USE globals, ONLY: rp,np,nnds,mnnds,V,base,ipiv
      USE basis, ONLY: tri_basis,quad_basis

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: pt
      INTEGER, INTENT(IN) :: eln
      REAL(rp), INTENT(IN) :: x(2)      
      REAL(rp), INTENT(OUT) :: r(2)      
      
      INTEGER :: it,et,p,n,i
      INTEGER :: info
      INTEGER :: maxit
      REAL(rp) :: tol
      REAL(rp) :: f,g,error
      REAL(rp) :: dfdr,dfds,dgdr,dgds,jac
      REAL(rp) :: phi(mnnds,1),dpdr(mnnds,1),dpds(mnnds,1)
      REAL(rp) :: l(mnnds,3)
        
      tol = 1d-9
      maxit = 100
      info = 0
      
      et = base%el_type(eln)
      p = np(et)  
        
      IF (mod(et,2) == 1) THEN
        r(1) = -1d0/3d0
        r(2) = -1d0/3d0
      ELSE IF (mod(et,2) == 0) THEN
        r(1) = 0d0
        r(2) = 0d0
      ENDIF

      DO it = 1,maxit     
           
        IF (mod(et,2) == 1) THEN
          CALL tri_basis(p,n,1,r(1),r(2),phi,dpdr,dpds)
        ELSE IF (mod(et,2) == 0) THEN
          CALL quad_basis(p,n,1,r(1),r(2),phi,dpdr,dpds)
        ENDIF
        
        DO i = 1,n

          l(i,1) = phi(i,1)
          l(i,2) = dpdr(i,1)
          l(i,3) = dpds(i,1)                
       
        ENDDO     
          

        CALL DGETRS("N",n,3,V(1,1,et),mnnds,ipiv(1,et),l,mnnds,info)
!         IF (info /= 0 ) PRINT*, "LAPACK ERROR"      
        
        dfdr = 0d0
        dfds = 0d0
        dgdr = 0d0
        dgds = 0d0
        f = 0d0
        g = 0d0        
        
        DO i = 1,n

          dfdr = dfdr + l(i,2)*base%elxy(i,eln,1)
          dfds = dfds + l(i,3)*base%elxy(i,eln,1)
          dgdr = dgdr + l(i,2)*base%elxy(i,eln,2)
          dgds = dgds + l(i,3)*base%elxy(i,eln,2)
          
          f = f + l(i,1)*base%elxy(i,eln,1)
          g = g + l(i,1)*base%elxy(i,eln,2)
        ENDDO
        
        jac = dfdr*dgds - dgdr*dfds
        
        f = f - x(1)
        g = g - x(2)
        
        r(1) = r(1) - (1d0/jac)*( dgds*f - dfds*g)
        r(2) = r(2) - (1d0/jac)*(-dgdr*f + dfdr*g)   
!         PRINT("(3(A,F20.15))"), "   f = ",f, "   g = ", g, "  jac = ", jac        
!         PRINT("(2(A,F20.15))"), "   r = ",r(1), "   s = ", r(2)
              
        
        IF (ABS(f) < tol .AND. ABS(g) < tol) THEN
          EXIT
        ENDIF        
       
        
      ENDDO
      
      error = max(abs(f),abs(g))
!       IF (it >= maxit) THEN
!         PRINT("(A,E22.15)"), "   MAX ITERATIONS EXCEEDED, error = ",error
!         PRINT("(2(A,F20.15))"), "   r = ",r(1), "   s = ", r(2)
!       ELSE       
!         PRINT("(A,I7,A,E22.15)"), "   iterations: ",it, "  error = ",error
!         PRINT("(2(A,F20.15))"), "   r = ",r(1), "   s = ", r(2)
!       ENDIF


      RETURN
      END SUBROUTINE newton
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      


      END MODULE find_element