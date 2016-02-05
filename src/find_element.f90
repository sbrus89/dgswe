      MODULE find_element

      USE globals, ONLY: rp
      USE kdtree2_module

      CONTAINS
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      

      SUBROUTINE in_element(seg,dt,pt,xt,t,bed)

      USE globals, ONLY: base,tree_xy,srchdp,closest,nverts,rsre,fbnds

      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: seg  
      INTEGER, INTENT(IN) :: pt
      REAL(rp), INTENT(IN) :: dt          
      REAL(rp), INTENT(IN) :: xt(2)
      
      INTEGER, INTENT(OUT) :: bed
      REAL(rp), INTENT(OUT) :: t
            
      INTEGER :: srch,i
      INTEGER :: el,eln,et,nvert,clnd
      INTEGER :: found,el_found,ed_found      
      INTEGER :: n1,n2
      INTEGER :: min_el
      REAL(rp) :: x(3),y(3),r(1),s(1)
      REAL(rp) :: sarea,area,diff,min_diff
      REAL(rp) :: tol
      
      PRINT*, "FINDING ELEMENT FOR POINT NODE: ",pt      
      
      tol = 1d-5 
        CALL kdtree2_n_nearest(tp=tree_xy,qv=xt,nn=srchdp,results=closest) ! find what element xt is in               
        
        ! Test elements to see which element point is located in    
        found = 0      
        min_diff = 999d0
search: DO srch = 1,srchdp

          clnd = fbnds(closest(srch)%idx)

   elem:  DO el = 1,base%nepn(clnd) 
 
            eln = base%epn(el,clnd)

            et = base%el_type(eln)
            nvert = nverts(et)   
            
!             PRINT("(A,I5)"), "   testing: ", eln
!             PRINT*, " "            
                         
          
            ! Compute the local (r,s) coordinates of the (x,y) station location
            CALL newton(pt,xt(1),xt(2),eln,r,s)
          
            ! Find reference element area
            IF (mod(et,2) == 1) THEN
              area = 2d0
            ELSE IF (mod(et,2) == 0) THEN
              area = 4d0
            ENDIF          
          
            ! Compute sum of sub-triangle areas
            sarea = 0d0
            DO i = 1,nvert
              n1 = mod(i+0,nvert)+1
              n2 = mod(i+1,nvert)+1           

              x(1) = rsre(1,n1,et)
              y(1) = rsre(2,n1,et)
            
              x(2) = rsre(1,n2,et)
              y(2) = rsre(2,n2,et)
            
              x(3) = r(1)
              y(3) = s(1)
            
              sarea = sarea + .5d0*abs((x(2)-x(1))*(y(3)-y(1)) - (x(3)-x(1))*(y(2)-y(1)))
            ENDDO
          
!               PRINT("(A,F20.15,A,F20.15)"), "   area = ",area, "   sarea = ", sarea
!               PRINT*, " "
          
            diff = abs(area-sarea)
            
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
              ed_found = 0               
              CALL find_edge(seg,eln,dt,ed_found,t,bed)
            
              IF (ed_found == 1) THEN
                EXIT search            
              ENDIF
            ENDIF    
            
          ENDDO elem
        
        ENDDO search    
        

        
        IF (found == 0) THEN
          eln = min_el        

          PRINT*, "ELEMENT NOT FOUND FOR POINT: ",pt  
          PRINT*, "USING ELEMENT ", eln, "(AREA DIFF = ",min_diff, ")"       
          
          CALL find_edge(seg,eln,dt,ed_found,t,bed)
          
        ENDIF     
 

      RETURN
      END SUBROUTINE in_element
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE find_edge(seg,el_in,dt,found,t,base_bed)
      
      USE globals, ONLY: base,nverts
      
      IMPLICIT NONE
      
     
      INTEGER, INTENT(IN) :: seg 
      REAL(rp), INTENT(IN) :: dt      
      INTEGER, INTENT(IN) :: el_in      

      INTEGER, INTENT(OUT) :: found
      INTEGER, INTENT(OUT) :: base_bed
      REAL(rp), INTENT(OUT) :: t      
      
      INTEGER :: nvert      
      INTEGER :: led,j
      INTEGER :: n1bed,n2bed,n1ed1,n2ed1


      
       nvert = nverts(base%el_type(el_in))            
            
       found = 0
              
 bseg: DO j = 1,base%fbseg(1,seg)-1
              
         n1bed = base%fbnds(j,seg)
         n2bed = base%fbnds(j+1,seg)   
          
         DO led = 1,nvert
 
           n1ed1 = base%vct(mod(led+0,nvert)+1,el_in)
           n2ed1 = base%vct(mod(led+1,nvert)+1,el_in)           
                                                        
           IF(((n1ed1 == n1bed).AND.(n2ed1 == n2bed)).OR. &
              ((n1ed1 == n2bed).AND.(n2ed1 == n1bed))) THEN
              PRINT*, "n1bed = ",n1bed, "n2bed = ",n2bed    
!             PRINT*, n1bed, base%xy(1,n1bed), base%xy(2,n1bed)

              found = 1
                   
              t = real(j-1,rp)*dt
              base_bed = j

              EXIT bseg
           ENDIF        
              
        ENDDO
              
      ENDDO bseg
            
      IF (found == 0) THEN
        PRINT*, "Boundary edge not found"
!         PAUSE
      ENDIF          
      
      RETURN
      END SUBROUTINE find_edge

      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
      SUBROUTINE newton(pt,x,y,eln,r,s)

      USE globals, ONLY: rp,np,nnds,mnnds,V,base,ipiv
      USE basis, ONLY: tri_basis,quad_basis

      IMPLICIT NONE
      INTEGER :: it,eln,et,p,n,i,pt
      INTEGER :: info
      INTEGER :: maxit
      REAL(rp) :: tol
      REAL(rp) :: x,y
      REAL(rp) :: r(1),s(1)
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
        s(1) = -1d0/3d0
      ELSE IF (mod(et,2) == 0) THEN
        r(1) = 0d0
        s(1) = 0d0
      ENDIF

      DO it = 1,maxit     
           
        IF (mod(et,2) == 1) THEN
          CALL tri_basis(p,n,1,r,s,phi,dpdr,dpds)
        ELSE IF (mod(et,2) == 0) THEN
          CALL quad_basis(p,n,1,r,s,phi,dpdr,dpds)
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
        
        f = f - x
        g = g - y
        
        r(1) = r(1) - (1d0/jac)*( dgds*f - dfds*g)
        s(1) = s(1) - (1d0/jac)*(-dgdr*f + dfdr*g)   
!         PRINT("(3(A,F20.15))"), "   f = ",f, "   g = ", g, "  jac = ", jac        
!         PRINT("(2(A,F20.15))"), "   r = ",r(1), "   s = ", s(1)
              
        
        IF (ABS(f) < tol .AND. ABS(g) < tol) THEN
          EXIT
        ENDIF        
       
        
      ENDDO
      
!       error = max(abs(f),abs(g))
!       IF (it >= maxit) THEN
!         PRINT("(A,E22.15)"), "   MAX ITERATIONS EXCEEDED, error = ",error
!         PRINT("(2(A,F20.15))"), "   r = ",r(1), "   s = ", s(1)
!       ELSE       
!         PRINT("(A,I7,A,E22.15)"), "   iterations: ",it, "  error = ",error
!         PRINT("(2(A,F20.15))"), "   r = ",r(1), "   s = ", s(1)
!       ENDIF


      RETURN
      END SUBROUTINE newton
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      


      END MODULE find_element