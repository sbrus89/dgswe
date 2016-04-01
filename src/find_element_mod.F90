      MODULE find_element

      USE globals, ONLY: rp
      USE kdtree2_module
      USE lapack_interfaces
      
      IMPLICIT NONE
      
      SAVE
      
      INTEGER, PARAMETER :: srchdp = 20
      REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: rsre          
      TYPE(kdtree2), POINTER :: tree_xy      
      TYPE(kdtree2_result), ALLOCATABLE, DIMENSION(:) :: closest   
      INTEGER, DIMENSION(:), ALLOCATABLE :: nnd2el
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: nd2el
      INTEGER, DIMENSION(:), ALLOCATABLE :: nverts
      INTEGER, DIMENSION(:), ALLOCATABLE :: np
      INTEGER, DIMENSION(:), ALLOCATABLE :: nnds      
      REAL(rp), DIMENSION(:), ALLOCATABLE :: elnx,elny
      REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: V
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: ipiv
      INTEGER :: mnnds
      

      CONTAINS
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   

      SUBROUTINE find_element_init(nel_type,nvertex,nporder,nnodes,nn,xy,nepn,epn)
      
      USE vandermonde, ONLY: vandermonde_area
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: nel_type
      INTEGER, DIMENSION(:), INTENT(IN) :: nvertex
      INTEGER, DIMENSION(:), INTENT(IN) :: nporder
      INTEGER, DIMENSION(:), INTENT(IN) :: nnodes      
      INTEGER, INTENT(IN) :: nn
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: xy      
      INTEGER, DIMENSION(:), INTENT(IN) :: nepn
      INTEGER, DIMENSION(:,:), INTENT(IN) :: epn
            
      INTEGER :: nd,el,et,n,p
      INTEGER :: mnepn,mnp
      INTEGER :: info
      
      ! initialize module varibles to increase flexibility and decrease number
      ! of subroutine calling arguments
      
      tree_xy => kdtree2_create(xy(1:2,1:nn), rearrange=.true., sort=.true.)
      ALLOCATE(closest(srchdp))      
            
            
      CALL ref_elem_coords(nel_type,rsre)
            
            
      mnepn = MAXVAL(nepn)
      ALLOCATE(nnd2el(nn),nd2el(mnepn,nn))      
      DO nd = 1,nn
        nnd2el(nd) = nepn(nd)
        DO el = 1,nnd2el(nd)
          nd2el(el,nd) = epn(el,nd)
        ENDDO      
      ENDDO
  
  
      ALLOCATE(nverts(nel_type),np(nel_type),nnds(nel_type))      
      DO et = 1,nel_type
        nverts(et) = nvertex(et)
        np(et) = nporder(et)
        nnds(et) = nnodes(et)
      ENDDO
      
      
      mnp = maxval(np)            
      mnnds = (mnp+1)**2
      ALLOCATE(elnx(mnnds),elny(mnnds))
      
      
      ALLOCATE(V(mnnds,mnnds,nel_type),ipiv(mnnds,nel_type))
      DO et = 1,nel_type
        p = np(et)
        CALL vandermonde_area(et,p,n,V(:,:,et))
        CALL DGETRF(n,n,V(1,1,et),mnnds,ipiv(1,et),info)           
      ENDDO
          
      
      RETURN
      END SUBROUTINE find_element_init
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      

      SUBROUTINE in_element(xt,el_type,elxy,el_found,leds)

      IMPLICIT NONE
              
      REAL(rp), INTENT(IN) :: xt(2)
      INTEGER, DIMENSION(:), INTENT(IN) :: el_type
      REAL(rp), DIMENSION(:,:,:), INTENT(IN) :: elxy
      INTEGER, INTENT(OUT) :: el_found          
      INTEGER, INTENT(OUT), OPTIONAL :: leds(4)
            
      INTEGER :: srch,nd
      INTEGER :: el,eln,clnd,et,n
      INTEGER :: found     
      INTEGER :: min_el
      REAL(rp) :: diff,min_diff
      REAL(rp) :: tol
            
      
      tol = 1d-5 
        CALL kdtree2_n_nearest(tp=tree_xy,qv=xt,nn=srchdp,results=closest) ! find what element xt is in               
        
        ! Test elements to see which element point is located in    
        found = 0      
        min_diff = 999d0
search: DO srch = 1,srchdp

          clnd = closest(srch)%idx
          
          PRINT("(A,I5)"), "   closest node: ", clnd

   elem:  DO el = 1,nnd2el(clnd) 
 
            eln = nd2el(el,clnd)
            
            PRINT("(A,I5)"), "   testing: ", eln                               

            ! Compute sum of sub-triangle areas
            CALL sub_element(xt,eln,el_type,elxy,diff,leds)            
            
            IF (diff < min_diff) THEN  ! keep track of element with minimum difference in sum of sub-triangle areas                                  
              min_diff = diff          ! to return if tolerance is not met
              min_el = eln
            ENDIF
          
            ! The station is in the element if the reference element area and sum of sub triangle are the same
            IF (diff < tol) THEN
              PRINT("(A,I5)"), "   element found", eln                            
                      
              el_found = eln        
              found = 1                        

              EXIT search            
            ENDIF    
            
!             PRINT*, " "  
          ENDDO elem
        
        ENDDO search    
        

        
        IF (found == 0) THEN
          el_found = min_el        

          PRINT*, "ELEMENT NOT FOUND"  
          PRINT*, "USING ELEMENT ", el_found, "(AREA DIFF = ",min_diff, ")"       

          CALL sub_element(xt,el_found,el_type,elxy,diff,leds)           
          
        ENDIF     
 

      RETURN
      END SUBROUTINE in_element
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 

      SUBROUTINE sub_element(xt,eln,el_type,elxy,diff,closest_ed)   
      
      USE transformation, ONLY: xy2rs
      
      IMPLICIT NONE
      
      REAL(rp), INTENT(IN) :: xt(2)        
      INTEGER, INTENT(IN) :: eln
      INTEGER, DIMENSION(:), INTENT(IN) :: el_type
      REAL(rp), DIMENSION(:,:,:), INTENT(IN) :: elxy          
      REAL(rp), INTENT(OUT) :: diff      
      INTEGER, INTENT(OUT) :: closest_ed(4)
      
      INTEGER :: i,j    
      INTEGER :: nv,et,n,npt,nd,p
      INTEGER :: n1,n2
      INTEGER :: etemp
      REAL(rp) :: area,area_sum
      REAL(rp) :: stri_area(4),stri_min
      REAL(rp) :: dist
      REAL(rp) :: x(1),y(1)
      REAL(rp) :: r(1),s(1)
      REAL(rp) :: x1,x2,x3,y1,y2,y3
      REAL(rp) :: atemp
      
      
      et = el_type(eln)
      nv = nverts(et)
      p = np(et)
      n = nnds(et)
      npt = 1

      
      DO nd = 1,n
        elnx(nd) = elxy(nd,eln,1)
        elny(nd) = elxy(nd,eln,2)
      ENDDO
      
      x(1) = xt(1)
      y(1) = xt(2)
                                                        
      ! Compute the local (r,s) coordinates of the (x,y) point
!       CALL xy2rs(et,p,elnx,elny,npt,x,y,r,s)
      CALL rs_coords(et,p,elnx,elny,x,y,r,s)
          
      ! Find reference element area
      IF (mod(et,2) == 1) THEN
        area = 2d0
      ELSE IF (mod(et,2) == 0) THEN
        area = 4d0
      ENDIF          
          
      ! Compute sum of sub-triangle areas
      area_sum = 0d0
      stri_min = 999d0
      DO i = 1,nv
        n1 = mod(i+0,nv)+1
        n2 = mod(i+1,nv)+1           

        x1 = rsre(1,n1,et)
        y1 = rsre(2,n1,et)
            
        x2 = rsre(1,n2,et)
        y2 = rsre(2,n2,et)
            
        x3 = r(1)
        y3 = s(1)
                        
        stri_area(i) = .5d0*abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1))
        closest_ed(i) = i
                                            
        area_sum = area_sum + stri_area(i)
      ENDDO
          
!     PRINT("(A,F20.15,A,F20.15)"), "   area = ",area, "   area sum = ", area_sum
!     PRINT*, " "
          
      diff = abs(area-area_sum) 
      
      DO i = 1,nv           ! keep track of minimum sub-triangle area to determine
        DO j = i+1,nv       ! which edge the point lies on, or is closest to
          IF (stri_area(j) < stri_area(i)) THEN
            atemp = stri_area(i)
            stri_area(i) = stri_area(j)
            stri_area(j) = atemp
            
            etemp = closest_ed(i)
            closest_ed(i) = closest_ed(j)
            closest_ed(j) = etemp            
          ENDIF
        ENDDO
      ENDDO
            
      
      RETURN
      END SUBROUTINE sub_element    
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
      
      
      SUBROUTINE rs_coords(et,p,elx,ely,x,y,r,s)

      USE basis, ONLY: element_basis,tri_basis,quad_basis

      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: et
      INTEGER, INTENT(IN) :: p
      REAL(rp), DIMENSION(:), INTENT(IN) :: elx
      REAL(rp), DIMENSION(:), INTENT(IN) :: ely
      REAL(rp), INTENT(IN) :: x(1)  
      REAL(rp), INTENT(IN) :: y(1)        
      REAL(rp), INTENT(OUT) :: r(1)    
      REAL(rp), INTENT(OUT) :: s(1)      
      
      INTEGER :: it,n,i,j
      INTEGER :: info
      INTEGER :: maxit
      REAL(rp) :: tol
      REAL(rp) :: f,g,error
      REAL(rp) :: dfdr,dfds,dgdr,dgds,jac_recip
      REAL(rp) :: phi(mnnds,1),dpdr(mnnds,1),dpds(mnnds,1)
      REAL(rp) :: l(mnnds,3)
        
      tol = 1d-9
      maxit = 100
      info = 0
     
        
      IF (mod(et,2) == 1) THEN
        r(1) = -1d0/3d0
        s(1) = -1d0/3d0
      ELSE IF (mod(et,2) == 0) THEN
        r(1) = 0d0
        s(1) = 0d0
      ENDIF

      DO it = 1,maxit                        
              
        CALL element_basis(et,p,n,1,r,s,phi,dpdr,dpds)
        
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

          dfdr = dfdr + l(i,2)*elx(i)
          dfds = dfds + l(i,3)*elx(i)
          dgdr = dgdr + l(i,2)*ely(i)
          dgds = dgds + l(i,3)*ely(i)
          
          f = f + l(i,1)*elx(i)
          g = g + l(i,1)*ely(i)

        ENDDO
                
        jac_recip = 1d0/(dfdr*dgds - dgdr*dfds)

        f = f - x(1)
        g = g - y(1)
        
        r(1) = r(1) - ( dgds*f - dfds*g)*jac_recip
        s(1) = s(1) - (-dgdr*f + dfdr*g)*jac_recip   
!         PRINT("(3(A,F20.15))"), "   f = ",f, "   g = ", g, "  jac = ", jac        
!         PRINT("(2(A,F20.15))"), "   r = ",r(1), "   s = ", s(1)
              
        
        IF (ABS(f) < tol .AND. ABS(g) < tol) THEN
          EXIT
        ENDIF        
       
        
      ENDDO
      
      error = max(abs(f),abs(g))
      
!       IF (it >= maxit) THEN
!         PRINT("(A,E22.15)"), "   MAX ITERATIONS EXCEEDED, error = ",error
!         PRINT("(2(A,F20.15))"), "   r = ",r(1), "   s = ", s(1)
!       ELSE       
!         PRINT("(A,I7,A,E22.15)"), "   iterations: ",it, "  error = ",error
!         PRINT("(2(A,F20.15))"), "   r = ",r(1), "   s = ", s(1)
!       ENDIF

      RETURN
      END SUBROUTINE rs_coords      
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   


      SUBROUTINE ref_elem_coords(nel_type,rsre)
      
      USE basis, ONLY: element_nodes
      
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nel_type
      REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: rsre
      
      INTEGER :: et,pt,n
      INTEGER :: np,space
      REAL(rp) :: r(4),s(4)
      
      ALLOCATE(rsre(2,4,nel_type))
      
      space = 1
      np = 1
      
      DO et = 1,nel_type      
        
        CALL element_nodes(et,space,np,n,r,s)
        
        DO pt = 1,n
          rsre(1,pt,et) = r(pt)
          rsre(2,pt,et) = s(pt)
        ENDDO
        
      ENDDO
      
      END SUBROUTINE ref_elem_coords
      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   


      SUBROUTINE check_elem(xt,el_found)      
      
      IMPLICIT NONE
      

      REAL(rp), INTENT(IN) :: xt(2)
      INTEGER, INTENT(IN) :: el_found
            
      INTEGER :: el,eln,clnd
      INTEGER :: found      

                  
      CALL kdtree2_n_nearest(tp=tree_xy,qv=xt,nn=1,results=closest)              
        
  
      found = 0      
      clnd = closest(1)%idx          

elem: DO el = 1,nnd2el(clnd) 
 
        eln = nd2el(el,clnd)
            
        IF (eln == el_found) THEN
          found = 1
          EXIT elem
        ENDIF
            
      ENDDO elem
        
        
      IF (found == 0) THEN
        PRINT*, "ERROR FINDING ELEMENT FOR VERTEX NODE"
        STOP
      ENDIF     
 
      
      
      RETURN
      END SUBROUTINE check_elem

      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       


      END MODULE find_element