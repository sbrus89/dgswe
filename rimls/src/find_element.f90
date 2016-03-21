      MODULE find_element

      USE globals, ONLY: rp
      USE kdtree2_module

      CONTAINS
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      

      SUBROUTINE in_element(pt,xt,eln,found)

      USE globals, ONLY: base,tree_xy,srchdp,closest,nverts,rsre

      IMPLICIT NONE
      
      INTEGER :: srch,i,pt
      INTEGER :: eln,et,nvert
      INTEGER :: found,el_found      
      INTEGER :: n1,n2
      INTEGER :: n,local_el(srchdp)
      INTEGER :: k(1)
      REAL(rp) :: xt(2),x(3),y(3),r(1),s(1)
      REAL(rp) :: sarea,area,diff(srchdp)
      REAL(rp) :: tol
      
      tol = 1d-5 
        CALL kdtree2_n_nearest(tp=tree_xy,qv=xt,nn=srchdp,results=closest) ! find what element xt is in               
        
        ! Test elements to see which element point is located in    
        found = 0      
        n = 0
        diff = 999d0
search: DO srch = 1,srchdp
          eln = closest(srch)%idx    

            et = base%el_type(eln)
            nvert = nverts(et)   
            
            n = n+1
            local_el(n) = eln ! keep track of elements (and sum of sub-triangle areas) 
                              ! to find minimum if tolerance is not met                      
          
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
          
!               PRINT("(A,I5,A,F20.15,A,F20.15)"), "   testing: ", eln, "   area = ",area, "   sarea = ", sarea
!               PRINT*, " "
          
            diff(n) = abs(area-sarea)
          
            ! The station is in the element if the reference element area and sum of sub triangle are the same
            IF (diff(n) < tol) THEN
!               PRINT("(A,I5)"), "   element found", eln
                      
              el_found = eln        
              found = 1                        
            
              EXIT search                        
            ENDIF       
        
        ENDDO search    
        

        
        IF (found == 0) THEN
          k = minloc(diff)
          eln = local_el(k(1))        

          PRINT*, "ELEMENT NOT FOUND FOR POINT: ",pt  
          PRINT*, "USING ELEMENT ", eln, "(AREA DIFF = ",diff(k(1)), ")"       
          
        ELSE         
         eln = el_found       
        ENDIF     
 

      RETURN
      END SUBROUTINE in_element
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
      SUBROUTINE newton(pt,x,y,eln,r,s)

      USE globals, ONLY: rp,np,nverts,nnds,mnnds,V,base,ipiv
      USE basis, ONLY: tri_basis,quad_basis
      USE shape_functions_mod, ONLY: shape_functions_area_eval
      USE transformation, ONLY: element_transformation

      IMPLICIT NONE
      INTEGER :: it,eln,et,p,n,i,pt,nv
      INTEGER :: info
      INTEGER :: maxit
      REAL(rp) :: tol
      REAL(rp) :: x,y
      REAL(rp) :: r(1),s(1)
      REAL(rp) :: f,g,error
      REAL(rp) :: drdf,drdg,dsdf,dsdg,jac
      REAL(rp) :: l(mnnds,1),dldr(mnnds,1),dlds(mnnds,1)
        
      tol = 1d-9
      maxit = 100
      info = 0
      
      et = base%el_type(eln)
      p = np(et)  
      nv = nverts(et)
        
      IF (mod(et,2) == 1) THEN
        r(1) = -1d0/3d0
        s(1) = -1d0/3d0
      ELSE IF (mod(et,2) == 0) THEN
        r(1) = 0d0
        s(1) = 0d0
      ENDIF

      DO it = 1,maxit     
           
        CALL shape_functions_area_eval(nv,p,n,1,r,s,l,dldr,dlds)     
        
        CALL element_transformation(n,base%elxy(:,eln,1),base%elxy(:,eln,2),l(:,1),f,g, &
                                    dldr(:,1),dlds(:,1),drdf,drdg,dsdf,dsdg,jac)               
        
        f = f - x
        g = g - y
        
        r(1) = r(1) - (drdf*f + drdg*g)
        s(1) = s(1) - (dsdf*f + dsdg*g)        
        
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
      END SUBROUTINE newton
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      


      END MODULE find_element