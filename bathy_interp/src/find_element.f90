      MODULE find_element

      USE globals, ONLY: rp
      USE kdtree2_module

      CONTAINS
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      

      SUBROUTINE in_element(xt,eln,rt,error,exceed)

      USE globals, ONLY: base,tree_xy,srchdp,closest,mnepn

      IMPLICIT NONE
      
      INTEGER, INTENT(OUT) :: eln
      INTEGER :: srch,i,n
      INTEGER :: k(1)
      INTEGER :: ndc,elc
      INTEGER :: et,nvert
      INTEGER :: found,el_found      
      INTEGER :: n1,n2
      INTEGER :: local_el(srchdp*mnepn)
      INTEGER :: exceed
      REAL(rp), INTENT(IN) :: xt(2)
      REAL(rp), INTENT(OUT) :: rt(2)
      REAL(rp) :: x(3),y(3)
      REAL(rp) :: r(srchdp*mnepn),s(srchdp*mnepn)
      REAL(rp) :: sarea,area,diff(srchdp*mnepn)
      REAL(rp) :: tol
      REAL(rp) :: error
      
      tol = 1d-5 
      CALL kdtree2_n_nearest(tp=tree_xy,qv=xt,nn=srchdp,results=closest) ! find what element xt is in               
        
        ! Test elements to see which element point is located in    
        found = 0    
        n = 0
        diff = 999d0
search: DO srch = 1,srchdp
          ndc = base%vxyn(closest(srch)%idx)
          
          DO elc = 1,base%nepn(ndc)          

            eln = base%epn(elc,ndc)
            et = base%el_type(eln)
            nvert = base%nverts(et)    
            
            n = n+1
            local_el(n) = eln ! keep track of elements (and sum of sub-triangle areas) 
                              ! to find minimum if tolerance is not met
          
            ! Compute the local (r,s) coordinates of the (x,y) station location
            CALL newton(xt(1),xt(2),eln,r(n),s(n),error,exceed)
          
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

              x(1) = base%rsre(1,n1,et)
              y(1) = base%rsre(2,n1,et)
            
              x(2) = base%rsre(1,n2,et)
              y(2) = base%rsre(2,n2,et)
            
              x(3) = r(n)
              y(3) = s(n)
            
              sarea = sarea + .5d0*abs((x(2)-x(1))*(y(3)-y(1)) - (x(3)-x(1))*(y(2)-y(1)))
            ENDDO
          
!               PRINT("(A,I8,A,F20.15,A,F20.15)"), "   testing: ", eln, "   area = ",area, "   sarea = ", sarea
!               PRINT*, " "
          
            diff(n) = abs(area-sarea)
          
            ! The station is in the element if the reference element area and sum of sub triangle are the same            
            IF (diff(n) < tol) THEN
!               PRINT("(A,I8)"), "   element found", eln
                      
              el_found = eln        
              found = 1                        
            
              EXIT search                        
            ENDIF    
            
          ENDDO
        
        ENDDO search    
        

        
        IF (found == 0) THEN      
          k = minloc(diff)
          eln = local_el(k(1))
          rt(1) = r(k(1))
          rt(2) = s(k(1))              
          PRINT*, "ELEMENT NOT FOUND"   
          PRINT*, "USING ELEMENT ", eln, "(AREA DIFF = ",diff(k(1)), ")"                 
        ELSE         
          eln = el_found 
          rt(1) = r(n)
          rt(2) = s(n)
        ENDIF     
 

      RETURN
      END SUBROUTINE in_element
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
      SUBROUTINE newton(x,y,eln,r,s,error,exceed)

      USE globals, ONLY: base
      USE basis, ONLY: tri_basis,quad_basis

      IMPLICIT NONE
      INTEGER :: it,eln,et,p,n,i
      INTEGER :: info
      INTEGER :: maxit,exceed
      REAL(rp) :: tol
      REAL(rp) :: x,y
      REAL(rp) :: r(1),s(1)
      REAL(rp) :: f,g,error
      REAL(rp) :: dfdr,dfds,dgdr,dgds,jac
      REAL(rp) :: phi(base%mnnds,1),dpdr(base%mnnds,1),dpds(base%mnnds,1)
      REAL(rp) :: l(base%mnnds,3)
        
      tol = 1d-11
      maxit = 100
      exceed = 0
      info = 0
      
      et = base%el_type(eln)
      p = base%np(et)  
        
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
          

        CALL DGETRS("N",n,3,base%V(1,1,et),base%mnnds,base%ipiv(1,et),l,base%mnnds,info)
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
        
        IF (ABS(f) < tol .AND. ABS(g) < tol) THEN
          EXIT
        ENDIF        
       
        
      ENDDO
      
      error = max(abs(f),abs(g))
      
      IF (it >= maxit) THEN
!         PRINT("(A,E22.15)"), "   MAX ITERATIONS EXCEEDED, error = ",error
!         PRINT("(2(A,F20.15))"), "   r = ",r(1), "   s = ", s(1)
          exceed = 1
      ELSE       
!         PRINT("(A,I7,A,E22.15)"), "   iterations: ",it, "  error = ",error
!         PRINT("(2(A,F20.15))"), "   r = ",r(1), "   s = ", s(1)
      ENDIF


      RETURN
      END SUBROUTINE newton
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      


      END MODULE find_element