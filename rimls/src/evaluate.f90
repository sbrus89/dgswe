      MODULE evaluate
      
      USE globals, ONLY: rp,r,sigma_n
      USE lapack_interfaces
      
      IMPLICIT NONE
      
      INTEGER :: maxit,maxptit
      REAL(rp) :: sigma_r,threshold,percent_max      
      
      CONTAINS
      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE compute_surface()
      
      USE globals, ONLY: refinement,base,eval,hmin,hmax
      USE kdtree2_module      
      
      IMPLICIT NONE

      
      maxit = 100
      maxptit = 100
      threshold = 1d-12
      percent_max = 100d0
      sigma_r = 0.5d0 
!       sigma_n = 1.5d0   ! 0.5 - 1.5
!       r = 15d0           ! 1.5 - 4.0      
      
      hmin = minval(base%depth)
      hmax = maxval(base%depth)
       

      
!       CALL filter_normals()
      

        PRINT("(A)"), "Computing rimls surface: verticies"
        CALL mls_surface(eval%nn,eval%npts_vertex,eval%xyh_vertex)      
        PRINT("(A)"), "Computing rimls surface: edges"      
        CALL mls_surface(eval%ned,eval%npts_edge,eval%xyh_edge)
        PRINT("(A)"), "Computing rimls surface: interior"
        CALL mls_surface(eval%ne,eval%npts_interior,eval%xyh_interior)      

!         PRINT("(A)"), "Computing rimls surface: verticies"
!         CALL rimls_surface(eval%nn,1,1,eval%xyhv)      
!         PRINT("(A)"), "Computing rimls surface: edges"      
!         CALL rimls_surface(eval%ned,eval%np(6)-1,eval%mnnds,eval%xyhe)
!         PRINT("(A)"), "Computing rimls surface: interior"
!         CALL rimls_surface(eval%ne,eval%mninds,eval%mnnds,eval%xyhi)      
! 
!         PRINT("(A)"), "Computing rimls surface: verticies"
!         CALL function_surface(eval%nn,1,1,eval%xyhv)      
!         PRINT("(A)"), "Computing rimls surface: edges"      
!         CALL function_surface(eval%ned,eval%np(6)-1,eval%mnnds,eval%xyhe)
!         PRINT("(A)"), "Computing rimls surface: interior"
!         CALL function_surface(eval%ne,eval%mninds,eval%mnnds,eval%xyhi)      

      RETURN
      END SUBROUTINE
      
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
      
      SUBROUTINE mls_surface(n,npts,xyh_eval)
      
      USE globals, ONLY: base,hmin,hmax,lsp,basis_opt
      USE find_element, ONLY: in_element,tree_xy,closest
      USE kdtree2_module
      USE basis, ONLY: element_basis,simple_basis
      USE transformation, ONLY: xy2rs
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: n
      INTEGER, DIMENSION(:), INTENT(IN) :: npts
      REAL(rp), DIMENSION(:,:,:), INTENT(INOUT) :: xyh_eval
      
      
      INTEGER :: i,j,pt,nd,l,m,k       
      INTEGER :: ne,el,elin,et
      INTEGER :: ndf,nnd,p,nhnd
      INTEGER :: neighbors(base%ne),nnd_neigh,nel_neigh,el_flag(base%ne)
      INTEGER :: small_flag
      INTEGER :: info,lwork,lda,ldb
      REAL(rp) :: elx(base%mnnds),ely(base%mnnds)
      REAL(rp) :: xpt(3),xnd(3),xnd_nrm(3),xbar(2),rs(2)      
      REAL(rp) :: w,rhs,lhs
      REAL(rp) :: s,t,sv(1),tv(1),x(1),y(1)
      REAL(rp) :: hpt,search_r,tol      
      REAL(rp) :: nwork(1)
      REAL(rp), ALLOCATABLE :: A(:,:),b(:),work(:)
      REAL(rp), ALLOCATABLE :: phi(:,:)
      
      ne = base%ne
      tol = 1d-16
      ndf = (lsp+1)**2
      
      lda = ne*base%mninds
      ldb = lda
      
      ALLOCATE(A(lda,ndf),b(ldb))
      CALL DGELS('N',ne,ndf,1,A,lda,b,ldb,nwork,-1,info)  ! find the optimal work array size
                                                          ! using ne is likely a huge overestimate     
                                                        
      ALLOCATE(phi(ndf,1))                                                        
      
      IF (info == 0) THEN
        lwork = INT(nwork(1))
        ALLOCATE(work(lwork))
      ELSE 
        STOP
      ENDIF
      
      DO i = 1,n
      
        IF (mod(i,1000) == 0) THEN       
          PRINT*,i,"/",n
        ENDIF   
        
        DO pt = 1,npts(i) 
        
          xpt(1) = xyh_eval(pt,i,1)
          xpt(2) = xyh_eval(pt,i,2)
          xpt(3) = xyh_eval(pt,i,3)
          
!           PRINT*, "pt = ", pt, "x = ",xpt(1), "y = ",xpt(2)         
          
          CALL in_element(xpt(1:2),base%el_type,base%elxy,elin,rs)          
          et = base%el_type(elin)       
          p = base%np(et)
          
          xbar(1) = base%elxy_center(1,elin,1)
          xbar(2) = base%elxy_center(1,elin,2)          
          
          DO nd = 1,base%nnds(et)
            elx(nd) = base%elxy(nd,elin,1)
            ely(nd) = base%elxy(nd,elin,2)
          ENDDO
   
          hpt = r*base%h(elin)           
          search_r = SQRT(-hpt**2*LOG(tol))
          CALL kdtree2_r_nearest(tp=tree_xy,qv=xpt(1:2),r2=search_r**2,nfound=nnd_neigh,nalloc=base%nn,results=closest)
          
          el_flag = 0
          nel_neigh = 0
          DO k = 1,nnd_neigh
            nd = closest(k)%idx
            DO j = 1,base%nepn(nd)
              el = base%epn(j,nd)
              
              IF (el_flag(el) == 0) THEN
                nel_neigh = nel_neigh + 1
                neighbors(nel_neigh) = el
                el_flag(el) = 1
              ENDIF
            ENDDO
          ENDDO

          
          CALL boundary_check2(i,elin,nel_neigh,neighbors)
          



          small_flag = 0
          nnd = 0
          DO j = 1,nel_neigh
                      
            el = neighbors(j)
            
            IF (mod(base%el_type(el),2) == 1) THEN
              nhnd = base%nnds(5)
            ELSE IF (mod(base%el_type(el),2) == 0) THEN
              nhnd = base%nnds(6)
            ENDIF
                        
            DO nd = 1,nhnd
            
              nnd = nnd + 1
              
              xnd(1) = base%elxyh(nd,el,1)
              xnd(2) = base%elxyh(nd,el,2)
              xnd(3) = base%elhb(nd,el)
            
!             xnd_nrm = base%nhb(:,el)  
            
              w = theta(xpt,xnd,hpt)
            
              IF (w <= tol*10d0) THEN
                small_flag = 1
              ENDIF
            
              x(1) = xnd(1)
              y(1) = xnd(2)            
              
              IF (basis_opt == 1) THEN
                CALL xy2rs(et,p,elx,ely,1,x,y,sv,tv)             
                CALL element_basis(et,lsp,ndf,1,sv,tv,phi)
              ELSE
                sv(1) = x(1) - xbar(1)
                tv(1) = y(1) - xbar(2)
                CALL simple_basis(lsp,ndf,1,sv,tv,phi)
              ENDIF
              
              DO m = 1,ndf
                A(nnd,m) = w*phi(m,1)
              ENDDO
            
              b(nnd) = w*xnd(3)            
            ENDDO
          ENDDO
          
          
!           IF (small_flag /= 1) THEN
!             PRINT*, "Increase search radius"
!             STOP
!           ENDIF
          
          CALL DGELS('N',nnd,ndf,1,A,lda,b,ldb,work,lwork,info)          

          x(1) = xpt(1)
          y(1) = xpt(2)
          
          IF (basis_opt == 1) THEN
            CALL xy2rs(et,p,elx,ely,1,x,y,sv,tv)      
            CALL element_basis(et,lsp,ndf,1,sv,tv,phi)          
          ELSE
            sv(1) = x(1) - xbar(1)
            tv(1) = y(1) - xbar(2)          
            CALL simple_basis(lsp,ndf,1,sv,tv,phi)          
          ENDIF
          
          xyh_eval(pt,i,3) = 0d0
          DO m = 1,ndf
              xyh_eval(pt,i,3) = xyh_eval(pt,i,3) + b(m)*phi(m,1)
          ENDDO    
              
!           IF (xyh_eval(pt,i,3) < hmin*.9d0) THEN
!             PRINT "(A,I7,A,I7,A,F15.7)", "Value exceeds minimum tolerance, point: ", i, " el: ", elin, " hb = ", xyh_eval(pt,i,3)
!             xyh_eval(pt,i,3) = hmin*.9d0
!           ENDIF
          
        ENDDO
      
        
      
      ENDDO
      
!       CALL invcpp(n,npts,mnpts,xyh_eval)      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!   Constant Least Squares Fit !!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           rhs = 0d0
!           lhs = 0d0
!           small_flag = 0
!           DO nd = 1,nneigh
!           
!             el = neighbors(nd)
!             xnd  = base%xyhc(:,el)
!             xnd_nrm = base%nhb(:,el)  
!             
!             w = theta(xpt,xnd,hpt)
!             
!             IF (w < 1d-10) THEN
!               small_flag = 1
!             ENDIF
!             
!             rhs = rhs + w**2*xnd(3)
!             lhs = lhs + w**2
!           
!           ENDDO
!           
!           
!           IF (small_flag /= 1) THEN
!             PRINT*, "Increase search radius"
!             STOP
!           ENDIF
!           
!           xyh_eval(pt,i,3) = rhs/lhs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      RETURN
      END SUBROUTINE mls_surface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

! 
! 
!       SUBROUTINE rimls_surface(n,npts,mnpts,xyh)
!       
!       USE globals, ONLY: tree_xy,tree_c,kdresults,base,hmin,hmax,phi0,lambda0,Erad
!       USE find_element, ONLY: in_element
!       USE kdtree2_module
!       
!       IMPLICIT NONE
!       
!       INTEGER :: i,j,pt,it,ptit
!       INTEGER :: n,npts,mnpts
!       INTEGER :: nneigh,nneigh_pre
!       INTEGER :: nd,el,elin
!       INTEGER :: neighbors(base%ne),leds(4)
!       REAL(rp) :: search_r
!       REAL(rp) :: xyh(mnpts,n,3)
!       REAL(rp) :: xpt(3),p(3),px(3),np(3),hpt,f,grad_f(3),grad_w(3),fgradf(3),npmgradf(3)
!       REAL(rp) :: alpha(base%ne),alpha_old(base%ne)
!       REAL(rp) :: sumW,sumF
!       REAL(rp) :: sumGw(3),sumGF(3),sumN(3)
!       REAL(rp) :: w,fx
!      
!       
! 
!       
!       DO i = 1,n
!        
!         IF (mod(i,1000) == 0) THEN       
!           PRINT*,i,"/",n
!         ENDIF
!         DO pt = 1,npts  
!         
!           xpt(1) = xyh(pt,i,1)
!           xpt(2) = xyh(pt,i,2)
!           xpt(3) = xyh(pt,i,3)          
!                 
!           
!           CALL in_element(xpt(1:2),base%el_type,base%elxy,elin,leds)           
!    
!           hpt = r*base%h(elin) 
!           search_r = SQRT(-hpt**2*LOG(threshold))          
!           CALL kdtree2_r_nearest(tp=tree_xy,qv=xpt(1:2),r2=search_r**2,nfound=nneigh,nalloc=base%ne,results=kdresults)
!           
! 
!           CALL boundary_check2(i,elin,nneigh,neighbors)
!                  
!           
!           grad_f(:) = 1d0
!           f = 1d0          
! 
!           ptit = 0
!           fgradf = f*grad_f
!    wloop: DO WHILE (norm(fgradf) > threshold)
!             it = 0
!                         
!             DO j = 1,nneigh
!               alpha(j) = 1d0
!               alpha_old(j) = 0d0
!             ENDDO
!             
!             DO WHILE (it < maxit)
!               
!               sumW = 0d0
!               sumF = 0d0
!               sumGw(:) = 0d0
!               sumGF(:) = 0d0
!               sumN(:) = 0d0
!             
!               DO nd = 1,nneigh
!               
! 
!                 el = neighbors(nd)
!                 p  = base%xyhc(:,el)
!                 np = base%nhb(:,el)
!                 
!                 px = xpt - p
!                 
!                 fx = DOT_PRODUCT(px,np)
!                 
!                 IF (it > 0) THEN
!                   npmgradf = np-grad_f
!                   alpha(nd) = exp(-((fx-f)/(sigma_r*hpt))**2)*exp(-(norm(npmgradf)/sigma_n)**2)
!                 ELSE
!                   alpha(nd) = 1d0
!                 ENDIF                           
!                 
!                 w = alpha(nd)*phi(px,hpt)
!                 grad_w = alpha(nd)*dphi(px,hpt)               
!                 
!                 sumW = sumW + w
!                 sumGw = sumGw + grad_w
!                 sumF = sumF + w*fx
!                 sumGF = sumGF + grad_w*fx
!                 sumN = sumN + w*np                               
!                 
!               ENDDO
!               
!               f = sumF/sumW
!               grad_f = (sumGF - f*sumGw + sumN)/sumW
!                             
!               it = it + 1              
!               
!               
!               DO j = 1,nneigh
!                 alpha_old(j) = alpha(j)
!               ENDDO
!               
!             ENDDO
!             
!             fgradf = f*grad_f
!             xpt = xpt - fgradf
!             ptit = ptit + 1
!             
!             IF (ptit == maxptit) THEN
!               EXIT wloop
!             ENDIF
!             
!           ENDDO wloop
!                 
!           IF (ptit < maxptit) THEN
!           
!             IF (xpt(1) /= xpt(1) .or. xpt(2) /= xpt(2) .or. xpt(3) /= xpt(3)) THEN
!           
!               PRINT*, "WARNING: Nan detected, i = ",i
!               PRINT*, "    # neighbors = ", nneigh
!                          
!               DO nd = 1,nneigh
!                 PRINT*, "    ", neighbors(nd)
!               ENDDO
!               
!             ELSE IF (xpt(3) < hmin) THEN 
!             
!               xyh(pt,i,3) = hmin 
! !               PRINT*, "WARNING: Limiting hmin"
!               
!             ELSE IF (xpt(3) > hmax) THEN
!             
!               xyh(pt,i,3) = hmax
! !               PRINT*, "WARNING: Limiting hmax"              
!               
! !             ELSE IF (abs((xpt(3)-xyh(pt,i,3))/xyh(pt,i,3))*100d0 > percent_max) THEN
! !             
! !               PRINT*, "WARNING: Percent change exceedes tolerance, i = ",i
!               
!             ELSE
!             
! !               print*, norm(fgradf)
!               xyh(pt,i,1) = xpt(1)
!               xyh(pt,i,2) = xpt(2)
!               xyh(pt,i,3) = xpt(3)               
! 
!             ENDIF
!             
!           ELSE
!             PRINT*, "WARNING: max point iterations reached, i = ",i
!             STOP
!           ENDIF
!       
!         ENDDO
!       ENDDO
!       
! !       CALL invcpp(n,npts,mnpts,xyh)
!       
!       
!       RETURN
!       END SUBROUTINE rimls_surface     



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      SUBROUTINE function_surface(n,npts,mnpts,xyh)
      
      USE globals, ONLY: pi
      
      IMPLICIT NONE
      
      INTEGER :: i,pt
      INTEGER :: n,npts,mnpts       
      REAL(rp) :: xyh(mnpts,n,3) 
      REAL(rp) :: x(3)
      
      DO i = 1,n
      
        IF (mod(i,1000) == 0) THEN       
          PRINT*,i,"/",n
        ENDIF   
        
        DO pt = 1,npts  
        
          x(1) = xyh(pt,i,1)
          x(2) = xyh(pt,i,2)
          x(3) = xyh(pt,i,3)   
          
          xyh(pt,i,3)  = 10d0 - 5d0*cos(2d0*pi/500d0*x(2))  
          
        ENDDO
        
      ENDDO
      
      
      END SUBROUTINE      
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
      FUNCTION theta(x,p,h) RESULT(w)
      
      IMPLICIT NONE
      
      REAL(rp) :: w
      REAL(rp), INTENT(IN) :: x(3),p(3),h
      REAL(rp) :: d
      
      d = sqrt((x(1)-p(1))**2 + (x(2)-p(2))**2)
      w = exp(-d**2/h**2)      
      
      END FUNCTION theta              
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      FUNCTION norm(a) RESULT(n)
      
      IMPLICIT NONE
      
      REAL(rp) :: n
      REAL(rp), INTENT(IN) :: a(3)
      
      n = sqrt(a(1)**2 + a(2)**2 + a(3)**2)
      
      END FUNCTION norm
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      FUNCTION phi(a,h) RESULT(p)
      
      IMPLICIT NONE
      
      REAL(rp) :: p
      REAL(rp), INTENT(IN) :: a(3),h
      
!       p = (1d0-(norm(a)**2)/h**2)**4
      p = exp(-(norm(a)**2)/h**2)
      
      END FUNCTION phi          
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      FUNCTION dphi(a,h) RESULT(gp)
      
      IMPLICIT NONE
      
      REAL(rp) :: gp(3)
      REAL(rp), INTENT(IN) :: a(3),h
      
!       gp = (-8d0*a*(1d0-norm(a)**2/h**2)**3)/h**2
      gp = ((-2d0*a)/h**2)*exp(-(norm(a)**2)/h**2)
      
      END FUNCTION dphi          

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      
      END MODULE evaluate