      SUBROUTINE compute_random()

      USE globals, ONLY: rp,pi,nrpt,eps,xy_rand,h_rand,eval,np,mnnds,mninds,out_direc,base, &
                         tree_xy,tree_xy_rand,srchdp,closest,kdresults
      USE evaluate, ONLY: grid_size                         
      USE kdtree2_module

      IMPLICIT NONE
      
      INTEGER :: i,j
      INTEGER :: pt
      REAL(rp) :: random,x,y,h
      REAL(rp) :: davg
      REAL(rp) :: f1,f2  
      REAL(rp) :: pm

     

!       nrpt = 40000
!       eps = 0d0
      
      
      PRINT("(A)"), "Computing random points"
!       PRINT*, "  nrpt", nrpt
!       PRINT*, "  eps", eps      
!       PRINT*, " "
      
      ALLOCATE(xy_rand(3,nrpt))
      
      CALL RANDOM_SEED(put=(/3/))
      
      DO i = 1,nrpt
      
        CALL RANDOM_NUMBER(random)        
        x = 6000d0*random   
        
        CALL RANDOM_NUMBER(random)        
        y = (f2(x)-f1(x))*random + f1(x)
        
        CALL RANDOM_NUMBER(random)
        IF (random >= .5d0) THEN
          pm = 1d0
        ELSE 
          pm = -1d0
        ENDIF
        
        CALL RANDOM_NUMBER(random)       
        h = pm*eps*random
        
        xy_rand(1,i) = x
        xy_rand(2,i) = y
    
        y = 500d0/(f2(x)-f1(x))*y - 500d0/(f2(x)-f1(x))*f1(x) ! modified coordinate to adjust bathymetry profile 
                                                              ! to fit the constriction region
        xy_rand(3,i) = 10d0 - 5d0*cos(2d0*pi/500d0*y) + h
      
      
      ENDDO
           
           
   
   
      tree_xy_rand => kdtree2_create(xy_rand(1:2,:), rearrange=.true., sort=.true.)
      tree_xy => kdtree2_create(base%xyhc(1:2,:)  , rearrange=.true., sort=.true.)      
      
      srchdp = 10      
      ALLOCATE(kdresults(nrpt))       
      ALLOCATE(closest(srchdp))      
      
   
   
      ALLOCATE(h_rand(nrpt))
      DO i = 1,nrpt
        CALL kdtree2_n_nearest(tp=tree_xy_rand,qv=xy_rand(1:2,i),nn=srchdp,results=closest)
        
        davg = 0d0
        DO j = 1,srchdp
          pt = closest(j)%idx          
          davg = davg + sqrt((xy_rand(1,i)-xy_rand(1,j))**2 + (xy_rand(2,i)-xy_rand(2,j))**2)
        ENDDO        
        davg = davg/real(srchdp,rp)
        
        h_rand(i) = davg        
      ENDDO
      
      
      OPEN(unit=12,file=TRIM(out_direc) // 'random_pts.d')
      WRITE(12,*) nrpt
      DO pt = 1,nrpt
        WRITE(12,*) (xy_rand(i,pt),i=1,3)
      ENDDO
      
      CLOSE(12)      
      
      
      
      
      PRINT("(A)"), "Computing rimls surface: verticies"
      CALL mls_random(eval%nn,1,1,eval%xyhv)      
      PRINT("(A)"), "Computing rimls surface: edges"      
      CALL mls_random(eval%ned,np(3)-1,mnnds,eval%xyhe)
      PRINT("(A)"), "Computing rimls surface: interior"
      CALL mls_random(eval%ne,mninds,mnnds,eval%xyhi)      
       
     
      RETURN
      END SUBROUTINE compute_random
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
      
      SUBROUTINE mls_random(n,npts,mnpts,xyh)
      
      USE globals, ONLY: rp,nrpt,xy_rand,h_rand,tree_xy,tree_xy_rand,closest,kdresults,r,base,lsp
      USE kdtree2_module
      USE find_element, ONLY: in_element,newton
      USE lapack_interfaces
      USE basis, ONLY: tri_basis,quad_basis
      
      IMPLICIT NONE
      
      INTEGER :: i,pt,cpt,l,m,k,rpt      
      INTEGER :: n,npts,mnpts  
      INTEGER :: ne,el,elin,et
      INTEGER :: ndf
      INTEGER :: nneigh,found
      INTEGER :: small_flag
      INTEGER :: info,lwork
      REAL(rp) :: xyh(mnpts,n,3)  
      REAL(rp) :: x(3),p(3),px(3),np(3),xbar(3)      
      REAL(rp) :: d,w,rhs,lhs
      REAL(rp) :: s,t,sv(1),tv(1)
      REAL(rp) :: hpt,search_r,tol
      REAL(rp) :: nwork(1)
      REAL(rp) :: theta
      REAL(rp), ALLOCATABLE :: A(:,:),b(:),work(:)
      REAL(rp), ALLOCATABLE :: phi(:,:)      
      
      ne = nrpt
      tol = 1d-16
      ndf = (lsp+1)**2
      
      ALLOCATE(A(ne,ndf),b(ne))
      CALL DGELS('N',ne,ndf,1,A,ne,b,ne,nwork,-1,info)  ! find the optimal work array size
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
        
        DO pt = 1,npts  
        
          x(1) = xyh(pt,i,1)
          x(2) = xyh(pt,i,2)
          x(3) = xyh(pt,i,3)       
          

!           CALL kdtree2_n_nearest(tp=tree_xy,qv=x(1:2),nn=1,results=closest)   
!           hpt = r*h_rand(closest(1)%idx)  
!           xbar = x
          
          CALL in_element(i,x(1:2),elin,found) 
          et = base%el_type(elin) 
          
          hpt = r*base%h(elin) 
          xbar = base%xyhc(:,elin)
                  
          search_r = SQRT(-hpt**2*LOG(tol))
          CALL kdtree2_r_nearest(tp=tree_xy_rand,qv=x(1:2),r2=search_r**2,nfound=nneigh,nalloc=nrpt,results=kdresults)         




          small_flag = 0
          DO cpt = 1,nneigh
          
            rpt = kdresults(cpt)%idx
            p  = xy_rand(:,rpt)
            
            w = theta(x,p,hpt)
            
            IF (w <= tol*10d0) THEN
              small_flag = 1
            ENDIF
            
            CALL newton(cpt,p(1),p(2),elin,sv,tv)
            
            IF (mod(et,2) == 1) THEN
              CALL tri_basis(lsp,ndf,1,sv,tv,phi)
            ELSE IF (mod(et,2) == 0) THEN
              CALL quad_basis(lsp,ndf,1,sv,tv,phi)
            ENDIF
            
            DO m = 1,ndf
              A(cpt,m) = w*phi(m,1)
            ENDDO
            
            b(cpt) = w*p(3)               
            
!!!!!!!  Standard (non-orthogonal) Basis !!!!!!!                 
!             s = p(1)-xbar(1)
!             t = p(2)-xbar(2)
!             
!             k = 0
!             DO l = 0,lsp
!               DO m = 0,lsp-l
!                 k = k + 1            
!                 A(cpt,k) = w*(s**l*t**m)
!               ENDDO
!             ENDDO
!             
!             b(cpt) = w*p(3)

          ENDDO
          
          
!           IF (small_flag /= 1) THEN
!             PRINT*, "Increase search radius"
!             STOP
!           ENDIF
          
          CALL DGELS('N',nneigh,ndf,1,A,ne,b,ne,work,lwork,info)
          
          CALL newton(cpt,x(1),x(2),elin,sv,tv)

          IF (mod(et,2) == 1) THEN
            CALL tri_basis(lsp,ndf,1,sv,tv,phi)
          ELSE IF (mod(et,2) == 0) THEN
            CALL quad_basis(lsp,ndf,1,sv,tv,phi)
          ENDIF   
          
          xyh(pt,i,3) = 0d0
          DO m = 1,ndf
              xyh(pt,i,3) = xyh(pt,i,3) + b(m)*phi(m,1)
          ENDDO             
          
!!!!!!!  Standard (non-orthogonal) Basis !!!!!!!               
!           s = (x(1)-xbar(1))
!           t = (x(2)-xbar(2))
!           
!           k = 0
!           xyh(pt,i,3) = 0d0
!           DO l = 0,lsp
!             DO m = 0,lsp-l
!               k = k + 1
!               xyh(pt,i,3) = xyh(pt,i,3) + b(k)*(s**l*t**m)
!             ENDDO
!           ENDDO
          
          
        ENDDO
      
      ENDDO
      
!       CALL invcpp(n,npts,mnpts,xyh)      
      
      RETURN
      END SUBROUTINE mls_random      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
      FUNCTION f1(x) RESULT(y)
      
      USE globals, ONLY: rp
      
      IMPLICIT NONE
      
      REAL(rp) :: x
      REAL(rp) :: y
      
      y = 0d0 + 100d0*(1d0/(COSH(4d0*(x-2000d0)/500d0)))
      
      END FUNCTION
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      

      FUNCTION f2(x) RESULT(y)
      
      USE globals, ONLY: rp
      
      IMPLICIT NONE
      
      REAL(rp) :: x
      REAL(rp) :: y
      
      y = 500d0 - 100d0*(1d0/(COSH(4d0*(x-2000d0)/500d0)))
      
      END FUNCTION
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
      FUNCTION theta(x,p,h) RESULT(w)
      
      USE globals, ONLY: rp
      
      IMPLICIT NONE
      
      REAL(rp) :: w
      REAL(rp), INTENT(IN) :: x(3),p(3),h
      REAL(rp) :: d
      
      d = sqrt((x(1)-p(1))**2 + (x(2)-p(2))**2)
      w = exp(-d**2/h**2)      
      
      END FUNCTION theta              
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
