      MODULE evaluate
      
      USE globals, ONLY: rp,r,sigma_n
      USE lapack_interfaces
      
      IMPLICIT NONE
      
      INTEGER :: maxit,maxptit
      REAL(rp) :: sigma_r,threshold,percent_max      
      
      CONTAINS
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!             
      
      SUBROUTINE ref_elem_coords()
      
      USE globals, ONLY: np,mnnds,nel_type,rsre
      USE basis, ONLY: tri_nodes,quad_nodes
      
      IMPLICIT NONE
      INTEGER :: et,pt,n
      REAL(rp) :: r(mnnds),s(mnnds)
      
      
      DO et = 1,nel_type

        IF (mod(et,2) == 1) THEN
          CALL tri_nodes(1,np(1),n,r,s)
        ELSE IF (mod(et,2) == 0) THEN
          CALL quad_nodes(1,np(2),n,r,s)
        ENDIF
        
        DO pt = 1,n
          rsre(1,pt,et) = r(pt)
          rsre(2,pt,et) = s(pt)
        ENDDO
        
      ENDDO
      
      
      
      RETURN      
      END SUBROUTINE ref_elem_coords
 
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
      
      SUBROUTINE normals(mesh)
      
      USE globals, ONLY: np,nnds,mnnds,V,l,dlds,dldr,ipiv,grid     

      IMPLICIT NONE
      INTEGER :: it,el,et,p,n,i,np1,nnd
      REAL(rp) :: x,y
      REAL(rp) :: dxdr,dxds,dydr,dyds,jac
      REAL(rp) :: drdx,drdy,dsdx,dsdy
      REAL(rp) :: dhdx,dhdy
      REAL(rp) :: nx,ny,nz,nrm
      REAL(rp) :: x1,x2,y1,y2,z1,z2
      
      TYPE(grid) :: mesh
      
      PRINT "(A)", "Computing element normals..."      
            
      ALLOCATE(mesh%nhb(3,mesh%ne))
      ALLOCATE(mesh%xyhc(3,mesh%ne))
      
      mesh%xyhc = 0d0
      
      DO el = 1,mesh%ne
      
        et = mesh%el_type(el)
        p = np(et)  
        n = nnds(et)
          IF (mod(et,2) == 1) THEN
            nnd = nnds(3)
          ELSE IF (mod(et,2) == 0) THEN
            nnd = nnds(4)
          ENDIF             
        np1 = nnd+1             
        
        dxdr = 0d0
        dxds = 0d0
        dydr = 0d0
        dyds = 0d0  
        
        DO i = 1,n         
          dxdr = dxdr + dldr(i,np1,et)*mesh%elxy(i,el,1)
          dxds = dxds + dlds(i,np1,et)*mesh%elxy(i,el,1)
          dydr = dydr + dldr(i,np1,et)*mesh%elxy(i,el,2)
          dyds = dyds + dlds(i,np1,et)*mesh%elxy(i,el,2)
          
          mesh%xyhc(1,el) = mesh%xyhc(1,el) + l(i,np1,et)*mesh%elxy(i,el,1)
          mesh%xyhc(2,el) = mesh%xyhc(2,el) + l(i,np1,et)*mesh%elxy(i,el,2)
          mesh%xyhc(3,el) = mesh%xyhc(3,el) + l(i,np1,et)*mesh%elhb(i,el)        
        ENDDO
        
        jac = dxdr*dyds - dydr*dxds

        drdx =  dyds/jac
        drdy = -dxds/jac
        dsdx = -dydr/jac
        dsdy =  dxdr/jac
      
        dhdx = 0d0
        dhdy = 0d0      

        DO i = 1,n
          dhdx = dhdx + (dldr(i,np1,et)*drdx + dlds(i,np1,et)*dsdx)*mesh%elhb(i,el)
          dhdy = dhdy + (dldr(i,np1,et)*drdy + dlds(i,np1,et)*dsdy)*mesh%elhb(i,el)        
        ENDDO
      
        nrm = sqrt(dhdx**2 + dhdy**2 + 1d0)
      
        mesh%nhb(1,el) = dhdx/nrm
        mesh%nhb(2,el) = dhdy/nrm
        mesh%nhb(3,el) = -1d0/nrm
      
!       !Cross product normals as a check
!       x1 = elxy(2,el,1) - elxy(1,el,1)
!       x2 = elxy(3,el,1) - elxy(1,el,1)
!       
!       y1 = elxy(2,el,2) - elxy(1,el,2)
!       y2 = elxy(3,el,2) - elxy(1,el,2)
!       
!       z1 = elhb(2,el) - elhb(1,el)
!       z2 = elhb(3,el) - elhb(1,el)
!       
!       nx =   y1*z2 - y2*z1
!       ny = -(x1*z2 - x2*z1)
!       nz =   x1*y2 - x2*y1
!       
!       nrm = sqrt(nx**2 + ny**2 + nz**2)
!       
!       nx = nx/nrm
!       ny = ny/nrm
!       nz = nz/nrm
!       
!       PRINT*, (nhb(i,el), i=1,3)
!       PRINT*, nx,ny,nz
!       PRINT*, " " 
        
      ENDDO  
      
      RETURN 
      END SUBROUTINE normals
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE coordinates(mesh)
      
      USE globals, ONLY: rp,mnnds,nnds,np,nverts,l, &
                         grid
      
      IMPLICIT NONE
      
      INTEGER :: el,i,pt,led,j,ed
      INTEGER :: et,n,pn,nnd,nv,bed
      
      REAL(rp) :: xpt,ypt,ytest
      
      TYPE(grid) :: mesh
           
      
      ALLOCATE(mesh%xyhi(mnnds,mesh%ne,3))      
      ALLOCATE(mesh%xyhe(mnnds,mesh%ned,3))
      ALLOCATE(mesh%bnd_flag(mnnds,mesh%ned))

      
      mesh%bnd_flag(:,:) = 0
      
      mesh%xyhe = 0d0
      mesh%xyhi = 0d0
      
      PRINT "(A)", "Computing extra edge nodes..."        
      
      DO ed = 1,mesh%ned
      
        el = mesh%ged2el(1,ed)
        led = mesh%ged2led(1,ed)
        et = mesh%el_type(el)
        nv = nverts(et)
        n = nnds(et)
        bed = mesh%ed_type(ed)
        
        IF (mod(et,2) == 1) THEN   
          pn = np(3)
        ELSE IF (mod(et,2) == 0) THEN
          pn = np(4)
        ENDIF    
        
        DO pt = 1,pn-1
          
          j = mod(led,nv)*pn + pt + 1
          
          DO i = 1,n                   
            mesh%xyhe(pt,ed,1) = mesh%xyhe(pt,ed,1) + l(i,j,et)*mesh%elxy(i,el,1)
            mesh%xyhe(pt,ed,2) = mesh%xyhe(pt,ed,2) + l(i,j,et)*mesh%elxy(i,el,2)
            mesh%xyhe(pt,ed,3) = mesh%xyhe(pt,ed,3) + l(i,j,et)*mesh%elhb(i,el)
          ENDDO 
          
          IF (bed == 10) THEN
            ytest = mesh%xyhe(pt,ed,2)
            xpt = mesh%xyhe(pt,ed,1)
            
            IF (ytest < 250d0) THEN
              ypt = 0d0 + 100d0*(1d0/(COSH(4d0*(xpt-2000d0)/500d0)))
            ELSE IF (ytest > 250d0) THEN
              ypt = 500d0 - 100d0*(1d0/(COSH(4d0*(xpt-2000d0)/500d0)))
            ENDIF
          
            mesh%xyhe(pt,ed,2) = ypt
          ENDIF
          
        ENDDO                               
      
      ENDDO
      
      PRINT "(A)", "Computing extra interior element nodes..."          
      
      
      DO el = 1,mesh%ne
      
        et = mesh%el_type(el)
        n = nnds(et)
        nv = nverts(et)
        
        IF (mod(et,2) == 1) THEN   
          pn = np(3)
          nnd = nnds(3)
        ELSE IF (mod(et,2) == 0) THEN
          pn = np(4)
          nnd = nnds(4)          
        ENDIF    

        pt = 0
        DO j = nv*(pn-1)+nv+1,nnd
        
          pt = pt + 1        
          
          DO i = 1,n                   
            mesh%xyhi(pt,el,1) = mesh%xyhi(pt,el,1) + l(i,j,et)*mesh%elxy(i,el,1)
            mesh%xyhi(pt,el,2) = mesh%xyhi(pt,el,2) + l(i,j,et)*mesh%elxy(i,el,2)
            mesh%xyhi(pt,el,3) = mesh%xyhi(pt,el,3) + l(i,j,et)*mesh%elhb(i,el)
          ENDDO        
        ENDDO
                
      ENDDO
      
      RETURN
      END SUBROUTINE coordinates
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE transformation()
      
      USE globals, ONLY: mnnds,nel_type,nnds,np,nverts,l,dldr,dlds
      USE basis, ONLY: tri_nodes,quad_nodes
      USE shape_functions_mod
      
      IMPLICIT NONE
      
      INTEGER :: et,pt,i,j
      INTEGER :: n,p,np1,pn,nnd,nv
      INTEGER :: info
      REAL(rp), DIMENSION(mnnds+1) :: r,s
      
      PRINT "(A)", "Computing interpolating polynomials..."        
      
      ALLOCATE(l(mnnds,mnnds+1,nel_type))
      ALLOCATE(dldr(mnnds,mnnds+1,nel_type))
      ALLOCATE(dlds(mnnds,mnnds+1,nel_type))
      
      DO et = 1,nel_type
        n = nnds(et)
        p = np(et)
        nv = nverts(et)
      
        IF (mod(et,2) == 1) THEN    
          pn = np(3)
          CALL tri_nodes(1,pn,nnd,r,s)
        ELSE IF (mod(et,2) == 0) THEN
          pn = np(4)
          CALL quad_nodes(1,pn,nnd,r,s)
        ENDIF     
        
        np1 = nnd+1
        
        IF (mod(et,2) == 1) THEN
          r(np1) = -1d0/3d0
          s(np1) = -1d0/3d0
        ELSE IF (mod(et,2) == 0) THEN
          r(np1) = 0d0
          s(np1) = 0d0
        ENDIF        
        
        CALL shape_functions_area_eval(nv,p,n,np1,r,s,l(:,:,et),dldr(:,:,et),dlds(:,:,et))
        
      ENDDO
      
      
      RETURN
      END SUBROUTINE transformation
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE compute_surface()
      
      USE globals, ONLY: refinement,np,mninds,mnnds,base,eval,kdresults,closest,tree_xy,tree_c,hmin,hmax,srchdp
      USE kdtree2_module      
      
      IMPLICIT NONE
      
      srchdp = 10
      
      maxit = 100
      maxptit = 100
      threshold = 1d-12
      percent_max = 100d0
      sigma_r = 0.5d0 
!       sigma_n = 1.5d0   ! 0.5 - 1.5
!       r = 15d0           ! 1.5 - 4.0      
      
      hmin = minval(base%depth)
      hmax = maxval(base%depth)
      
      ! Build kd-tree           
!       tree_xy => kdtree2_create(vxy , rearrange=.true., sort=.true.)
      tree_xy => kdtree2_create(base%xyhc(1:2,:)  , rearrange=.true., sort=.true.)
      tree_c  => kdtree2_create(base%xyhc, rearrange=.true., sort=.true.)
      
      ALLOCATE(kdresults(base%ne))       
      ALLOCATE(closest(srchdp))        
!       
      CALL grid_size(base)
      
!       CALL filter_normals()
      
      
      IF (refinement) THEN
        PRINT("(A)"), "Computing rimls surface: verticies"
        CALL mls_surface(eval%nn,1,1,eval%xyhv)      
        PRINT("(A)"), "Computing rimls surface: edges"      
        CALL mls_surface(eval%ned,np(3)-1,mnnds,eval%xyhe)
        PRINT("(A)"), "Computing rimls surface: interior"
        CALL mls_surface(eval%ne,mninds,mnnds,eval%xyhi)      
      ELSE 
        PRINT("(A)"), "Computing rimls surface: verticies"
        CALL mls_surface(base%nn,1,1,base%xyhv)      
        PRINT("(A)"), "Computing rimls surface: edges"      
        CALL mls_surface(base%ned,np(3)-1,mnnds,base%xyhe)
        PRINT("(A)"), "Computing rimls surface: interior"
        CALL mls_surface(base%ne,mninds,mnnds,base%xyhi)         
      ENDIF
      
!       IF (refinement) THEN
!         PRINT("(A)"), "Computing rimls surface: verticies"
!         CALL rimls_surface(eval%nn,1,1,eval%xyhv)      
!         PRINT("(A)"), "Computing rimls surface: edges"      
!         CALL rimls_surface(eval%ned,np(3)-1,mnnds,eval%xyhe)
!         PRINT("(A)"), "Computing rimls surface: interior"
!         CALL rimls_surface(eval%ne,mninds,mnnds,eval%xyhi)      
!       ELSE 
!         PRINT("(A)"), "Computing rimls surface: verticies"
!         CALL rimls_surface(base%nn,1,1,base%xyhv)      
!         PRINT("(A)"), "Computing rimls surface: edges"      
!         CALL rimls_surface(base%ned,np(3)-1,mnnds,base%xyhe)
!         PRINT("(A)"), "Computing rimls surface: interior"
!         CALL rimls_surface(base%ne,mninds,mnnds,base%xyhi)         
!       ENDIF
      
!       IF (refinement) THEN
!         PRINT("(A)"), "Computing rimls surface: verticies"
!         CALL function_surface(eval%nn,1,1,eval%xyhv)      
!         PRINT("(A)"), "Computing rimls surface: edges"      
!         CALL function_surface(eval%ned,np(3)-1,mnnds,eval%xyhe)
!         PRINT("(A)"), "Computing rimls surface: interior"
!         CALL function_surface(eval%ne,mninds,mnnds,eval%xyhi)      
!       ELSE 
!         PRINT("(A)"), "Computing rimls surface: verticies"
!         CALL function_surface(base%nn,1,1,base%xyhv)      
!         PRINT("(A)"), "Computing rimls surface: edges"      
!         CALL function_surface(base%ned,np(3)-1,mnnds,base%xyhe)
!         PRINT("(A)"), "Computing rimls surface: interior"
!         CALL function_surface(base%ne,mninds,mnnds,base%xyhi)         
!       ENDIF      
!       
      RETURN
      END SUBROUTINE
      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE filter_normals()
      
      USE globals, ONLY: base,tree_c,kdresults
      USE kdtree2_module          
      
      IMPLICIT NONE
      
      INTEGER :: j,i,it
      INTEGER :: nneigh,el
      REAL(rp) :: xi(3),pj(3),ni(3),nj(3),top(3),bottom
      REAL(rp) :: hpt,phiw,phi0
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: nrm
      
      ALLOCATE(nrm(3,base%ne))
      
      DO j = 1,base%ne
      
        pj = base%xyhc(:,j)       
      
        CALL kdtree2_n_nearest(tp=tree_c,qv=pj,nn=1,results=kdresults)
          
        el = kdresults(1)%idx           
        hpt = (r*base%h(el))**2   
          
        CALL kdtree2_r_nearest(tp=tree_c,qv=pj,r2=hpt,nfound=nneigh,nalloc=base%ne,results=kdresults)
          
        hpt = sqrt(hpt)  
        
        top = 0d0
        bottom = 0d0
        DO i = 1,nneigh
          xi  = base%xyhc(:,kdresults(i)%idx)
          ni = base%nhb(:,kdresults(i)%idx)  
          
          phi0 = (1d0-norm(pj-xi)**2/hpt**2)**4
          
          top = top + phi0*ni
          bottom = bottom + phi0
          
        ENDDO       
        
        nrm(:,j) = top/bottom
        
  iter: DO it = 1,maxptit
          top = 0d0
          bottom = 0d0
          nj = nrm(:,j)
          DO i = 1,nneigh
          
            xi  = base%xyhc(:,kdresults(i)%idx)
            ni = base%nhb(:,kdresults(i)%idx)   
                        
            phiw = (1d0-norm(pj-xi)**2/hpt**2)**4*exp(-(norm(nj-ni)/sigma_n)**2)
            
            top = top + phiw*ni
            bottom = bottom + phiw
          ENDDO
          nrm(:,j) = top/bottom
          IF (norm(nrm(:,j)-nj) < threshold) THEN
            PRINT*, "iter = ", it
            EXIT iter
          ENDIF
        ENDDO iter
!         IF (it == maxptit+1) THEN
!           PRINT*, "max iterations reached"
!         ENDIF
      ENDDO

      base%nhb = nrm
      
      RETURN
      END SUBROUTINE filter_normals

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
      
      SUBROUTINE rimls_surface(n,npts,mnpts,xyh)
      
      USE globals, ONLY: tree_xy,tree_c,kdresults,base,hmin,hmax,phi0,lambda0,Erad
      USE find_element, ONLY: in_element
      USE kdtree2_module
      
      IMPLICIT NONE
      
      INTEGER :: i,j,pt,it,ptit,cpt
      INTEGER :: n,npts,mnpts
      INTEGER :: nneigh,nneigh_pre
      INTEGER :: nd,el
      INTEGER :: neighbors(base%ne),found
      REAL(rp) :: search_r
      REAL(rp) :: xyh(mnpts,n,3)
      REAL(rp) :: x(3),p(3),px(3),np(3),hpt,f,grad_f(3),grad_w(3),fgradf(3),npmgradf(3)
      REAL(rp) :: alpha(base%ne),alpha_old(base%ne)
      REAL(rp) :: sumW,sumF
      REAL(rp) :: sumGw(3),sumGF(3),sumN(3)
      REAL(rp) :: w,fx
     
      

      
      DO i = 1,n
       
        IF (mod(i,1000) == 0) THEN       
          PRINT*,i,"/",n
        ENDIF
        DO pt = 1,npts  
        
          x(1) = xyh(pt,i,1)
          x(2) = xyh(pt,i,2)
          x(3) = xyh(pt,i,3)          
      
!           CALL kdtree2_n_nearest(tp=tree_xy,qv=(/x(1),x(2)/),nn=1,results=kdresults)
!           
! !           nd = vxyn(kdresults(1)%idx)
!            nd = kdresults(1)%idx           
!            hpt = (r*base%h(base%epn(1,nd)))**2   

!           CALL kdtree2_n_nearest(tp=tree_xy,qv=x(1:2),nn=1,results=kdresults)
!           el = kdresults(1)%idx        
          
          CALL in_element(i,x(1:2),el,found)                    
   
          hpt = r*base%h(el) 
          search_r = SQRT(-hpt**2*LOG(threshold))          
          CALL kdtree2_r_nearest(tp=tree_xy,qv=x(1:2),r2=search_r**2,nfound=nneigh,nalloc=base%ne,results=kdresults)
          
!           IF (nneigh == 0) THEN
!             PRINT*, x(1)/(Erad*cos(phi0))+lambda0,x(2)/Erad
!             PRINT*,el, base%ged2el(1,i),base%ged2el(2,i)
!           ENDIF
          
          
!           IF (i == 8595) PRINT*, el, nneigh
          
!           CALL boundary_check(i,x,nneigh,neighbors)
          CALL boundary_check2(i,el,nneigh,neighbors)
                 
          
          grad_f(:) = 1d0
          f = 1d0          

          ptit = 0
          fgradf = f*grad_f
   wloop: DO WHILE (norm(fgradf) > threshold)
            it = 0
                        
            DO j = 1,nneigh
              alpha(j) = 1d0
              alpha_old(j) = 0d0
            ENDDO
            
            DO WHILE (it < maxit)
              
              sumW = 0d0
              sumF = 0d0
              sumGw(:) = 0d0
              sumGF(:) = 0d0
              sumN(:) = 0d0
            
              DO cpt = 1,nneigh
              
!                 el = kdresults(cpt)%idx
                el = neighbors(cpt)
                p  = base%xyhc(:,el)
                np = base%nhb(:,el)
!                 hpt = r*base%h(el)   ! try using grid spacing of neighbor point
                
                px = x - p
                
                fx = DOT_PRODUCT(px,np)
                
                IF (it > 0) THEN
                  npmgradf = np-grad_f
                  alpha(cpt) = exp(-((fx-f)/(sigma_r*hpt))**2)*exp(-(norm(npmgradf)/sigma_n)**2)
                ELSE
                  alpha(cpt) = 1d0
                ENDIF                           
                
                w = alpha(cpt)*phi(px,hpt)
                grad_w = alpha(cpt)*dphi(px,hpt)               
                
                sumW = sumW + w
                sumGw = sumGw + grad_w
                sumF = sumF + w*fx
                sumGF = sumGF + grad_w*fx
                sumN = sumN + w*np                               
                
              ENDDO
              
              f = sumF/sumW
              grad_f = (sumGF - f*sumGw + sumN)/sumW
                            
              it = it + 1              
              
              
              DO j = 1,nneigh
                alpha_old(j) = alpha(j)
              ENDDO
              
            ENDDO
            
            fgradf = f*grad_f
            x = x - fgradf
            ptit = ptit + 1
            
            IF (ptit == maxptit) THEN
              EXIT wloop
            ENDIF
            
          ENDDO wloop
                
          IF (ptit < maxptit) THEN
          
            IF (x(1) /= x(1) .or. x(2) /= x(2) .or. x(3) /= x(3)) THEN
          
              PRINT*, "WARNING: Nan detected, i = ",i
              PRINT*, "    # neighbors = ", nneigh
                         
              DO cpt = 1,nneigh
                PRINT*, "    ", neighbors(cpt)
              ENDDO
              
            ELSE IF (x(3) < hmin) THEN 
            
              xyh(pt,i,3) = hmin 
!               PRINT*, "WARNING: Limiting hmin"
              
            ELSE IF (x(3) > hmax) THEN
            
              xyh(pt,i,3) = hmax
!               PRINT*, "WARNING: Limiting hmax"              
              
!             ELSE IF (abs((x(3)-xyh(pt,i,3))/xyh(pt,i,3))*100d0 > percent_max) THEN
!             
!               PRINT*, "WARNING: Percent change exceedes tolerance, i = ",i
              
            ELSE
            
!               print*, norm(fgradf)
              xyh(pt,i,1) = x(1)
              xyh(pt,i,2) = x(2)
              xyh(pt,i,3) = x(3)               

            ENDIF
            
          ELSE
            PRINT*, "WARNING: max point iterations reached, i = ",i
            STOP
          ENDIF
      
        ENDDO
      ENDDO
      
      CALL invcpp(n,npts,mnpts,xyh)
      
      
      RETURN
      END SUBROUTINE rimls_surface        
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
      
      SUBROUTINE mls_surface(n,npts,mnpts,xyh)
      
      USE globals, ONLY: tree_xy,tree_c,kdresults,base,hmin,hmax,phi0,lambda0,Erad,lsp
      USE find_element, ONLY: in_element,newton
      USE kdtree2_module
      USE basis, ONLY: tri_basis,quad_basis
      
      IMPLICIT NONE
      
      INTEGER :: i,pt,cpt,l,m,k      
      INTEGER :: n,npts,mnpts  
      INTEGER :: ne,el,elin,et
      INTEGER :: ndf
      INTEGER :: neighbors(base%ne),nneigh,found
      INTEGER :: small_flag
      INTEGER :: info,lwork
      REAL(rp) :: xyh(mnpts,n,3)  
      REAL(rp) :: x(3),p(3),px(3),np(3),xbar(3)      
      REAL(rp) :: d,w,rhs,lhs
      REAL(rp) :: s,t,sv(1),tv(1)
      REAL(rp) :: hpt,search_r,tol
      REAL(rp) :: nwork(1)
      REAL(rp), ALLOCATABLE :: A(:,:),b(:),work(:)
      REAL(rp), ALLOCATABLE :: phi(:,:)
      
      ne = base%ne
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
          

          CALL in_element(i,x(1:2),elin,found)  
          et = base%el_type(elin)          
   
          hpt = r*base%h(elin)           
          search_r = SQRT(-hpt**2*LOG(tol))
          CALL kdtree2_r_nearest(tp=tree_xy,qv=x(1:2),r2=search_r**2,nfound=nneigh,nalloc=base%ne,results=kdresults)
          
          CALL boundary_check2(i,elin,nneigh,neighbors)
          
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!   Constant Least Squares Fit !!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           rhs = 0d0
!           lhs = 0d0
!           small_flag = 0
!           DO cpt = 1,nneigh
!           
!             el = neighbors(cpt)
!             p  = base%xyhc(:,el)
!             np = base%nhb(:,el)  
!             
!             w = theta(x,p,hpt)
!             
!             IF (w < 1d-10) THEN
!               small_flag = 1
!             ENDIF
!             
!             rhs = rhs + w**2*p(3)
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
!           xyh(pt,i,3) = rhs/lhs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


          xbar = base%xyhc(:,elin)

          small_flag = 0
          DO cpt = 1,nneigh
          
            el = neighbors(cpt)
            p  = base%xyhc(:,el)
            np = base%nhb(:,el)  
            
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
      
      CALL invcpp(n,npts,mnpts,xyh)      
      
      RETURN
      END SUBROUTINE mls_surface

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

      SUBROUTINE grid_size(mesh)
      
      USE globals, ONLY: grid
      
      IMPLICIT NONE
      
      INTEGER :: el
      REAL(rp) :: x1,x2,x3
      REAL(rp) :: y1,y2,y3
      REAL(rp) :: a,b,c,s,r
      
      TYPE(grid) :: mesh
      
      ALLOCATE(mesh%h(mesh%ne))
                 
      DO el = 1,mesh%ne
        x1 = mesh%xy(1,mesh%ect(1,el))
        x2 = mesh%xy(1,mesh%ect(2,el))
        x3 = mesh%xy(1,mesh%ect(3,el))
        
        y1 = mesh%xy(2,mesh%ect(1,el))
        y2 = mesh%xy(2,mesh%ect(2,el))
        y3 = mesh%xy(2,mesh%ect(3,el))
        
        a = sqrt((x1-x2)**2+(y1-y2)**2)
        b = sqrt((x2-x3)**2+(y2-y3)**2)
        c = sqrt((x3-x1)**2+(y3-y1)**2)
        
        s = .5d0*(a+b+c)
        r = sqrt((s-a)*(s-b)*(s-c)/s)
        
        mesh%h(el) = 2d0*r
      ENDDO
      
      RETURN
      END SUBROUTINE grid_size
      
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

      SUBROUTINE invcpp(n,np,mnp,xyh,xyh2)
      
      USE globals, ONLY: Erad,phi0,lambda0
      
      IMPLICIT NONE
      
      INTEGER :: i,pt
      INTEGER :: n,np,mnp
      REAL(rp) :: xyh(mnp,n,3)
      REAL(rp), OPTIONAL :: xyh2(mnp,n,3)      
      
      IF(PRESENT(xyh2)) THEN
      
        DO i = 1,n
          DO pt = 1,np
            xyh2(pt,i,1) = xyh(pt,i,1)/(Erad*cos(phi0))+lambda0
            xyh2(pt,i,2) = xyh(pt,i,2)/Erad
            xyh2(pt,i,3) = xyh(pt,i,3)           
          ENDDO
        ENDDO      
      
      ELSE
      
        DO i = 1,n
          DO pt = 1,np
            xyh(pt,i,1) = xyh(pt,i,1)/(Erad*cos(phi0))+lambda0
            xyh(pt,i,2) = xyh(pt,i,2)/Erad
            xyh(pt,i,3) = xyh(pt,i,3)           
          ENDDO
        ENDDO
        
      ENDIF
      
      RETURN
      END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      
      END MODULE evaluate