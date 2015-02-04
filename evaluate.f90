      MODULE evaluate
      
      USE globals, ONLY: pres
      
      IMPLICIT NONE
      
      INTEGER :: maxit,maxptit
      REAL(pres) :: r,sigma_n,sigma_r,threshold      
      
      CONTAINS
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!             
      
      SUBROUTINE vandermonde()
      
      USE globals, ONLY: ctp,np,nnds,mnnds,nel_type,V,ipiv
      USE basis, ONLY: tri_nodes,tri_basis,quad_nodes,quad_basis
      
      IMPLICIT NONE
      INTEGER :: et,pt,dof,i,n
      REAL(pres) :: r(mnnds),s(mnnds)
      REAL(pres) :: phi(mnnds*mnnds)
      INTEGER :: info  
      
      PRINT "(A)", "Computing Vandermode matrix..."
      
      ALLOCATE(V(mnnds,mnnds,nel_type))
      ALLOCATE(ipiv(mnnds,nel_type))
      
      DO et = 1,nel_type
        n = nnds(et)
        IF (mod(et,2) == 1) THEN
          CALL tri_nodes(1,np(et),n,r,s)
          CALL tri_basis(np(et),n,n,r,s,phi)       
        ELSE IF (mod(et,2) == 0) THEN
          CALL quad_nodes(1,np(et),n,r,s)
          CALL quad_basis(np(et),n,n,r,s,phi)
        ENDIF
        
        DO pt = 1,n
          DO dof = 1,n
            i = (dof-1)*n + pt
            V(dof,pt,et) = phi(i)
          ENDDO
        ENDDO
        
        CALL DGETRF(n,n,V(1,1,et),mnnds,ipiv(1,et),info)        
!         DO pt = 1,n
!             PRINT("(100(e15.5))"), (V(dof,pt,et), dof = 1,n)
!         ENDDO        
!         PRINT*, " "
        
      ENDDO
      
      
      
      RETURN      
      END SUBROUTINE vandermonde
 
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
      
      SUBROUTINE normals(mesh)
      
      USE globals, ONLY: np,nnds,mnnds,V,l,ipiv,grid     

      IMPLICIT NONE
      INTEGER :: it,el,et,p,n,i,np1,nnd
      REAL(pres) :: x,y
      REAL(pres) :: dxdr,dxds,dydr,dyds,jac
      REAL(pres) :: drdx,drdy,dsdx,dsdy
      REAL(pres) :: dhdx,dhdy
      REAL(pres) :: nx,ny,nz,nrm
      REAL(pres) :: x1,x2,y1,y2,z1,z2
      
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
          dxdr = dxdr + l(i,2*np1,et)*mesh%elxy(i,el,1)
          dxds = dxds + l(i,3*np1,et)*mesh%elxy(i,el,1)
          dydr = dydr + l(i,2*np1,et)*mesh%elxy(i,el,2)
          dyds = dyds + l(i,3*np1,et)*mesh%elxy(i,el,2)
          
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
          dhdx = dhdx + (l(i,2*np1,et)*drdx + l(i,3*np1,et)*dsdx)*mesh%elhb(i,el)
          dhdy = dhdy + (l(i,2*np1,et)*drdy + l(i,3*np1,et)*dsdy)*mesh%elhb(i,el)        
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
      
      USE globals, ONLY: mnnds,nnds,np,nverts,l, &
                         grid
      
      IMPLICIT NONE
      
      INTEGER :: el,i,pt,led,j,ed
      INTEGER :: et,n,pn,nnd,nv,bed
      
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
        ENDDO
        
        bed = mesh%bed_flag(ed)
        
        IF (bed == 1) THEN
          DO pt = 1,pn-1
            mesh%bnd_flag(pt,ed) = 1
          ENDDO
        ENDIF
        
      
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
      
      USE globals, ONLY: mnnds,nel_type,nnds,np,l,V,ipiv
      USE basis, ONLY: tri_nodes,quad_nodes,tri_basis,quad_basis
      
      IMPLICIT NONE
      
      INTEGER :: et,pt,i,j
      INTEGER :: n,p,np1,pn,nnd
      INTEGER :: info
      REAL(pres), DIMENSION(mnnds+1,nel_type) :: r,s
      REAL(pres), DIMENSION(mnnds*(mnnds+1)) :: phi,dpdr,dpds
      
      PRINT "(A)", "Computing interpolating polynomials..."        
      
      ALLOCATE(l(mnnds,3*(mnnds+1),nel_type))
      
      DO et = 1,nel_type
        n = nnds(et)
        p = np(et)
      
        IF (mod(et,2) == 1) THEN
          nnd = nnds(3)      
          pn = np(3)
          CALL tri_nodes(1,pn,nnd,r(1,et),s(1,et))
        ELSE IF (mod(et,2) == 0) THEN
          nnd = nnds(4) 
          pn = np(4)
          CALL quad_nodes(1,pn,nnd,r(1,et),s(1,et))
        ENDIF     
        
        np1 = nnd+1
        
        IF (mod(et,2) == 1) THEN
          r(np1,et) = -1d0/3d0
          s(np1,et) = -1d0/3d0
        ELSE IF (mod(et,2) == 0) THEN
          r(np1,et) = 1d0
          s(np1,et) = 1d0
        ENDIF        
        
        IF (mod(et,2) == 1) THEN
          CALL tri_basis(p,n,np1,r(1,et),s(1,et),phi,dpdr,dpds)
        ELSE IF (mod(et,2) == 0) THEN
          CALL quad_basis(p,n,np1,r(1,et),s(1,et),phi,dpdr,dpds)
        ENDIF     
        
!         DO i = 1,np1
!           PRINT*, r(i,et),s(i,et)
!         ENDDO
!         PRINT*, "" 
        
        DO pt = 1,np1
          DO i = 1,n
            j = (i-1)*np1+pt
          
            l(i,pt,et) = phi(j)
            l(i,np1+pt,et) = dpdr(j)
            l(i,2*np1+pt,et) = dpds(j)          
          ENDDO         
        ENDDO

        CALL DGETRS("N",n,3*np1,V(1,1,et),mnnds,ipiv(1,et),l(1,1,et),mnnds,info)
      ENDDO
      
      
      RETURN
      END SUBROUTINE transformation
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE compute_surface()
      
      USE globals, ONLY: refinement,np,mninds,mnnds,base,fine,kdresults,tree_xy,tree_c
      USE kdtree2_module      
      
      IMPLICIT NONE
      
      maxit = 10
      maxptit = 100
      threshold = 1d-4
      sigma_r = 0.5d0 
      sigma_n = 0.5d0   ! 0.5 - 1.5
      r = 2.5d0           ! 1.5 - 4.0      
      
      ! Build kd-tree           
!       tree_xy => kdtree2_create(vxy , rearrange=.true., sort=.true.)
      tree_xy => kdtree2_create(base%xy  , rearrange=.true., sort=.true.)
      tree_c  => kdtree2_create(base%xyhc, rearrange=.true., sort=.true.)
      
      ALLOCATE(kdresults(base%ne))       
!       
      CALL grid_size(base)
      
!       CALL filter_normals()
      
      
      IF (refinement) THEN
        PRINT("(A)"), "Computing rimls surface: verticies"
        CALL rimls_surface(fine%nn,1,1,fine%xyhv)      
        PRINT("(A)"), "Computing rimls surface: edges"      
        CALL rimls_surface(fine%ned,np(3)-1,mnnds,fine%xyhe)
        PRINT("(A)"), "Computing rimls surface: interior"
        CALL rimls_surface(fine%ne,mninds,mnnds,fine%xyhi)      
      ELSE 
        PRINT("(A)"), "Computing rimls surface: verticies"
        CALL rimls_surface(base%nn,1,1,base%xyhv)      
        PRINT("(A)"), "Computing rimls surface: edges"      
        CALL rimls_surface(base%ned,np(3)-1,mnnds,base%xyhe)
        PRINT("(A)"), "Computing rimls surface: interior"
        CALL rimls_surface(base%ne,mninds,mnnds,base%xyhi)         
      ENDIF
      
      RETURN
      END SUBROUTINE
      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!       SUBROUTINE filter_normals()
!       
!       USE globals, ONLY: base,tree_xy,tree_c,kdresults
!       
!       IMPLICIT NONE
!       
!       INTEGER :: j,i,it
!       INTEGER :: nneigh
!       REAL(pres) :: pj(3),top(3),bottom
!       REAL(pres) :: hpt,phiw
!       REAL(pres), ALLOCATABLE, DIMENSION(:,:) :: nrm
!       
!       ALLOCATE(nrm(3,base%ne))
!       
!       DO j = 1,base%ne
!         x(1) = base%xyc(1,j)
!         x(2) = base%xyc(2,j)
!         x(3) = base%xyc(3,j)          
!       
!         CALL kdtree2_n_nearest(tp=tree_xy,qv=(/x(1),x(2)/),nn=1,results=kdresults)
!           
! !       nd = vxyn(kdresults(1)%idx)
!         nd = kdresults(1)%idx           
!         hpt = (r*base%h(base%epn(1,nd)))**2   
!           
!         CALL kdtree2_r_nearest(tp=tree_c,qv=x,r2=hpt,nfound=nneigh,nalloc=base%ne,results=kdresults)
!           
!         hpt = sqrt(hpt)    
!         
!         DO it = 1,maxptit
!           top = 0d0
!           bottom = 0d0
!           DO i = 1,nneigh
!           
!             p  = base%xyhc(:,kdresults(i)%idx)
!             np = base%nhb(:,kdresults(i)%idx)          
!             top = top + phi2*nrm(j)
!           ENDDO
!           nrm = top/bottom
!         ENDDO
!       ENDDO
! 
!       
!       
!       RETURN
!       END SUBROUTINE filter_normals

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
      
      SUBROUTINE rimls_surface(n,npts,mnpts,xyh)
      
      USE globals, ONLY: tree_xy,tree_c,kdresults,base
      USE kdtree2_module
      
      IMPLICIT NONE
      
      INTEGER :: i,j,pt,it,ptit,cpt
      INTEGER :: n,npts,mnpts
      INTEGER :: nneigh
      INTEGER :: nd,el
      REAL(pres) :: xyh(mnpts,n,3)
      REAL(pres) :: x(3),p(3),px(3),np(3),hpt,f,grad_f(3),grad_w(3),fgradf(3),npmgradf(3)
      REAL(pres) :: alpha(base%ne),alpha_old(base%ne)
      REAL(pres) :: sumW,sumF
      REAL(pres) :: sumGw(3),sumGF(3),sumN(3)
      REAL(pres) :: w,fx
     
      

      
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

          CALL kdtree2_n_nearest(tp=tree_c,qv=x,nn=1,results=kdresults)
          
          el = kdresults(1)%idx           
          hpt = (r*base%h(el))**2   
          
          CALL kdtree2_r_nearest(tp=tree_c,qv=x,r2=hpt,nfound=nneigh,nalloc=base%ne,results=kdresults)
          
          hpt = sqrt(hpt)           
          
          grad_f(:) = 1d0
          f = 1d0          

          ptit = 0
          fgradf = f*grad_f
          DO WHILE (norm(fgradf) > threshold .and. ptit < maxptit)
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
              
                el = kdresults(cpt)%idx
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
                
                w = alpha(cpt)*(1d0-(norm(px)**2)/hpt**2)**4
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
            
          ENDDO
                
          IF (ptit < maxptit) THEN
          
            xyh(pt,i,1) = x(1)
            xyh(pt,i,2) = x(2)
            xyh(pt,i,3) = x(3)    
            
          ELSE
            PRINT*, "WARNING: max point iterations reached, i=:",i
            STOP
          ENDIF
      
        ENDDO
      ENDDO
      
      CALL invcpp(n,npts,mnpts,xyh)
      
      
      RETURN
      END SUBROUTINE rimls_surface      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE grid_size(mesh)
      
      USE globals, ONLY: grid
      
      IMPLICIT NONE
      
      INTEGER :: el
      REAL(pres) :: x1,x2,x3
      REAL(pres) :: y1,y2,y3
      REAL(pres) :: a,b,c,s,r
      
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
      
      REAL(pres) :: n
      REAL(pres), INTENT(IN) :: a(3)
      
      n = sqrt(a(1)**2 + a(2)**2 + a(3)**2)
      
      END FUNCTION norm
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      FUNCTION dphi(a,h) RESULT(gp)
      
      IMPLICIT NONE
      
      REAL(pres) :: gp(3)
      REAL(pres), INTENT(IN) :: a(3),h
      
      gp = (-8d0*a*(1d0-norm(a)**2/h**2)**3)/h**2
      
      END FUNCTION dphi    
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE invcpp(n,np,mnp,xyh,xyh2)
      
      USE globals, ONLY: Erad,phi0,lambda0
      
      IMPLICIT NONE
      
      INTEGER :: i,pt
      INTEGER :: n,np,mnp
      REAL(pres) :: xyh(mnp,n,3)
      REAL(pres), OPTIONAL :: xyh2(mnp,n,3)      
      
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