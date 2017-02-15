      MODULE calc_spline
      
      USE globals, ONLY: rp      
      
      IMPLICIT NONE
      
      
      
      CONTAINS
                  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE spline_init(num,nmax)
      
      USE globals, ONLY: rp,base,eval,ax,bx,cx,dx,ay,by,cy,dy,dt, &
                         nfbnds,fbnds,fbnds_xy,nfbnd2el,fbnd2el, &
                         nverts         
      
      IMPLICIT NONE
      
      INTEGER :: num
      INTEGER :: nbou
      INTEGER :: nmax
      INTEGER :: nd,bou,i,j,skip
      INTEGER :: el,et,nv,ged,led,edt
      INTEGER :: btype

      
      ! Count no normal flow boundaries and nodes
      
      num = 0       ! number of no normal flow boundaries
      nmax = 0      ! max number of nodes in any no normal flow boundary
      nfbnds = 0    ! number of total no normal flow boundary nodes
      DO bou = 1,base%nbou
        btype = base%fbseg(2,bou)
        IF( btype == 0 .OR. btype == 10 .OR. btype == 20  .OR. &   ! land boundaries
            btype == 1 .OR. btype == 11 .OR. btype == 21 ) THEN    ! island boundaries
        
          num = num + 1
          
          IF (base%fbseg(1,bou) > nmax) THEN
            nmax = base%fbseg(1,bou)
          ENDIF
          
          DO j = 1,base%fbseg(1,bou)
            nfbnds = nfbnds + 1            
          ENDDO
          
          
        ENDIF
      ENDDO
      
      ! Create a list of all no normal flow boundary nodes and their coordinates
      ! used to create the k-d tree to find spline coefficients to evaluate eval grid points
      
      ALLOCATE(fbnds(nfbnds),fbnds_xy(2,nfbnds))
      
      nfbnds = 0
      DO bou = 1,base%nbou
        btype = base%fbseg(2,bou)
        IF( btype == 0 .OR. btype == 10 .OR. btype == 20  .OR. &   ! land boundaries
            btype == 1 .OR. btype == 11 .OR. btype == 21 ) THEN    ! island boundaries
                            
          DO j = 1,base%fbseg(1,bou)
            nd = base%fbnds(j,bou)
            
            skip = 0           ! skip duplicate nodes
            DO i = 1,nfbnds
              IF (fbnds(i) == nd) THEN
                skip = 1
              ENDIF
            ENDDO
          
            IF (skip == 0) THEN
              nfbnds = nfbnds + 1            

              fbnds(nfbnds) = nd
              fbnds_xy(1:2,nfbnds) = base%xy(1:2,nd)
            ENDIF
          ENDDO
                    
        ENDIF
      ENDDO
      
      
      ALLOCATE(nfbnd2el(nfbnds),fbnd2el(base%mnepn,nfbnds))
      
      nfbnd2el = 0
      
      DO i = 1,nfbnds
        nd = fbnds(i)
        
        DO j = 1,base%nepn(nd)
          el = base%epn(j,nd)
          et = base%el_type(el)
          nv = nverts(et)
          
          DO led = 1,nv
            ged = base%el2ged(el,led)
            edt = base%ed_type(ged)
            
            IF (edt == 10) THEN
              nfbnd2el(i) = nfbnd2el(i) + 1
              fbnd2el(nfbnd2el(i),i) = el
            ENDIF
          ENDDO
          
        ENDDO
      ENDDO      
      
      

      PRINT "(A)", " "
      PRINT "(A,I5)", "Total number of type 0 normal flow boundaries ",num
      PRINT "(A,I5)", "Max number of nodes in a flow boundary ",nmax
      PRINT "(A)", " "
      
      nbou = base%nbou
      
      ALLOCATE(ax(nmax,nbou),cx(nmax,nbou),bx(nmax-1,nbou),dx(nmax-1,nbou))
      ALLOCATE(ay(nmax,nbou),cy(nmax,nbou),by(nmax-1,nbou),dy(nmax-1,nbou)) 
      ALLOCATE(dt(nmax,nbou))
      

      
      RETURN
      END SUBROUTINE spline_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      
      SUBROUTINE calc_cubic_spline(coord,bou,n,btype,sig,a,b,c,d,dt)

      USE globals, ONLY: base
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: coord
      INTEGER, INTENT(IN) :: bou
      INTEGER, INTENT(IN) :: n
      INTEGER, INTENT(IN) :: btype
      INTEGER :: i      
      INTEGER :: info
      
      REAL(rp), INTENT(IN) :: sig
      REAL(rp), INTENT(OUT) , DIMENSION(:) :: a
      REAL(rp), INTENT(OUT), DIMENSION(:) :: b,c,d
      REAL(rp), INTENT(OUT), DIMENSION(:) :: dt
      REAL(rp) :: mult
      REAL(rp) :: x1,y1,x2,y2
      REAL(rp), DIMENSION(n) :: Ml,Md,Mu,v
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: M
      INTEGER, DIMENSION(n) :: ipiv
     
     
      DO i = 1,n-1
      
        dt(i) = 1d0/(real(n,rp)-1d0)
        
!         x1 = base%xy(1,base%fbnds(i,bou))
!         y1 = base%xy(2,base%fbnds(i,bou))
!         
!         x2 = base%xy(1,base%fbnds(i+1,bou))
!         y2 = base%xy(2,base%fbnds(i+1,bou))
!         
!         dt(i) = sqrt((x2-x1)**2 + (y2-y1)**2)
      
      ENDDO
      

      ! Load nodal boundary coordinates 
      DO i = 1,n
        a(i) = base%xy(coord,base%fbnds(i,bou))
      ENDDO      
     



      IF (sig < 1d-12) THEN   ! No tension      
      
!           ! Set up matrix 
          
          Ml(1) = 0d0
          Md(1) = 1d0
          Mu(1) = 0d0
          DO i = 2,n-1     
            Ml(i) = dt(i-1)
            Md(i) = 2d0*(dt(i-1)+dt(i))
            Mu(i) = dt(i)
          ENDDO
          Ml(n) = 0d0
          Md(n) = 1d0
          Mu(n) = 0d0

          ! Set up RHS
          c(1) = 0d0
          DO i = 2,n-1
            c(i) = 3d0/dt(i)*(a(i+1)-a(i)) - 3d0/dt(i-1)*(a(i)-a(i-1)) 
          ENDDO
          c(n) = 0d0
                   
         IF (btype == 1 .OR. btype == 11 .OR. btype == 21) THEN    ! periodic b.c.'s for islands, no longer tridiagonal
         
           ! Set up matrix 
           ALLOCATE(M(n,n))
           M = 0d0
           DO i = 2,n-1
             M(i,i-1) = Ml(i)
             M(i,i) = Md(i)
             M(i,i+1) = Mu(i)
           ENDDO
           
           ! boundary conditions
            M(1,1) = 1d0                                                ! second derivatives equal
            M(1,n) = -1d0
        
            M(n,1) = -2d0*dt(n-1)/3d0 - 2d0*dt(1)/3d0                   ! first derivatives equal
            M(n,2) = -dt(1)/3d0
            M(n,n-1) = -dt(n-1)/3d0
            c(n) = (a(n)-a(n-1))/dt(n-1) - (a(2)-a(1))/dt(1)      
   
           ! Solve system for c coefficients      
           CALL DGESV(n,1,M,n,ipiv,c,n,info) 
           DEALLOCATE(M)
         ELSE
           ! Solve system for c coefficients
           CALL DGTSV(n,1,Ml,Md,Mu,c,n,info)            
         ENDIF        

      ELSE   ! With tension (had to multiply the LHS of Palucci's notes by 2 )

          ! Set up matrix 
          Ml(1) = 0d0
          Md(1) = 1d0
          Mu(1) = 0d0
          DO i = 2,n-1        
            Ml(i) = 2d0*(1d0/dt(i-1) - sig/sinh(sig*dt(i-1)))/sig**2
            Md(i) = 2d0*(sig*cosh(sig*dt(i-1))/sinh(sig*dt(i-1)) - 1d0/dt(i-1) + &
                     sig*cosh(sig*dt(i))/sinh(sig*dt(i)) - 1d0/dt(i))/sig**2
            Mu(i) = 2d0*(1d0/dt(i) - sig/sinh(sig*dt(i)))/sig**2            
          ENDDO
          Ml(n) = 0d0
          Md(n) = 1d0
          Mu(n) = 0d0

          ! Set up RHS
          c(1) = 0d0
          DO i = 2,n-1
            c(i) = (a(i+1)-a(i))/dt(i) - (a(i)-a(i-1))/dt(i-1)
          ENDDO
          c(n) = 0d0
          
        ! Solve system for c coefficients
        CALL DGTSV(n,1,Ml,Md,Mu,c,n,info)             
          
      ENDIF
      
      

      
!       ! Solve system for c coefficients, forward sweep
!       ! RHS is v 
!       DO i = 2,n
!         mult = Ml(i)/Md(i-1)
!         Md(i) = Md(i) - mult*Mu(i-1)
!         v(i) = v(i) - mult*v(i-1)
!       ENDDO
! 
!       ! Solve system for c coefficients, backward sweep
!       c(n) = v(n)/Md(n)
!       DO i = n-1,1,-1
!         c(i) = (v(i) - Mu(i)*c(i+1))/Md(i)
!       ENDDO
 
      ! Solve for other coefficients d and b
      DO i = 1,n-1      
        d(i) = (c(i+1)-c(i))/(3d0*dt(i))
        b(i) = (a(i+1)-a(i))/dt(i) - dt(i)*(2d0*c(i)+c(i+1))/3d0
      ENDDO
      
      
      
      
      RETURN
      END SUBROUTINE calc_cubic_spline
                  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   

      SUBROUTINE eval_cubic_spline(t,ti,a,b,c,d,f,fp,fpp)
      
      IMPLICIT NONE
      
      REAL(rp), INTENT(IN) :: t,ti
      REAL(rp), INTENT(IN) :: a,b,c,d
      REAL(rp), INTENT(OUT) :: f
      REAL(rp), INTENT(OUT), OPTIONAL :: fp,fpp
      
      f = a + b*(t-ti) + c*(t-ti)**2 + d*(t-ti)**3
      
      IF (PRESENT(fp)) THEN
        fp = b + 2d0*c*(t-ti) + 3d0*d*(t-ti)**2
      ENDIF 
      
      IF (PRESENT(fpp)) THEN
        fpp = 2d0*c + 6d0*d*(t-ti)
      ENDIF
      
      
      
      RETURN
      END SUBROUTINE eval_cubic_spline
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE evaluate(r,dt,ti,xr,ax,bx,cx,dx,ay,by,cy,dy,x,y,error_flag)
                 
      IMPLICIT NONE      
      
      REAL(rp), INTENT(INOUT) :: r
      REAL(rp), INTENT(IN) :: dt
      REAL(rp), INTENT(IN) :: ti
      REAL(rp), DIMENSION(:), INTENT(IN) :: xr
      REAL(rp), INTENT(IN) :: ax,bx,cx,dx,ay,by,cy,dy
      REAL(rp), INTENT(OUT) :: x,y
      INTEGER, INTENT(OUT) :: error_flag
      
      REAL(rp) :: r0
      REAL(rp) :: dr
      REAL(rp), PARAMETER :: it_tol = 1d-5      
      INTEGER :: try
      INTEGER :: ntry
      
      r0 = r
      PRINT*, "R0 = ",r0
      CALL newton(r,dt,ti,xr,ax,bx,cx,dx, &
                             ay,by,cy,dy, &
                             x,y,error_flag)
              


      ! Try new initial guess if minimum was not found in (-1,1) interval              
      IF (abs(r)-1d0 > it_tol) THEN
        PRINT "(A,F24.17)", "WARNING: R VALUE NOT FOUND IN INTERVAL, R = ", r
        PRINT "(A)", "  trying negative of r0 value..."     
        r = r0*-1d0
        CALL newton(r,dt,ti,xr,ax,bx,cx,dx, &
                               ay,by,cy,dy, &
                               x,y,error_flag)           
            IF (abs(r)-1d0 < it_tol .and. error_flag == 0) THEN
              
            ELSE 
              PRINT "(A,F24.17)", "WARNING: R VALUE NOT FOUND IN INTERVAL, R = ", r                
            ENDIF                               
      ENDIF
              
      IF (error_flag) THEN
        ntry = 20
        r0 = 0d0
        dr = 1d0/(real(ntry,rp)-1d0)
try_loop: DO try = 1,2*ntry+1
            PRINT "(A,F24.17)", "  trying another r0 value...",r0        
            r = r0
            CALL newton(r,dt,ti,xr,ax,bx,cx,dx, &
                                   ay,by,cy,dy, &
                                   x,y,error_flag)                   
            r0 = r0*-1d0
            IF (mod(try,2) == 1) THEN
              r0 = r0+dr
            ENDIF
            IF (abs(r)-1d0 < it_tol .and. error_flag == 0) THEN
              EXIT try_loop
            ELSE 
              PRINT "(A,F24.17)", "WARNING: R VALUE NOT FOUND IN INTERVAL, R = ", r                
            ENDIF
          ENDDO try_loop
              
        ENDIF

              
!         ! Evaluate spline at specified parameter value (no distance minimiztion)              
!         tpt = .5d0*dt*(r + 1d0) + ti               
!         CALL eval_cubic_spline(tpt,ti,ax,bx,cx,dx,x)
!         CALL eval_cubic_spline(tpt,ti,ay,by,cy,dy,y)              
                                                                    
      
      
      RETURN
      END SUBROUTINE 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      SUBROUTINE newton(r,dt,ti,xr,ax,bx,cx,dx,ay,by,cy,dy,x,y,error_flag)
      
      IMPLICIT NONE      
      
      INTEGER :: it,maxit
      REAL(rp), INTENT(INOUT) :: r
      REAL(rp), INTENT(IN) :: dt
      REAL(rp), INTENT(IN) :: ti
      REAL(rp), INTENT(IN) :: xr(2)
      REAL(rp), INTENT(IN) :: ax,bx,cx,dx,ay,by,cy,dy
      REAL(rp), INTENT(OUT) :: x,y
      INTEGER, INTENT(OUT), OPTIONAL :: error_flag
      REAL(rp) :: t
      REAL(rp) :: f,fp,fpp,g,gp,gpp
      REAL(rp) :: d,dp,tol
      
      tol = 1d-8
      maxit = 1000
      
      t = .5d0*dt*(r + 1d0) + ti          ! initial guess for iteration 
      
iter: DO it = 1,maxit

        CALL eval_cubic_spline(t,ti,ax,bx,cx,dx,f,fp,fpp)               
        CALL eval_cubic_spline(t,ti,ay,by,cy,dy,g,gp,gpp)
        
        d = 2d0*(f-xr(1))*fp + 2d0*(g-xr(2))*gp
        dp = 2d0*fp**2 + 2d0*(f-xr(1))*fpp + 2d0*gp**2 + 2d0*(g-xr(2))*gpp
        
        t = t - d/dp
        
        IF (ABS(d) < tol) THEN
!           PRINT*, "iterations", it
!           PRINT*, d
          EXIT iter
        ENDIF
        
        
      ENDDO iter
      
      CALL eval_cubic_spline(t,ti,ax,bx,cx,dx,x)
      CALL eval_cubic_spline(t,ti,ay,by,cy,dy,y)
      
      r =  2d0/dt*(t-ti)-1d0
      
      IF (it >= maxit) THEN
        PRINT "(A,E28.16)", "MAX ITERATIONS EXCEEDED IN FINDING EVALUATION PARAMETER, ERROR: ", ABS(d)        
      ENDIF  
      
      IF (PRESENT(error_flag)) THEN
        IF (ABS(d) > 1d0 .or. abs(r)-1d0 > 1d-5) THEN
          error_flag = 1
        ELSE
          error_flag = 0
        ENDIF
      ENDIF

      
      RETURN 
      END SUBROUTINE newton      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      

      SUBROUTINE update_elxy_spline(mesh,nverts,bou,bnd,dt,ti,ax,bx,cx,dx,ay,by,cy,dy)  
      
      USE globals, ONLY: grid      
      USE curvilinear_nodes_mod, ONLY: edge_coordinates_curved
      USE transformation, ONLY: element_transformation
      
      IMPLICIT NONE
      
      TYPE(grid), INTENT(INOUT) :: mesh
      INTEGER, DIMENSION(:), INTENT(IN) :: nverts
      INTEGER, INTENT(IN) :: bou
      INTEGER, INTENT(IN) :: bnd
      REAL(rp), INTENT(IN) :: dt
      REAL(rp), INTENT(IN) :: ti
      REAL(rp), INTENT(INOUT) :: ax,bx,cx,dx,ay,by,cy,dy 
            
      INTEGER :: ed,pt,nd
      INTEGER :: ged,el,led,et
      INTEGER :: n1,n2
      INTEGER :: found,error_flag
      REAL(rp) :: n1x,n1y,n2x,n2y
      REAL(rp) :: r,xr(2)
      REAL(rp) :: x,y
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: segxy
      
      ALLOCATE(segxy(2,mesh%ctp))
      
            
      found = 0      
edges:DO ed = 1,mesh%nnfbed
      
        IF (mesh%nfbednn(ed,1) == bou .AND. mesh%nfbednn(ed,3) == bnd) THEN
          ged = mesh%nfbedn(ed)
          el = mesh%ged2el(1,ged)
          led = mesh%ged2led(1,ged)
          
          n1 = mesh%fbnds(bnd,bou)
          n2 = mesh%fbnds(bnd+1,bou)
                   
          n1x = mesh%xy(1,n1)
          n1y = mesh%xy(2,n1)
          
          n2x = mesh%xy(1,n2)
          n2y = mesh%xy(2,n2)                   
          
          found = 1
          EXIT edges
        ENDIF
      
      ENDDO edges
      
      IF (found == 1) THEN
      
        DO nd = 1,mesh%ctp-1
          r = mesh%rpts(nd+1) 
        
          xr(1) = .5d0*(1d0-r)*n1x + .5d0*(1d0+r)*n2x
          xr(2) = .5d0*(1d0-r)*n1y + .5d0*(1d0+r)*n2y        
        
          CALL evaluate(r,dt,ti,xr,ax,bx,cx,dx, &
                                   ay,by,cy,dy, &
                                   x,y,error_flag)    
                                   
          IF (abs(r)-1d0 > 1d-4) THEN
            PRINT "(A,F24.17)", "ERROR: R VALUE NOT FOUND IN INTERVAL, R = ", r
            PRINT "(2(F24.17))", abs(r)-1d0
            STOP
          ENDIF
                                 
          segxy(1,nd) = x
          segxy(2,nd) = y
          
          
        
        ENDDO
       
        CALL edge_coordinates_curved(el,mesh%ctp,led,mesh%nnds,nverts,mesh%el_type_spline,mesh%xy,mesh%ect,segxy,mesh%psiv,mesh%elxy_spline)              
      
      ELSE
       
        PRINT*, "ELEMENT NOT FOUND: cannot update element coordinates"
        STOP
        
       
      ENDIF
      
      RETURN
      END SUBROUTINE update_elxy_spline

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
      
      END MODULE calc_spline