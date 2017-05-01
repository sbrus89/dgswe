      MODULE triangulation
      
      USE globals, ONLY: rp
      
      INTEGER :: ncall = 0
      
      CONTAINS
      

      SUBROUTINE reference_element_delaunay(et,p,npt,ri,si,ntri,ect)

      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: et
      INTEGER, INTENT(IN) :: p
      INTEGER, INTENT(IN) :: npt
      REAL(rp), DIMENSION(:), INTENT(IN) :: ri
      REAL(rp), DIMENSION(:), INTENT(IN) :: si
      INTEGER, INTENT(OUT) :: ntri
      INTEGER, DIMENSION(:,:), INTENT(OUT) :: ect

      
      INTEGER :: i,j,m,k,ed
      INTEGER :: remove
      INTEGER :: n
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: r,s      
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: x,y        
      INTEGER, ALLOCATABLE, DIMENSION(:) :: list
      INTEGER, ALLOCATABLE, DIMENSION(:) :: lptr
      INTEGER, ALLOCATABLE, DIMENSION(:) :: lend
      INTEGER :: lnew
      INTEGER, ALLOCATABLE, DIMENSION(:) :: near
      INTEGER, ALLOCATABLE, DIMENSION(:) :: next    
      INTEGER :: ier
      INTEGER :: nit
      INTEGER :: nx
      REAL(rp) :: nx_rp
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: dist
      REAL(rp) :: wx1,wx2,wy1,wy2
      INTEGER :: ncc 
      INTEGER :: lcc(1),lct(1)
      INTEGER :: nrow
      INTEGER :: nt
      INTEGER :: bcnt
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ltri   
      INTEGER :: lwk
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: iwk
      INTEGER, DIMENSION(:), ALLOCATABLE :: tri_flag
      INTEGER, DIMENSION(:), ALLOCATABLE :: k2m
      INTEGER :: n1,n2      
      REAL(rp) :: rdiff,sdiff,tol
      CHARACTER(3) :: nout

      n = npt
        
      
      IF (mod(et,2) == 1) THEN
      
        ALLOCATE(r(n+1),s(n+1),x(n+1),y(n+1))
        DO i = n,1,-1    ! add dummy node so that first three aren't collinear
          r(i+1) = ri(i)
          s(i+1) = si(i)
        ENDDO
        r(1) = -1.1d0
        s(1) = -2d0
        n = n+1  
                
        DO i = 1,n       ! transform to equilateral triangle for better triangluation results
          y(i) = (6d0*s(i) + 2d0)/(4d0*sqrt(3d0))
          x(i) = (6d0*r(i)+2d0+2d0*sqrt(3d0)*y(i))/6d0
        ENDDO 
      ELSE
        ALLOCATE(r(n),s(n),x(n),y(n))
        ALLOCATE(k2m(n))

        nx_rp = sqrt(real(n,rp))
        nx = int(nx_rp)    
        
        ! Put nodes in the order expected by trmshr, left to right starting at top  
        m = 1         
        ! Edge 4
        DO i = 1,nx
          j = 1
        
          k = (nx-j)*nx + i
        
          r(k) = ri(m)
          s(k) = si(m)
          k2m(k) = m
           
          m = m+1
        ENDDO     
      
        ! Edge 1
        DO j = 2,nx
          i = nx
        
          k = (nx-j)*nx + i
        
          r(k) = ri(m)
          s(k) = si(m)
          k2m(k) = m        
          
          m = m+1
        ENDDO     
      
        ! Edge 2
        DO i = nx-1,1,-1
          j = nx
        
          k = (nx-j)*nx + i        
        
          r(k) = ri(m)
          s(k) = si(m)
          k2m(k) = m        
          
          m = m+1      
        ENDDO      
      
        ! Edge 3
        DO j = nx-1,2,-1
          i = 1
        
          k = (nx-j)*nx + i         
        
          r(k) = ri(m)
          s(k) = si(m)
          k2m(k) = m        
        
          m = m+1
        ENDDO    
      
        ! Interior nodes
        DO i = 2,nx-1
          DO j = 2,nx-1
        
            k = (nx-j)*nx + i              
        
            r(k) = ri(m)
            s(k) = si(m)
            k2m(k) = m          
          
            m = m+1
          ENDDO
        ENDDO       
        
        DO i = 1,n       ! transform to equilateral triangle for better triangluation results
          y(i) = (6d0*s(i) + 2d0)/(4d0*sqrt(3d0))
          x(i) = (6d0*r(i)+2d0+2d0*sqrt(3d0)*y(i))/6d0
        ENDDO
      ENDIF
      
      
      ALLOCATE(list(6*n),lptr(6*n),lend(n))
      ALLOCATE(near(n),next(n),dist(n))   
      
      IF (mod(et,2) == 1) THEN
        CALL trmesh(n,x,y,list,lptr,lend,lnew,near,next,dist,ier)
      ELSE      
        nit = 100
        CALL trmshr(n,nx,x,y,nit,list,lptr,lend,lnew,ier)      
      ENDIF

      ncc = 0 
      lcc = 0
      nrow = 6      
      
!       lwk = n
!       ALLOCATE(iwk(2,n))      
!       call delnod (1,ncc,lcc,n,r,s,list,lptr,lend,lnew,lwk,iwk,ier )      

      
      wx1 = r(1)
      wx2 = wx1
      wy1 = s(1)
      wy2 = wy1
      DO i = 2,N
        IF (r(i) .LT. wx1) wx1 = r(i)
        IF (r(i) .GT. wx2) wx2 = r(i)
        IF (s(i) .LT. wy1) wy1 = s(i)
        IF (s(i) .GT. wy2) wy2 = s(i)
      ENDDO 

!       ncall = ncall + 1
!       WRITE(nout,"(I3.3)") ncall
!       OPEN(UNIT=90,file="tri"//nout//".ps")
!       CALL trplot(90,7.5d0,wx1,wx2,wy1,wy2,ncc,lcc,n,r,s,list,lptr,lend,"(test)",.true.,ier)
!       CLOSE(90)    

      

      ALLOCATE(ltri(nrow,12*n))
      CALL trlist (ncc,lcc,n,list,lptr,lend,nrow,nt,ltri,lct,ier)

      
      IF (mod(et,2) == 0) THEN ! map quad nodes back to original order
        DO i = 1,nt 
          DO j = 1,3            
            ltri(j,i) = k2m(ltri(j,i))
          ENDDO          
        ENDDO        
      ENDIF
      
      
      ntri = 0
      DO i = 1,nt
      
        remove = 0        
        
        IF (mod(et,2) == 1) THEN  
        
          bcnt = 0           ! remove triangle if all its nodes are from the hypotenuse
          DO j = 1,3
            DO k = 1,p+1
              m = p+1+k
              IF (ltri(j,i) == m) THEN
                bcnt = bcnt + 1
              ENDIF
            ENDDO          
          ENDDO
          
          IF (bcnt == 3) THEN
            remove = 1
          ENDIF          
          
   nodes: DO j = 1,3        ! remove triangle if it contains the dummy node
            IF (ltri(j,i) == 1) THEN
              remove = 1
              EXIT nodes                        
            ENDIF 
          ENDDO nodes     
          
        ELSE IF (mod(et,2) == 0) THEN
        
          DO ed = 1,4       ! remove triangle if all its nodes are from a single edge
            bcnt = 0
            DO j = 1,3
              DO k = 1,p+1
                m = mod(ed,4)*p + k
                IF (m == 4*p+1) THEN
                  m = 1
                ENDIF
              
                IF (ltri(j,i) ==  m)THEN
                  bcnt = bcnt + 1
                ENDIF
              ENDDO          
            ENDDO     
            IF (bcnt == 3) THEN
              remove = 1
            ENDIF   
          ENDDO
          
        ENDIF
        

            
        
        IF (remove == 0) THEN
          ntri = ntri + 1
          DO j = 1,3
            IF (mod(et,2) == 1) THEN
              ect(j,ntri) = ltri(j,i)-1    ! -1 to account for extra node
            ELSE IF (mod(et,2) == 0) THEN
              ect(j,ntri) = ltri(j,i)
            ENDIF
          ENDDO
        ENDIF
        
      ENDDO
      


      RETURN
      END SUBROUTINE reference_element_delaunay
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
      
      SUBROUTINE delaunay_triangulation(n,x,y,nbou,fbseg,fbnds,ntri,ect,ncbou,cbseg,cbnds)

      IMPLICIT NONE
      
      INTEGER, INTENT(INOUT) :: n
      REAL(rp), DIMENSION(n), INTENT(INOUT) :: x
      REAL(rp), DIMENSION(n), INTENT(INOUT) :: y
      INTEGER :: nbou
      INTEGER, DIMENSION(:,:), INTENT(IN) :: fbseg
      INTEGER, DIMENSION(:,:), INTENT(INOUT) :: fbnds
      INTEGER, INTENT(OUT) :: ntri
      INTEGER, DIMENSION(:,:), INTENT(OUT) :: ect
      INTEGER, INTENT(IN), OPTIONAL :: ncbou
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: cbseg
      INTEGER, DIMENSION(:,:), INTENT(INOUT), OPTIONAL :: cbnds      


      
      INTEGER :: i,j,k,ed
      INTEGER :: remove,found
      INTEGER :: segtype   
      INTEGER :: constraint
      INTEGER, ALLOCATABLE, DIMENSION(:) :: list
      INTEGER, ALLOCATABLE, DIMENSION(:) :: lptr
      INTEGER, ALLOCATABLE, DIMENSION(:) :: lend
      INTEGER :: lnew
      INTEGER, ALLOCATABLE, DIMENSION(:) :: near
      INTEGER, ALLOCATABLE, DIMENSION(:) :: next
      INTEGER :: ier
      INTEGER :: lwk
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: iwk
      INTEGER, ALLOCATABLE, DIMENSION(:) :: iwork      
      INTEGER :: bou,nd
      INTEGER :: nbnds
      INTEGER :: nd1,nd2
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: dist
      INTEGER :: ncc 
      INTEGER, ALLOCATABLE, DIMENSION(:) :: lcc,lct
      INTEGER :: nrow
      INTEGER :: nt
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ltri  
      INTEGER :: vflag(3),vnd(3)
      INTEGER :: vmin,vmax
      INTEGER :: n1,n2  
      INTEGER :: v1,v2
      INTEGER :: el1,el2
      REAL(rp) :: xc,yc
      REAL(rp) :: dc(3),dmax,dmin



      ncc = 0
      ALLOCATE(lcc(1))
      
      constraint = 0
      IF (PRESENT(ncbou)) THEN
        constraint = 1
      ENDIF

      
      PRINT*, "Triangulating nodes..."
      ALLOCATE(list(6*n),lptr(6*n),lend(n))
      ALLOCATE(near(n),next(n),dist(n))
      CALL trmesh(n,x,y,list,lptr,lend,lnew,near,next,dist,ier)
      
      IF (ier < 0) THEN
        PRINT*, "Error in trmesh, ier = ",ier
        STOP
      ENDIF      
      
      
      
      
      IF (constraint == 1) THEN
        PRINT*, "Enforcing feature constraints..."
        lwk = 6*n      
        ALLOCATE(iwk(2,lwk))
        DO bou = 1,ncbou
          nbnds = cbseg(bou)        
          DO i = 1,nbnds-1
             nd1 = cbnds(i,bou)
             nd2 = cbnds(i+1,bou)
           
             CALL edge(nd1,nd2,x,y,lwk,iwk,list,lptr,lend,ier)
!            PRINT*, ier
          ENDDO
        ENDDO 
      ENDIF
      
      PRINT*, "Enforcing boundary constraints..."
      DO bou = 1,nbou
        nbnds = fbseg(1,bou)        
        DO i = 1,nbnds-1
           nd1 = fbnds(i,bou)
           nd2 = fbnds(i+1,bou)
           
           CALL edge(nd1,nd2,x,y,lwk,iwk,list,lptr,lend,ier)
!            PRINT*, ier
        ENDDO
      ENDDO       
      
      
    
               

      PRINT*, "Converting triangulation data ..."
      nrow = 6          
      ALLOCATE(ltri(nrow,12*n))
      ALLOCATE(lct(ncc))
      CALL trlist ( ncc, lcc, n, list, lptr, lend, nrow, nt, ltri, lct, ier )
      PRINT*, ier
      

      
      



      ntri = 0
      DO i = 1,nt 
        remove = 0
        
  bchk: DO bou = 1,nbou
          nbnds = fbseg(1,bou)
          vflag = 0
          vnd = 0
          DO j = 1,3
            nd = ltri(j,i)
      nchk: DO k = 1,nbnds
              IF (nd == fbnds(k,bou)) THEN
                vflag(j) = 1
                vnd(j) = k
                EXIT nchk
              ENDIF
            ENDDO nchk
          ENDDO
                 
          IF (vflag(1) + vflag(2) + vflag(3) == 3) THEN
            vmax = fbnds(maxval(vnd),bou)
            vmin = fbnds(minval(vnd),bou)
            
            IF ((vmin == ltri(1,i) .and. vmax == ltri(3,i)) .or. &
                (vmin == ltri(2,i) .and. vmax == ltri(1,i)) .or. &
                (vmin == ltri(3,i) .and. vmax == ltri(2,i))) THEN
              remove = 1
!               PRINT*, remove
              EXIT bchk
            ENDIF
          ENDIF
        ENDDO bchk     
        
        IF(remove == 0) THEN
          ntri = ntri + 1
          DO j = 1,3              
            ect(j,ntri) = ltri(j,i)             
          ENDDO  
        ENDIF
      ENDDO
      
      
!       remove = 0
!       DO el1 = 1,ntri
!         
!          xc = 0d0
!          yc = 0d0
!          DO j = 1,3
!            k = ect(j,el1)
!            xc = xc + x(k)
!            yc = yc + y(k)
!          ENDDO
!          xc = xc/3d0
!          yc = yc/3d0
!                   
!          DO j = 1,3
!            k = ect(j,el1)         
!            dc(j) = sqrt((x(k)-xc)**2+(y(k)-yc)**2)
!          ENDDO
!          
!          dmin = 1d9
!          dmax = 0d0
!          DO j = 1,3
!            IF (dc(j) < dmin) THEN
!              dmin = dc(j)
!              vmin = j
!            ENDIF
!            IF (dc(j) > dmax) THEN
!              dmax = dc(j)
!            ENDIF
!          ENDDO
!       
!         IF (dmin < dmax*.1d0) THEN
!         
!           n1 = ect(mod(vmin+0,3)+1,el1)
!           n2 = ect(mod(vmin+1,3)+1,el1)
!           nd1 = ect(vmin,el1)
!           
!           found = 0
!    tsrch: DO el2 = 1,ntri
!             DO ed = 1,3
!               v1 = ect(mod(ed+0,3)+1,el2)
!               v2 = ect(mod(ed+1,3)+1,el2)
!               
!               IF (((n1 == v1) .and. (n2 == v2)) .or. &
!                   ((n1 == v2) .and. (n2 == v1))) THEN
!                   nd2 = ect(ed,el2)
!                   found = 1
!                   EXIT tsrch                 
!               ENDIF
!             ENDDO
!           ENDDO tsrch
!           
!           IF (found == 1) THEN
!             ect(mod(vmin+0,3)+1,el1) = nd2            
!             ect(mod(ed+0,3)+1,el2) = nd1  
!             PRINT*, "Swapped element", nd1,nd2
!           ELSE
!             PRINT*, "Error finding neighbor element"
!           ENDIF
!           
!         ENDIF
!         
!       ENDDO
      

      RETURN
      END SUBROUTINE delaunay_triangulation      
      
      END MODULE triangulation