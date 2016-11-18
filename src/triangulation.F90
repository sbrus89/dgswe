      MODULE triangulation
      
      USE globals, ONLY: rp
      
      INTEGER :: ncall = 0
      
      CONTAINS
      

      SUBROUTINE reference_element_delaunay(et,npt,ri,si,ntri,ect)

      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: et
      INTEGER, INTENT(IN) :: npt
      REAL(rp), DIMENSION(:), INTENT(IN) :: ri
      REAL(rp), DIMENSION(:), INTENT(IN) :: si
      INTEGER, INTENT(OUT) :: ntri
      INTEGER, DIMENSION(:,:), INTENT(OUT) :: ect

      
      INTEGER :: i,j,m,k
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
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ltri   
      INTEGER :: lwk
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: iwk
      INTEGER, DIMENSION(:), ALLOCATABLE :: tri_flag
      INTEGER :: n1,n2      
      REAL(rp) :: rdiff,sdiff,tol
      CHARACTER(3) :: nout

      n = npt
        
      
      IF (mod(et,2) == 1) THEN
      
        ALLOCATE(r(n+1),s(n+1),x(n+1),y(n+1))
        DO i = n,1,-1    ! add extra node so that first three aren't collinear
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
          
        m = m+1
      ENDDO     
      
      ! Edge 1
      DO j = 2,nx
        i = nx
        
        k = (nx-j)*nx + i
        
        r(k) = ri(m)
        s(k) = si(m)
          
        m = m+1
      ENDDO     
      
      ! Edge 2
      DO i = nx-1,1,-1
        j = nx
        
        k = (nx-j)*nx + i        
        
        r(k) = ri(m)
        s(k) = si(m)
          
        m = m+1      
      ENDDO      
      
      ! Edge 3
      DO j = nx-1,2,-1
        i = 1
        
        k = (nx-j)*nx + i         
        
        r(k) = ri(m)
        s(k) = si(m)
        
        m = m+1
      ENDDO    
      
      ! Interior nodes
      DO i = 2,nx-1
        DO j = 2,nx-1
        
          k = (nx-j)*nx + i              
        
          r(k) = ri(m)
          s(k) = si(m)
          
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
      CALL trlist ( ncc, lcc, n, list, lptr, lend, nrow, nt, ltri, lct, ier )

      
      ntri = 0
      DO i = 1,nt
        remove = 0
 nodes: DO j = 1,3
          IF (ltri(j,i) == 1) THEN
            remove = 1
            EXIT nodes                        
          ENDIF 
        ENDDO nodes
        
        IF (remove == 0) THEN
          ntri = ntri + 1
          DO j = 1,3
            ect(j,ntri) = ltri(j,i)-1    ! -1 to account for extra node
          ENDDO
        ENDIF
      ENDDO



!       ntri = 0
!       DO i = 1,nt     
!           ntri = ntri + 1
!           DO j = 1,3
!             ect(j,ntri) = ltri(j,i)
!           ENDDO
!       ENDDO

      RETURN
      END SUBROUTINE reference_element_delaunay
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
      
      SUBROUTINE delaunay_triangulation(npt,ri,si,ntri,ect,ra,sa)

      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: npt
      REAL(rp), DIMENSION(:), INTENT(IN) :: ri
      REAL(rp), DIMENSION(:), INTENT(IN) :: si
      INTEGER, INTENT(OUT) :: ntri
      INTEGER, DIMENSION(:,:), INTENT(OUT) :: ect
      REAL(rp), INTENT(IN), OPTIONAL :: ra
      REAL(rp), INTENT(IN), OPTIONAL :: sa

      
      INTEGER :: i,j
      INTEGER :: remove
      INTEGER :: n
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: r,s      
      INTEGER, ALLOCATABLE, DIMENSION(:) :: list
      INTEGER, ALLOCATABLE, DIMENSION(:) :: lptr
      INTEGER, ALLOCATABLE, DIMENSION(:) :: lend
      INTEGER :: lnew
      INTEGER, ALLOCATABLE, DIMENSION(:) :: near
      INTEGER, ALLOCATABLE, DIMENSION(:) :: next
      INTEGER :: ier
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: dist
      REAL(rp) :: wx1,wx2,wy1,wy2
      INTEGER :: ncc 
      INTEGER :: lcc(1),lct(1)
      INTEGER :: nrow
      INTEGER :: nt
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ltri      
      CHARACTER(1) :: nout

      n = npt
      
      ALLOCATE(r(n+1),s(n+1))
      DO i = n,1,-1    ! add extra node so that first three aren't collinear
        r(i+1) = ri(i)
        s(i+1) = si(i)
      ENDDO
      IF (PRESENT(ra) .and. PRESENT(sa)) THEN
        r(1) = ra
        s(1) = sa
      ELSE
        r(1) = -1.1d0
        s(1) = -2d0
      ENDIF
      n = n+1      
      
      
      ALLOCATE(list(6*n),lptr(6*n),lend(n))
      ALLOCATE(near(n),next(n),dist(n))
      CALL trmesh(n,r,s,list,lptr,lend,lnew,near,next,dist,ier)
      
      IF (ier < 0) THEN
        PRINT*, "Error in trmesh, ier = ",ier
        STOP
      ENDIF
      
      ncc = 0 
      lcc = 0
      nrow = 6      

      
!       wx1 = r(1)
!       wx2 = wx1
!       wy1 = s(1)
!       wy2 = wy1
!       DO i = 2,N
!         IF (r(i) .LT. wx1) wx1 = r(i)
!         IF (r(i) .GT. wx2) wx2 = r(i)
!         IF (s(i) .LT. wy1) wy1 = s(i)
!         IF (s(i) .GT. wy2) wy2 = s(i)
!       ENDDO 

!       ncall = ncall + 1
!       WRITE(nout,"(I1)") ncall
!       OPEN(UNIT=90,file="tri"//nout//".ps")
!       CALL trplot(90,7.5d0,wx1,wx2,wy1,wy2,ncc,lcc,n,r,s,list,lptr,lend,"(test)",.true.,ier)
!       CLOSE(90)    

      ALLOCATE(ltri(nrow,12*n))
      CALL trlist ( ncc, lcc, n, list, lptr, lend, nrow, nt, ltri, lct, ier )

      
      ntri = 0
      DO i = 1,nt
        remove = 0
 nodes: DO j = 1,3
          IF (ltri(j,i) == 1) THEN
            remove = 1
            EXIT nodes                        
          ENDIF 
        ENDDO nodes
        
        IF (remove == 0) THEN
          ntri = ntri + 1
          DO j = 1,3
            ect(j,ntri) = ltri(j,i)-1    ! -1 to account for extra node
          ENDDO
        ENDIF
      ENDDO

      RETURN
      END SUBROUTINE delaunay_triangulation      
      
      END MODULE triangulation