      SUBROUTINE boundary_check2(pt,el_closest,nneigh,neighbors)
      
      USE globals, ONLY: rp,base,hmin,closest,tree_xy,kdresults,nverts
      
      IMPLICIT NONE
      
      INTEGER :: i,j,q,n,el,eln,pt
      INTEGER :: et,nv,ged,etype
      INTEGER :: nneigh,el_closest
      INTEGER :: nbel,qflag,nflag,nq,qpos
      INTEGER :: neighbors(base%ne),connected(base%ne),queue(base%ne),nex,nin

      nbel = 0
      
      DO i = 1,nneigh
        neighbors(i) = kdresults(i)%idx
      ENDDO      
      
      DO i = 1,nneigh                   ! check if any neighboring elements are on boundaries
        el = neighbors(i)
        et = base%el_type(el)
        nv = nverts(et)
 edges: DO j = 1,nv
          ged = base%el2ged(el,j)
          etype = base%ed_type(ged)
          IF (etype /= 0) THEN
            nbel = nbel + 1
            EXIT edges
          ENDIF
        ENDDO edges     
      ENDDO
      
      IF (nbel == 0) THEN               ! if there are no boundary elements, then use all neighbors     
        RETURN                          ! and return               
      ENDIF
 
      el = el_closest      
      nq = 1
      qpos = 1
      queue(nq) = el

search:DO    
        
        et = base%el_type(el)
        nv = nverts(et)
        
 elems: DO j = 1,nv                     
 
          eln = base%el2el(el,j)        ! find elements connected to el
          IF (eln == 0) THEN
            CYCLE elems
          ENDIF
          
          qflag = 0                     ! check if the element has been added to the queue        
   quen: DO q = 1,nq
            IF (queue(q) == eln) THEN
              qflag = 1
              EXIT quen
            ENDIF
          ENDDO quen
          
          IF (qflag == 0) THEN          ! check if the element is in the k-d radius
            nflag = 0
     neigh: DO n = 1,nneigh
              IF (neighbors(n) == eln) THEN
                nflag = 1
                EXIT neigh
              ENDIF
            ENDDO neigh
          ENDIF
          
          IF (qflag == 0 .and. nflag == 1) THEN  ! add element to the queue if it hasn't been already
            nq = nq + 1                          ! and it's in the kd-radius
            queue(nq) = eln
          ENDIF
          
        ENDDO elems
        
        IF (qpos == nq) THEN            ! exit when all elements have been added
          EXIT search
        ENDIF
        
        qpos = qpos + 1                 ! get ready to check the next element
        el = queue(qpos)
        
      ENDDO search
      
      
      
      
      
      DO i = 1,nq
        neighbors(i) = queue(i)       
      ENDDO
      
      IF (nq == 0) THEN
      
        DO i = 1,nneigh
          PRINT*, kdresults(i)%idx
        ENDDO
        
      ENDIF
      
      nneigh = nq
      
!       IF (pt == 8595) THEN
!         DO i = 1,nq
!           PRINT*,queue(i)      
!         ENDDO      
!       ENDIF

      
      
      
      
      RETURN
      END SUBROUTINE boundary_check2     
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE boundary_check(pt,x,nneigh,neighbors)
      
      USE globals, ONLY: rp,base,hmin,closest,tree_xy,kdresults,nverts
      USE find_element, ONLY: in_element
      
      IMPLICIT NONE
      
      INTEGER :: i,j,el,ed,tstp,elt,bed,ged,nd1,nd2,pt
      INTEGER :: et,nv,etype
      INTEGER :: nneigh,ntstp
      INTEGER :: nbel,flag,n,found
      INTEGER :: neighbors(base%ne),cross(base%ne),ncross,exclude(base%ne),nex,nin
      REAL(rp) :: tol
      REAL(rp) :: dt,t,u,xt(2)
      REAL(rp) :: det,b(2)
      REAL(rp) :: x(3),x1(2),y1(2),x2(2),y2(2)
      
      tol = 1d-12
      
      nbel = 0
      
      dt = .1d0
      ntstp = 2d0/dt
      
      DO i = 1,nneigh                   ! check if any neighboring elements are on boundaries
        el = neighbors(i)
        et = base%el_type(el)
        nv = nverts(et)
 edges: DO j = 1,nv
          ged = base%el2ged(el,j)
          etype = base%ed_type(ged)
          IF (etype /= 0) THEN
            nbel = nbel + 1
            EXIT edges
          ENDIF
        ENDDO edges     
      ENDDO
      
      IF (nbel == 0) THEN               ! if there are no boundary elements, then use all neighbors
                                        ! and return
        DO i = 1,nneigh
          neighbors(i) = kdresults(i)%idx
        ENDDO
        
        RETURN          
      ENDIF
      
      x1(1) = x(1)                      ! if there are boundary elements,
      y1(1) = x(2)                      ! set the xy coordinates of the starting point for the point-neighbor line
      
      nex = 0
      nin = 0
      cross = 0 
      
      DO i = 1,nneigh                   ! loop through all neighbors
        el = kdresults(i)%idx            
        x1(2) = base%xyhc(1,el)         ! set the xy coordinates for the end point of the point-neighbor line
        y1(2) = base%xyhc(2,el)
        
!         PRINT*, " pt ", pt, " neighbor ", el 
        
        ncross = 0
        DO tstp = 1,ntstp-1               ! increment across point-neighbor line
          t = -1d0 + real(tstp*dt,rp)
          xt(1) = .5d0*((1d0-t)*x1(1) + (1d0+t)*x1(2))
          xt(2) = .5d0*((1d0-t)*y1(1) + (1d0+t)*y1(2))     
          
          CALL in_element(i,xt,elt,found)            ! find what element xt is in       
               
          IF (found == 1) THEN     
!           IF (base%bel2bed(elt,1) > 0) THEN  ! determine if element has a boundary edge
!     bedges: DO bed = 1,base%bel2bed(elt,1)   ! determine if point-neighbor line intersects the boundary edge
!               ged = base%bel2bed(elt,bed+1)
!               nd1 = base%ged2nn(1,ged)
!               nd2 = base%ged2nn(2,ged)
!               x2(1) = base%xy(1,nd1)
!               x2(2) = base%xy(1,nd2)
!               y2(1) = base%xy(2,nd1)
!               y2(2) = base%xy(2,nd2)
!               
!               det = (x1(2)-x1(1))*(y2(1)-y2(2)) - (y1(2)-y1(1))*(x2(1)-x2(2))
!               IF (abs(det) < tol) THEN
!                 CYCLE bedges                 ! skip if lines are parallel
!               ENDIF
!               
!               b(1) = x2(1) + x2(2) - x1(1) - x1(2)
!               b(2) = y2(1) + y2(2) - y1(1) - y1(2)
!              
!               t = ((y2(1) - y2(2))*b(1) + (x2(2) - x2(1))*b(2))/det
!               u = ((y1(1) - y1(2))*b(1) + (x1(2) - x1(1))*b(2))/det
!               
!               IF (abs(t) < 1d0 .and. abs(u) < 1d0) THEN  ! determine if lines intersect 
!               
!                 flag = 0 
!                 DO j = 1,nneigh               ! determine if element has already been added to cross list
!                   IF (cross(j) == elt) THEN 
!                     flag = 1
!                   ENDIF                                    
!                 ENDDO
!                 
!                 IF (flag == 0) THEN           ! if element is new, count crossing
!                   ncross = ncross + 1
!                   cross(ncross) = elt
!                 ENDIF
!               ENDIF
!               
!             ENDDO bedges
!           ENDIF
          
          ELSE 
            ncross = ncross + 2
          
          ENDIF
        ENDDO
        
        IF (ncross >= 2) THEN                 ! exclude neighbor element if more than one crossing was found
!           PRINT*, "ncross = ",ncross
!           DO j = 1,ncross
!             PRINT*, cross(j)
!           ENDDO
          nex = nex + 1
          exclude(nex) = el          
        ELSE                                  ! include neighbor element if one or fewer crossings were found
          nin = nin + 1
          neighbors(nin) = el
        ENDIF
      
      ENDDO      
      
     
!       IF (nin == 0) THEN
!         PRINT*, "  All neighbors excluded", x(1),x(2)
!         DO i = 1,nex
!           PRINT*,  exclude(i),base%bel2bed(exclude(i),1)
!         ENDDO
!         
!         PRINT*, "  "
!         PRINT*, "Original neighbors:"
!         DO i = 1,nneigh                  
!           PRINT*, kdresults(i)%idx         
!         ENDDO
!       ENDIF
      
!       IF (nex /= 0 .or. pt == 8595) THEN
!       IF (pt == 8595) THEN
! 
!         PRINT*, "Point: ", pt
!         PRINT*, "  Some neighbors excluded", x(1),x(2)
!         DO i = 1,nex
!           PRINT*,  exclude(i),base%bel2bed(exclude(i),1)
!         ENDDO
!         PRINT*, "Incuded: "
!         DO i = 1,nin
!           PRINT*,  neighbors(i)
!         ENDDO        
!       ENDIF      
!       
!       IF (pt == 1555) STOP
     
      nneigh = nin
      
      
      
      RETURN
      END SUBROUTINE boundary_check


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!