      SUBROUTINE boundary_check2(pt,el_closest,nneigh,neighbors)
      
      USE globals, ONLY: base,nverts
      
      IMPLICIT NONE
      
      INTEGER :: i,j,q,n,el,eln,pt
      INTEGER :: et,nv,ged,etype
      INTEGER :: nneigh,el_closest
      INTEGER :: nbel,qflag,nflag,nq,qpos
      INTEGER :: neighbors(base%ne),connected(base%ne),queue(base%ne),nex,nin

      nbel = 0
      
   
      
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
      
      
      nneigh = nq
      
!       IF (pt == 1555) THEN
!         PRINT*, nq
!         DO i = 1,nq
!           PRINT*,queue(i)      
!         ENDDO      
!       ENDIF

      
      
      
      
      RETURN
      END SUBROUTINE boundary_check2     

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!