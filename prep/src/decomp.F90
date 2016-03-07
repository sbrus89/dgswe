      SUBROUTINE decomp()
      
      USE globals, ONLY: nn,part,nproc,nresnd,nd_g2l,nresel,el_l2g

      IMPLICIT NONE
      
      INTEGER :: nd,el
      INTEGER :: pe
      INTEGER :: n1,n2,n3
      
      INTEGER :: proc(3)
      
   
      ALLOCATE(nresnd(nproc))
      ALLOCATE(nd_g2l(2,nn))
      
      nresnd = 0
      DO nd = 1,nn
        pe = part(nd)
        nresnd(pe) = nresnd(pe) + 1 ! count the number of resident nodes
        
        nd_g2l(1,nd) = pe          ! global node nd belongs to partition pe
        nd_g2l(2,nd) = nresnd(pe)  ! global node nd's local node number is nresnd(pe)
        nd_l2g(nresnd(pe),pe) = nd ! local node number nresnd(pe) on partition pe has global number nd
      ENDDO
      
      DO pe = 1,nproc
        DO el = 1,ne
          n1 = ect(1,el)
          n2 = ect(2,el)
          n3 = ect(2,el)
          
          proc(1) = nd_g2l(1,n1)
          proc(2) = nd_g2l(1,n2)
          proc(3) = nd_g2l(1,n3)
          
          IF((proc(1)==pe).and.(proc(2)==pe).and.(proc(3)==pe))THEN ! element is resident if all nodes are in partition pe
            nresel(pe) = nresel(pe) + 1
            el_l2g(nresel(pe)),pe) = el
          ELSE 
          
            DO i = 1,3
              IF(proc(i) /= pe) THEN
                match = 0
                DO j = 1,nadjpe(pe)
                  IF(adjpen(j,pe) == proc(i)) THEN
                    match = 1
                  ENDIF
                ENDDO            
                IF(match = 0) THEN
                  nadjpe(pe) = nadjpe(pe) + 1
                  adjpen(nadjpe(pe),pe) = proc(i)
                ENDIF
              ENDIF 
            ENDDO
            
          ENDIF
        ENDDO
      ENDDO
      
      

      RETURN
      END SUBROUTINE decomp