      SUBROUTINE eval_pts()

      USE globals, ONLY: rp,nel_type,eval,mnept,nept,ept
      USE basis, ONLY: tri_nodes,quad_nodes

      IMPLICIT NONE 
      INTEGER :: i,pt,et
      
      mnept = (eval%hbp+1)**2
      ALLOCATE(ept(mnept,2,nel_type))      
      
      DO i = 1,nel_type     

        IF (mod(i,2) == 1) THEN
          CALL tri_nodes(1,eval%np(5),nept(i),ept(:,1,i),ept(:,2,i))
        ELSE
          CALL quad_nodes(1,eval%np(6),nept(i),ept(:,1,i),ept(:,2,i))
        ENDIF
                
      ENDDO
      
      mnept = maxval(nept)
               
      RETURN
      END SUBROUTINE eval_pts
      
     