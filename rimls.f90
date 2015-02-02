      PROGRAM rimls

      USE globals
      USE allocation
      USE evaluate
      USE basis
      USE write_results

      IMPLICIT NONE

      INTEGER :: el,nd,pt,i,j,ed
      
      CALL read_input()  
      
      CALL sizes()

      CALL read_grid(base)
      
      IF (refinement) THEN
        CALL read_grid(fine)
      ENDIF         
      
      CALL connect(base) 
      
      IF (refinement) THEN
        CALL connect(fine)       
      ENDIF
      
      
      
      
      CALL vandermonde()  
      
      CALL transformation()      
      
      CALL normals(base)
      
      IF (refinement) THEN
        CALL coordinates(fine)
      ELSE
        CALL coordinates(base)
      ENDIF
      

      
      
      
      
      CALL write_normals()

      IF (refinement) THEN
        CALL write_linear_nodes(fine)
      ELSE 
        CALL write_linear_nodes(base)
      ENDIF
      
       
      
      CALL compute_surface()

      
      
      IF (refinement) THEN
        CALL write_rimls_nodes(fine)
      ELSE
        CALL write_rimls_nodes(base)
      ENDIF
      
     
      END PROGRAM rimls