      PROGRAM rimls

      USE globals
      USE allocation
      USE evaluate
      USE basis
      USE write_results

      IMPLICIT NONE

      INTEGER :: el,nd,pt,i,j,ed
      
      CALL version()
      
      CALL read_input()  
      
      CALL sizes()

      CALL read_grid(base)
      
      IF (refinement) THEN
        CALL read_grid(eval)
      ENDIF         
      
      CALL connect(base) 
      
      IF (refinement) THEN
        CALL connect(eval)       
      ENDIF
      
      
      
      
      CALL vandermonde()  
      
      CALL transformation()      
      
      CALL normals(base)
      
      IF (refinement) THEN
        CALL coordinates(eval)
      ELSE
        CALL coordinates(base)
      ENDIF
      

      
      
      
      
      CALL write_normals()

      IF (refinement) THEN
        CALL write_linear_nodes(eval)
      ELSE 
        CALL write_linear_nodes(base)
      ENDIF
      
       
      
      CALL compute_surface()

      
      
      IF (refinement) THEN
        CALL write_rimls_nodes(eval)
        CALL rewrite_fort14(eval)  
        CALL write_elem_nodes(eval)
      ELSE
        CALL write_rimls_nodes(base)
        CALL rewrite_fort14(base)  
        CALL write_elem_nodes(base)
      ENDIF

      
     
      END PROGRAM rimls