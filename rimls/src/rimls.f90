      PROGRAM rimls

      USE globals
      USE allocation
      USE evaluate
      USE basis
      USE write_results
      USE grid_file_mod

      IMPLICIT NONE

      INTEGER :: el,nd,pt,i,j,ed
      INTEGER :: myrank
      
      myrank = 0      
      
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
      
      
      
      
      CALL ref_elem_coords()  
      
      CALL shape_functions_eval()     
      
      IF (refinement) THEN
        CALL read_curve_file(myrank,curve_file,ctp,eval%nbou,eval%xy,eval%bndxy)       
      ELSE
        CALL read_curve_file(myrank,curve_file,ctp,base%nbou,base%xy,base%bndxy)     
      ENDIF
      
      
      
      IF (refinement) THEN
        CALL curvilinear(eval)
      ELSE 
        CALL curvilinear(base)       
      ENDIF      
      
      
      IF (refinement) THEN
        CALL bathymetry_eval_nodes(eval)
      ELSE
        CALL bathymetry_eval_nodes(base)
      ENDIF
      
      
      CALL normals(base)
      
      IF (refinement) THEN
        CALL group_eval_coordinates(eval)
      ELSE
        CALL group_eval_coordinates(base)
      ENDIF
      

      
      
      
      
      CALL write_normals()

      IF (refinement) THEN
        CALL write_linear_nodes(eval)
      ELSE 
        CALL write_linear_nodes(base)
      ENDIF
      
       
            

      IF (nrpt > 0) THEN      
        CALL compute_random()
      ELSE 
        CALL compute_surface()
      ENDIF
      
      
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