      PROGRAM rimls

      USE globals, ONLY: nel_type,nverts,eval,base
      USE allocation, ONLY: sizes
      USE evaluate
      USE write_results
      USE grid_file_mod, ONLY: read_curve_file
      USE find_element, ONLY: find_element_init
      USE transformation, ONLY: init_vandermonde
      USE group_coordinates, ONLY: separate_coordinates, gather_coordinates
      USE version, ONLY: version_information

      IMPLICIT NONE

      INTEGER :: el,nd,pt,i,j,ed
      INTEGER :: myrank
   
      
      CALL version_information(unit=6)
      
      CALL read_input()  
      
      CALL sizes(base)
      CALL sizes(eval)

      CALL read_grid(base)      
      CALL read_grid(eval)                         
      
      CALL connect(base) 
      CALL connect(eval)       

      
      myrank = 0      
      CALL read_curve_file(myrank,base%curve_file,base%ctp,base%nbou,base%xy,base%bndxy)      
      CALL read_curve_file(myrank,eval%curve_file,eval%ctp,eval%nbou,eval%xy,eval%bndxy)       
      CALL curvilinear(base)     
      CALL curvilinear(eval)

      CALL bathymetry_eval_nodes(eval)
      CALL bathymetry_base_nodes(base)
      CALL element_center_nodes(base)      
    
            
      
      CALL find_element_init(nel_type,nverts,base%np,base%nnds,base%nn,base%xy,base%nepn,base%epn)           

      CALL separate_coordinates(eval)            
      CALL init_vandermonde(nel_type,base%np)
      
      CALL write_normals()
      CALL write_linear_nodes(eval)  
      CALL write_curved_coordinates(eval)      
      
      
      CALL compute_surface()

      
      
      CALL gather_coordinates(eval)
      CALL write_rimls_nodes(eval) 
      CALL write_elem_nodes(eval,base)
      CALL rewrite_fort14(eval) 
      CALL write_hb_file(eval)      

      
     
      END PROGRAM rimls