      SUBROUTINE connect(mesh)

      USE globals, ONLY: rp,nverts,nel_type,solution,exclude_bndel
      USE edge_connectivity_mod, ONLY: elements_per_node,find_edge_pairs,find_interior_edges, &
                                       find_open_edges,find_flow_edges,print_connect_info
      USE curvilinear_nodes_mod, ONLY: shape_functions_linear_at_ctp,eval_coordinates_curved   
                         

      IMPLICIT NONE
      
      INTEGER :: ged
      INTEGER :: et
      INTEGER :: el
      TYPE(solution) :: mesh
         


      CALL elements_per_node(mesh%ne,mesh%nn,nverts, &
                             mesh%el_type,mesh%ect, &
                             mesh%nepn,mesh%mnepn,mesh%epn)

      CALL find_edge_pairs(mesh%ne,nverts,mesh%el_type,mesh%ect, &
                           mesh%nepn,mesh%epn,mesh%ned, &
                           mesh%ged2el,mesh%ged2nn,mesh%ged2led)      

      CALL find_interior_edges(mesh%ned,mesh%ged2el, &
                               mesh%nied,mesh%iedn, &
                               mesh%ed_type,mesh%recv_edge)      
      
      CALL find_open_edges(mesh%nope,mesh%obseg,mesh%obnds,mesh%ged2nn, &
                           mesh%nobed,mesh%obedn, &
                           mesh%ed_type,mesh%recv_edge)
                           
      CALL find_flow_edges(mesh%nbou,mesh%fbseg,mesh%fbnds,mesh%ged2nn, &
                           mesh%nnfbed,mesh%nfbedn,mesh%nfbednn,mesh%nfbed,mesh%fbedn, &
                           mesh%recv_edge,mesh%ed_type)                           
        
      mesh%nred = 0
      CALL print_connect_info(mesh%mnepn,mesh%ned,mesh%nied,mesh%nobed,mesh%nfbed,mesh%nnfbed,mesh%nred)

      
      
      CALL shape_functions_linear_at_ctp(nel_type,mesh%np,mesh%psiv)                   
      CALL eval_coordinates_curved(mesh%ctp,mesh%nnds,nverts,mesh%el_type,mesh%xy,mesh%ect,mesh%fbseg,mesh%fbnds, &
                                   mesh%nnfbed,mesh%nfbedn,mesh%nfbednn,mesh%ged2el,mesh%ged2led, &
                                   mesh%psiv,mesh%bndxy,mesh%elxy)       
      
      
          
      ALLOCATE(mesh%bndel(mesh%ne))
      mesh%bndel = 0          
      
      IF (exclude_bndel) THEN
        DO ged = 1,mesh%ned
          et = mesh%ed_type(ged)
          el = mesh%ged2el(1,ged)
        
          IF (et == 1 .or. et == 10 .or. et == 12) THEN
            mesh%bndel(el) = 1
          ENDIF
        
        ENDDO
      ENDIF
      
      
          
          
      RETURN
      END SUBROUTINE connect