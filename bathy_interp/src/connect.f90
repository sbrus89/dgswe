      SUBROUTINE connect(mesh)

      USE globals, ONLY: rp,nverts,grid
      USE edge_connectivity_mod
                         

      IMPLICIT NONE
      
      TYPE(grid) :: mesh
         


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
                           
      CALL find_element_edges(mesh%ne,mesh%ned, &
                              mesh%ged2el,mesh%ged2led, &
                              mesh%el2ged)
      
      CALL find_neighbor_elements(mesh%ne,mesh%ned, &
                                  mesh%ged2el,mesh%ged2led, &
                                  mesh%el2el)
        
      mesh%nred = 0
      CALL print_connect_info(mesh%mnepn,mesh%ned,mesh%nied,mesh%nobed,mesh%nfbed,mesh%nnfbed,mesh%nred)
          

      RETURN
      END SUBROUTINE connect