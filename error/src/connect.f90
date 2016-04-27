      SUBROUTINE connect(mesh)

      USE globals, ONLY: rp,nverts,solution,exclude_bndel
      USE edge_connectivity_mod
                         

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