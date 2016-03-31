      SUBROUTINE connect(mesh)

      USE globals, ONLY: rp,nverts,grid
      USE edge_connectivity_mod           

      IMPLICIT NONE
      
      TYPE(grid) :: mesh
      INTEGER :: ed
      INTEGER :: n1,n2
      INTEGER :: el_in,el_ex
      REAL(rp) :: x1,x2,y1,y2
      

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


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! find minimum edge length
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      ALLOCATE(mesh%edlen(mesh%ned),mesh%minedlen(mesh%ne))
      
      mesh%minedlen = 1d6
      DO ed = 1,mesh%ned
        n1 = mesh%ged2nn(1,ed)
        n2 = mesh%ged2nn(2,ed)                
        
        x1 = mesh%xy(1,n1)
        x2 = mesh%xy(1,n2)
        
        y1 = mesh%xy(2,n1)
        y2 = mesh%xy(2,n2)
        
        mesh%edlen(ed) = sqrt((x2-x1)**2 + (y2-y1)**2)        
        
        el_in = mesh%ged2el(1,ed)
        el_ex = mesh%ged2el(2,ed)
        
        IF (mesh%edlen(ed) < mesh%minedlen(el_in)) THEN
          mesh%minedlen(el_in) = mesh%edlen(ed)
        ENDIF
        
        IF (el_ex > 0 ) THEN
          IF (mesh%edlen(ed) < mesh%minedlen(el_ex)) THEN
            mesh%minedlen(el_ex) = mesh%edlen(ed)
          ENDIF          
        ENDIF
        
        
      ENDDO
      
      
      PRINT "(A)", "---------------------------------------------"
      PRINT*, ""



      RETURN
      END SUBROUTINE connect