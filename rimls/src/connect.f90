      SUBROUTINE connect(mesh)

      USE globals, ONLY: rp,nverts,grid
      USE edge_connectivity_mod
                         

      IMPLICIT NONE
      INTEGER :: el,el1,el2,led1,led2,i,seg,m,ged,ed,ed1,ed2,nd,ne,nn,ned
      INTEGER :: n1,n2,nnds
      INTEGER :: n1ed1,n2ed1,n1ed2,n2ed2
      INTEGER :: n1bed,n2bed
      INTEGER :: segtype
      INTEGER :: alloc_status
      INTEGER :: found
      INTEGER :: nvert1,nvert2,nvert
      INTEGER :: el_in,el_ex
      REAL(rp) :: x1,x2,x3,y1,y2,y3
      
      TYPE(grid) :: mesh
      
      INTEGER, ALLOCATABLE, DIMENSION(:) :: nfbnd_temp      
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ged2nn_temp, ged2el_temp, ged2led_temp ! temporary arrays for edge connectivity     
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: edflag
      INTEGER :: mnepn ! maximum number of elements per node

      
      ne = mesh%ne
      nn = mesh%nn

      ALLOCATE(ged2nn_temp(2,3*ne),ged2el_temp(2,3*ne),ged2led_temp(2,3*ne),STAT = alloc_status)
      IF(alloc_status /= 0) THEN
        PRINT*, 'Allocation error: ged2nn_temp,ged2el_temp,ged2led_temp'
      ENDIF 

      CALL elements_per_node(mesh%ne,mesh%nn,nverts, &
                             mesh%el_type,mesh%ect, &
                             mesh%nepn,mesh%mnepn,mesh%epn)

      CALL find_edge_pairs(mesh%ne,nverts,mesh%el_type,mesh%ect, &
                           mesh%nepn,mesh%epn,mesh%ned, &
                           mesh%ged2el,mesh%ged2nn,mesh%ged2led)      

      CALL find_interior_edges(mesh%ned,mesh%ged2el, &
                               mesh%nbed,mesh%bedn,mesh%nied,mesh%iedn, &
                               mesh%ed_type,mesh%recv_edge)      
      
      CALL find_open_edges(mesh%nope,mesh%obseg,mesh%obnds,mesh%ged2nn, &
                           mesh%nbed,mesh%bedn,mesh%nobed,mesh%obedn, &
                           mesh%ed_type,mesh%recv_edge)
                           
      CALL find_flow_edges(mesh%nbou,mesh%fbseg,mesh%fbnds,mesh%ged2nn, &
                           mesh%nbed,mesh%bedn,mesh%nnfbed,mesh%nfbedn,mesh%nfbednn,mesh%nfbed,mesh%fbedn, &
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