      PROGRAM consolidate

      USE globals, ONLY: rp,ne,nn,nverts,ect,xy,depth,el_type, &
                         nope,neta,obseg,obnds,nvel,nbou,fbseg,fbnds,grid_name, &
                         nepn,mnepn,epn,ned,ged2el,ged2nn,ged2led
      USE grid_file_mod
      USE edge_connectivity_mod      
      USE triangulation, ONLY: reference_element_delaunay       

      IMPLICIT NONE
      
      CHARACTER(200) :: grid_file_in     
      CHARACTER(200) :: grid_file_out        
      
      INTEGER :: i,j,k,l,m
      INTEGER :: ed,nd,bou,el,el1,el2
      INTEGER :: n1,n2,n3
      INTEGER :: nodes_kept
      INTEGER :: nbnds
      INTEGER, ALLOCATABLE, DIMENSION(:) :: nadjnds
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: adjnds   
      INTEGER, ALLOCATABLE, DIMENSION(:) :: keep_node
      INTEGER, ALLOCATABLE, DIMENSION(:) :: nd_fine2coarse
      INTEGER, ALLOCATABLE, DIMENSION(:) :: bou_node
      INTEGER, ALLOCATABLE, DIMENSION(:) :: bou_node_coarse
      INTEGER :: ntri
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: tri            
      INTEGER :: nn_coarse,ne_coarse      
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ect_coarse    
      INTEGER, ALLOCATABLE, DIMENSION(:) :: et_coarse
      INTEGER :: neta_coarse
      INTEGER, ALLOCATABLE, DIMENSION(:) :: obseg_coarse
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: obnds_coarse
      INTEGER :: nvel_coarse
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: fbseg_coarse
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: fbnds_coarse
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: fbseg_orient
      INTEGER :: flag(3)
      INTEGER :: found
      INTEGER, ALLOCATABLE, DIMENSION(:) :: delete_el
      INTEGER :: v1,v2,v11,v21,v12,v22
      REAL(rp) :: xa,ya,xc,yc,x1,y1,x2,y2    
      REAL(rp) :: cross_product
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: xy_coarse   
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: depth_coarse     
      
      
      grid_file_in = "/home/sbrus/data-drive/galveston_SL18/grid_dev/v17_cart/galveston_SL18_cart.grd"
      grid_file_out = "coarse.grd"
         
      nverts(1) = 3
      nverts(2) = 4
      nverts(3) = 3
      nverts(4) = 4         
         

      CALL read_header(0,grid_file_in,grid_name,ne,nn)        
      CALL read_coords(nn,xy,depth)
      CALL read_connectivity(ne,ect,el_type)                 
      CALL read_open_boundaries(nope,neta,obseg,obnds)            
      CALL read_flow_boundaries(nbou,nvel,fbseg,fbnds)
      CALL print_grid_info(grid_file_in,grid_name,ne,nn)          
      
      CALL elements_per_node(ne,nn,nverts,el_type,ect,nepn,mnepn,epn)       
      CALL find_edge_pairs(ne,nverts,el_type,ect,nepn,epn,ned,ged2el,ged2nn,ged2led)           
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Find the nodes adjacent each node
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
      
      ALLOCATE(nadjnds(nn))
      ALLOCATE(adjnds(mnepn,nn))
      nadjnds = 0
      
      DO ed = 1,ned
        n1 = ged2nn(1,ed)  ! find the node numbers on each edge
        n2 = ged2nn(2,ed)
        
        nadjnds(n1) = nadjnds(n1) + 1 ! count the nodes adjacent to node n1
        nadjnds(n2) = nadjnds(n2) + 1 ! count the nodes adjacent to node n2
        
        adjnds(nadjnds(n1),n1) = n2 ! node n2 is adjacent to node n1
        adjnds(nadjnds(n2),n2) = n1 ! node n1 is adjacent to node n2                       
      ENDDO
      
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Eliminate node from input (fine) mesh
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!           
      
      ALLOCATE(keep_node(nn))
      ALLOCATE(bou_node(nn))
      ALLOCATE(fbseg_orient(nbou))
      keep_node = -1
      nodes_kept = 0    
      bou_node = 0       
            
      ! Open ocean boundaries
      
      DO bou = 1,nope
        nbnds = obseg(bou)   
        DO i = 1,nbnds   
          nd = obnds(i,bou)    
          bou_node(nd) = 1
        ENDDO
      ENDDO

      DO bou = 1,nope
        nbnds = obseg(bou)           
        
        nd = obnds(1,bou)                 ! keep first node
        keep_node(nd) = 1
        DO i = 1,nadjnds(nd)
          keep_node(adjnds(i,nd)) = 0
        ENDDO 
        nodes_kept = nodes_kept + 1     
        
        
        DO j = 2,nbnds-1
          nd = obnds(j,bou)
          IF (keep_node(nd) == -1) THEN
            keep_node(nd) = 1
            DO i = 1,nadjnds(nd)
              keep_node(adjnds(i,nd)) = 0
            ENDDO
            nodes_kept = nodes_kept + 1          
          ENDIF
        ENDDO
        
        nd = obnds(nbnds,bou)             ! keep last node
        keep_node(nd) = 1
        DO i = 1,nadjnds(nd)
          keep_node(adjnds(i,nd)) = 0
        ENDDO 
        nodes_kept = nodes_kept + 1       

      ENDDO    
      
      neta_coarse = nodes_kept
      
      ! Flow boundaries 
      
      DO bou = 1,nbou
        nbnds = fbseg(1,bou)   
        DO i = 1,nbnds   
          nd = fbnds(i,bou)    
          bou_node(nd) = 1
        ENDDO
      ENDDO      
      
      DO bou = 1,nbou
        nbnds = fbseg(1,bou)           
        
        nd = fbnds(1,bou)             ! keep first node
        keep_node(nd) = 1
        DO i = 1,nadjnds(nd)
          keep_node(adjnds(i,nd)) = 0
        ENDDO 
        nodes_kept = nodes_kept + 1    
        
        DO j = 2,nbnds-1
          nd = fbnds(j,bou)
          IF (keep_node(nd) == -1) THEN
            keep_node(nd) = 1
            DO i = 1,nadjnds(nd)
              IF (bou_node(adjnds(i,nd)) == 1) THEN
                IF (adjnds(i,nd) == fbnds(j-1,bou) .or. adjnds(i,nd) == fbnds(j+1,bou)) THEN ! needed for very narrow channels, so that nodes on the other                 
                  keep_node(adjnds(i,nd)) = 0                                                ! side aren't eliminated
                ENDIF
              ELSE
                keep_node(adjnds(i,nd)) = 0
              ENDIF
            ENDDO
            nodes_kept = nodes_kept + 1        
          ENDIF
        ENDDO
        
        nd = fbnds(nbnds,bou)         ! keep last node
        keep_node(nd) = 1
        DO i = 1,nadjnds(nd)
          keep_node(adjnds(i,nd)) = 0
        ENDDO 
        nodes_kept = nodes_kept + 1    
        
        n1 = fbnds(1,bou)                ! find boundary orientation
        n2 = fbnds(2,bou)

        DO el = 1,ne
          DO i = 1,3
            v1 = mod(i+0,3) + 1
            v2 = mod(i+1,3) + 1                 
            IF ((ect(v1,el) == n1 .and. ect(v2,el) == n2) .or. &
                (ect(v1,el) == n2 .and. ect(v2,el) == n1)) THEN
                
               xc = 0d0
               yc = 0d0
               DO j = 1,3
                 xc = xc + xy(1,ect(i,el)) 
                 yc = yc + xy(2,ect(i,el))
               ENDDO
               xc = xc/3d0
               yc = yc/3d0
               
               x1 = xy(1,n1)
               y1 = xy(2,n1)
               x2 = xy(1,n2)
               y2 = xy(2,n2)
               
               cross_product = (x2-x1)*(yc-y1) - (y2-y1)*(xc-x1)
               
               IF (cross_product < 0d0) THEN
                 fbseg_orient(bou) = -1d0    ! CCW
               ELSE 
                 fbseg_orient(bou) = 1d0     ! CW
               ENDIF
              
              EXIT  
            ENDIF
          ENDDO
        ENDDO        
      ENDDO
      
      nvel_coarse = nodes_kept - neta_coarse
      
            
      ! Interior nodes
      DO nd = 1,nn
        IF (keep_node(nd) == -1) THEN
          keep_node(nd) = 1
          DO i = 1,nadjnds(nd)
            keep_node(adjnds(i,nd)) = 0            
          ENDDO
          nodes_kept = nodes_kept + 1                  
        ENDIF
      ENDDO           
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Keep track of node infromation for coarse mesh
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!           
      
      ALLOCATE(xy_coarse(2,nodes_kept))
      ALLOCATE(depth_coarse(nodes_kept))
      ALLOCATE(bou_node_coarse(nodes_kept))
      ALLOCATE(nd_fine2coarse(nn))
      
      nn_coarse = 0
      nd_fine2coarse = 0
      DO nd = 1,nn
        IF (keep_node(nd) == 1) THEN
          nn_coarse = nn_coarse + 1
          xy_coarse(1,nn_coarse) = xy(1,nd)
          xy_coarse(2,nn_coarse) = xy(2,nd)
          depth_coarse(nn_coarse) = depth(nd)
          bou_node_coarse(nn_coarse) = bou_node(nd)
          nd_fine2coarse(nd) = nn_coarse
        ENDIF        
      ENDDO      
      
      ALLOCATE(obseg_coarse(nope))
      ALLOCATE(obnds_coarse(neta_coarse,nope))
      
      DO bou = 1,nope
        nbnds = 0
        DO i = 1,obseg(bou)
          nd = obnds(i,bou)
          IF (keep_node(nd) == 1) THEN
            nbnds = nbnds + 1
            obnds_coarse(nbnds,bou) = nd_fine2coarse(nd)
          ENDIF
        ENDDO
        obseg_coarse(bou) = nbnds
      ENDDO
      
      ALLOCATE(fbseg_coarse(2,nbou))
      ALLOCATE(fbnds_coarse(nvel_coarse,nbou))      
      
      DO bou = 1,nbou
        nbnds = 0
        DO i = 1,fbseg(1,bou)
          nd = fbnds(i,bou)
          IF (keep_node(nd) == 1) THEN
            nbnds = nbnds + 1
            fbnds_coarse(nbnds,bou) = nd_fine2coarse(nd)
          ENDIF
        ENDDO
        fbseg_coarse(1,bou) = nbnds
        fbseg_coarse(2,bou) = fbseg(2,bou)
      ENDDO
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Triangulate remaining nodes
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                
      
      xa = 0d0
      ya = 0d0      
            
      
      ALLOCATE(tri(3,3*nn_coarse))
      ALLOCATE(ect_coarse(3,3*nn_coarse))
      CALL reference_element_delaunay(nn_coarse,xy_coarse(1,:),xy_coarse(2,:),ntri,tri,xa,ya)      
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Create coarse element connectity table
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
      
      ALLOCATE(et_coarse(ntri))
      DO el = 1,ntri
        et_coarse(el) = 1
      ENDDO      
      
      CALL elements_per_node(ntri,nn_coarse,nverts,et_coarse,tri,nepn,mnepn,epn)
      
      ! Delete elements that are outside the mesh boundaries
      ALLOCATE(delete_el(ntri))
      delete_el = 0
      ne_coarse = 0
      DO el = 1,ntri
        flag = 0 
        xc = 0d0
        yc = 0d0
        DO i = 1,3
          flag(i) = bou_node_coarse(tri(i,el))
          xc = xc + xy_coarse(1,tri(i,el))
          yc = yc + xy_coarse(2,tri(i,el))
        ENDDO
        xc = xc/3d0
        yc = yc/3d0
        
        IF (flag(1) == 1 .and. flag(2) == 1 .and. flag(3) == 1) THEN      ! elements containing 3 boundary nodes are candidates for deletion            
           found = 0
   search: DO bou = 1,nbou
             DO nd = 1,fbseg_coarse(1,bou)-1
               n1 = fbnds_coarse(nd,bou)
               n2 = fbnds_coarse(nd+1,bou)
               DO i = 1,3
                 v1 = mod(i+0,3) + 1
                 v2 = mod(i+1,3) + 1                 
                 IF ((tri(v1,el) == n1 .and. tri(v2,el) == n2) .or. &
                     (tri(v1,el) == n2 .and. tri(v2,el) == n1)) THEN
                     found = 1
                     EXIT search
                 ENDIF
               ENDDO               
             ENDDO
           ENDDO search
           
           IF (found == 1) THEN
             x1 = xy_coarse(1,n1)
             y1 = xy_coarse(2,n1)
             x2 = xy_coarse(1,n2)
             y2 = xy_coarse(2,n2)
           
             cross_product = fbseg_orient(bou)*((x2-x1)*(yc-y1) - (y2-y1)*(xc-x1))
           
             IF (cross_product <= 0d0) THEN
               delete_el(el) = 1
             ENDIF
           ELSE
             delete_el(el) = 1
           ENDIF
        ENDIF
          
      ENDDO

      

      
      
      ! Fix situation where boundary elements do not connect adjacent boundary nodes (because elements extend across narrow islands)
      DO bou = 1,nbou
        nbnds = fbseg_coarse(1,bou)
        DO i = 1,nbnds-1
          n1 = fbnds_coarse(i,bou)
          n2 = fbnds_coarse(i+1,bou)
          found = 0
 search2: DO j = 1,nepn(n1)
            el = epn(j,n1)

            DO k = 1,3
              v1 = mod(k+0,3) + 1
              v2 = mod(k+1,3) + 1
               IF ((tri(v1,el) == n1 .and. tri(v2,el) == n2) .or. &
                   (tri(v1,el) == n2 .and. tri(v2,el) == n1)) THEN
                   found = 1
                   EXIT search2
               ENDIF
            ENDDO
          ENDDO search2
          
          IF (found == 0) THEN
 search3:   DO j = 1,nepn(n1)
              el1 = epn(j,n1)
              DO k = 1,3
                v11 = mod(k+0,3) + 1
                v21 = mod(k+1,3) + 1  
                DO l = 1,nepn(n2)
                  el2 = epn(l,n2)
                  DO m = 1,3
                    v12 = mod(m+0,3) + 1
                    v22 = mod(m+1,3) + 1
                    IF ((tri(v11,el1) == tri(v12,el2) .and. tri(v21,el1) == tri(v22,el2)) .or. &
                        (tri(v11,el1) == tri(v22,el2) .and. tri(v21,el1) == tri(v12,el2))) THEN
                        
                      delete_el(el2) = 1  
                      
                      IF (bou_node_coarse(tri(v11,el1)) == 0) THEN
                        n3 = tri(v11,el1)
                      ELSE IF (bou_node_coarse(tri(v21,el1)) == 0) THEN
                        n3 = tri(v21,el1) 
                      ENDIF                      
                      
                      tri(1,el1) = n1
                      IF (fbseg_orient(bou) < 0d0) THEN
                        tri(2,el1) = n3
                        tri(3,el1) = n2
                      ELSE
                        tri(2,el1) = n2
                        tri(3,el1) = n3
                      ENDIF
                      
                      EXIT search3
                                            
                    ENDIF
                  ENDDO
                ENDDO
              ENDDO
            ENDDO search3         
          ENDIF
          
        ENDDO
      ENDDO 
      
      
      ! Keep elements containing 3 boundary nodes if the are connected to three interior elements
      DO el1 = 1,ntri
        IF (delete_el(el1) == 1) THEN
        
          flag = 0
          DO i = 1,3
            v11 = mod(i+0,3) + 1
            v21 = mod(i+1,3) + 1
            DO j = 1,3
              nd = tri(j,el1)
              DO k = 1,nepn(nd)
                el2 = epn(k,nd)
                DO m = 1,3
                  v12 = mod(m+0,3) + 1
                  v22 = mod(m+1,3) + 1                  
                  IF ((tri(v11,el1) == tri(v12,el2) .and. tri(v21,el1) == tri(v22,el2)) .or. &
                      (tri(v11,el1) == tri(v22,el2) .and. tri(v21,el1) == tri(v12,el2))) THEN
                    IF (delete_el(el2) == 0) THEN   
                      flag(i) = 1
                    ENDIF
                  ENDIF
                ENDDO
              ENDDO
            ENDDO                                    
          ENDDO
          
          IF (flag(1) == 1 .and. flag(2) == 1 .and. flag(3) == 1) THEN
            delete_el(el1) = 0
          ENDIF      
          
        ENDIF
      ENDDO
      
      
      DO el = 1,ntri
        IF (delete_el(el) == 0) THEN
          ne_coarse = ne_coarse + 1          
          DO i = 1,3
            ect_coarse(i,ne_coarse) = tri(i,el)
          ENDDO
        ENDIF      
      ENDDO
        
      
      

      CALL print_grid_info(grid_file_out,grid_name,ne_coarse,nn_coarse) 
      PRINT "(A,F8.5)", "element reduction factor: ", real(ne,rp)/real(ne_coarse,rp)
      PRINT "(A,F8.5)", "node reduction factor: ", real(nn,rp)/real(nn_coarse,rp)
      PRINT "(A)", ""
      CALL write_header(grid_file_out,grid_name,ne_coarse,nn_coarse)
      CALL write_coords(nn_coarse,xy_coarse,depth_coarse)
      CALL write_connectivity(ne_coarse,ect_coarse,et_coarse,nverts)
      CALL write_open_boundaries(nope,neta_coarse,obseg_coarse,obnds_coarse)
      CALL write_flow_boundaries(nbou,nvel_coarse,fbseg_coarse,fbnds_coarse)

      END PROGRAM consolidate