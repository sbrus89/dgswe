      PROGRAM consolidate

      USE globals, ONLY: rp,ne,nn,nverts,ect,xy,depth,el_type, &
                         nope,neta,obseg,obnds,nvel,nbou,fbseg,fbnds,grid_name, &
                         nepn,mnepn,epn,ned,ged2el,ged2nn,ged2led
      USE grid_file_mod
      USE edge_connectivity_mod      
      USE triangulation, ONLY: delaunay_triangulation      
      USE kdtree2_module      

      IMPLICIT NONE
      
      CHARACTER(200) :: grid_file_in     
      CHARACTER(200) :: grid_file_out        
      
      INTEGER :: i,j,k,l,m
      INTEGER :: ed,nd,bou,el,el1,el2
      INTEGER :: n1,n2,n3
      INTEGER :: nodes_kept
      INTEGER :: nbnds
      INTEGER :: bou_type
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
      INTEGER :: nbou_coarse
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: fbseg_coarse
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: fbnds_coarse
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: fbseg_orient
      INTEGER :: ncbou
      INTEGER, ALLOCATABLE, DIMENSION(:) :: cbseg_coarse
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: cbnds_coarse
      INTEGER :: flag(3)
      INTEGER :: found
      INTEGER, ALLOCATABLE, DIMENSION(:) :: delete_el
      INTEGER :: v1,v2,v11,v21,v12,v22
      REAL(rp) :: xa,ya,xc,yc,x1,y1,x2,y2    
      REAL(rp) :: cross_product
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: xy_coarse   
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: depth_coarse  
      
      TYPE(kdtree2), POINTER :: tree_xy      
      TYPE(kdtree2_result), ALLOCATABLE, DIMENSION(:) :: closest     
      INTEGER :: srchdp
      REAL(rp) :: len_avg
      REAL(rp), DIMENSION(:), ALLOCATABLE :: hnd
      INTEGER :: it
      REAL(rp) :: px,py
      REAL(rp) :: lxy
      REAL(rp) :: hpt,l0
      REAL(rp) :: dt
      REAL(rp), DIMENSION(:,:), ALLOCATABLE ::  frc    
      
      
      grid_file_in = "/home/sbrus/data-drive/galveston_SL18/grid_dev/v25_cart/galveston_SL18_cart.grd"
      grid_file_out = "/home/sbrus/data-drive/galveston_SL18/grid_dev/v25_cart/coarse.grd"
!       grid_file_in = "/home/sbrus/data-drive/galveston_SL18/grid_dev/v23_cart/coarse/galveston_SL18_cart_coarse.grd"
!       grid_file_out = "/home/sbrus/data-drive/galveston_SL18/grid_dev/v23_cart/coarse/coarse_x2.grd"
         
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
      
      CALL find_adjacent_nodes(nn,mnepn,ned,ged2nn,nadjnds,adjnds)
      
!       ALLOCATE(nadjnds(nn))
!       ALLOCATE(adjnds(mnepn,nn))
!       nadjnds = 0
!       
!       DO ed = 1,ned
!         n1 = ged2nn(1,ed)  ! find the node numbers on each edge
!         n2 = ged2nn(2,ed)
!         
!         nadjnds(n1) = nadjnds(n1) + 1 ! count the nodes adjacent to node n1
!         nadjnds(n2) = nadjnds(n2) + 1 ! count the nodes adjacent to node n2
!         
!         adjnds(nadjnds(n1),n1) = n2 ! node n2 is adjacent to node n1
!         adjnds(nadjnds(n2),n2) = n1 ! node n1 is adjacent to node n2                       
!       ENDDO
      
      
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
      
      DO bou = 1,nope                     ! flag boundary nodes 
        nbnds = obseg(bou)   
        DO i = 1,nbnds   
          nd = obnds(i,bou)    
          bou_node(nd) = 1
        ENDDO
      ENDDO

      DO bou = 1,nope                     ! determine nodes to keep
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
      
      DO bou = 1,nbou                     ! flag boundary nodes
        bou_type = fbseg(2,bou)
        IF (bou_type /= 30) THEN
          nbnds = fbseg(1,bou)   
          DO i = 1,nbnds   
            nd = fbnds(i,bou)    
            bou_node(nd) = 1
          ENDDO
        ENDIF
      ENDDO      
      
      nbou_coarse = 0
 fbnd:DO bou = 1,nbou                     ! determine boundary nodes to keep
        nbnds = fbseg(1,bou) 
        bou_type = fbseg(2,bou)
        
        IF (bou_type == 30) THEN          ! skip feature constraint nodestrings
          CYCLE fbnd
        ENDIF
        
        nbou_coarse = nbou_coarse + 1
        
        IF (nbnds <= 8) THEN              ! keep small island boundaries
          DO j = 1,nbnds
            nd = fbnds(j,bou)
            IF (keep_node(nd) == -1) THEN
              keep_node(nd) = 1
              nodes_kept = nodes_kept + 1
              PRINT*, "keep node", nd
            ENDIF
          ENDDO
        ELSE
        
          nd = fbnds(1,bou)               ! keep first node
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
        
          nd = fbnds(nbnds,bou)           ! keep last node
          keep_node(nd) = 1
          DO i = 1,nadjnds(nd)
            keep_node(adjnds(i,nd)) = 0
          ENDDO 
          nodes_kept = nodes_kept + 1   
          
        ENDIF
        
        n1 = fbnds(1,bou)                 ! find boundary orientation
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
                 fbseg_orient(nbou_coarse) = -1d0    ! CCW
               ELSE 
                 fbseg_orient(nbou_coarse) = 1d0     ! CW
               ENDIF
              
              EXIT  
            ENDIF
          ENDDO
        ENDDO        
      ENDDO fbnd
      
      nvel_coarse = nodes_kept - neta_coarse
      

      ! Feature constraint nodestrings
      
      ncbou = 0
 cbnd:DO bou = 1,nbou
        nbnds = fbseg(1,bou)
        bou_type = fbseg(2,bou)
        
        IF (bou_type /= 30) THEN
          CYCLE
        ENDIF
        
        ncbou = ncbou + 1
        
        DO j = 1,nbnds
          nd = fbnds(j,bou)
          
          IF (keep_node(nd) == -1) THEN          
            keep_node(nd) = 1
            DO i = 1,nadjnds(nd)
              IF (keep_node(adjnds(i,nd)) /= 1) THEN
                keep_node(adjnds(i,nd)) = 0          
              ENDIF
            ENDDO
            nodes_kept = nodes_kept + 1
          ENDIF          
        ENDDO
        
      ENDDO cbnd
      
      
            
      ! Interior nodes
      
      DO nd = 1,nn
        IF (keep_node(nd) == -1) THEN
          keep_node(nd) = 1
          DO i = 1,nadjnds(nd)
            IF (keep_node(adjnds(i,nd)) /= 1) THEN
              keep_node(adjnds(i,nd)) = 0          
            ENDIF
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
      bou_node_coarse = 0
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
      
      ALLOCATE(fbseg_coarse(2,nbou_coarse))
      ALLOCATE(fbnds_coarse(nvel_coarse,nbou_coarse)) 
      ALLOCATE(cbseg_coarse(ncbou))
      ALLOCATE(cbnds_coarse(nvel_coarse,ncbou))
      
      nbou_coarse = 0
      ncbou = 0 
      DO bou = 1,nbou
        bou_type = fbseg(2,bou)
        nbnds = 0        
        IF (bou_type /= 30) THEN
          nbou_coarse = nbou_coarse + 1
          DO i = 1,fbseg(1,bou)
            nd = fbnds(i,bou)
            IF (keep_node(nd) == 1) THEN
              nbnds = nbnds + 1
              fbnds_coarse(nbnds,nbou_coarse) = nd_fine2coarse(nd)
            ENDIF
          ENDDO
          fbseg_coarse(1,nbou_coarse) = nbnds
          fbseg_coarse(2,nbou_coarse) = fbseg(2,bou)
        ELSE
          ncbou = ncbou + 1
          DO i = 1,fbseg(1,bou)
            nd = fbnds(i,bou)
            IF (keep_node(nd) == 1) THEN
              nbnds = nbnds + 1
              cbnds_coarse(nbnds,ncbou) = nd_fine2coarse(nd)
            ENDIF
          ENDDO
          cbseg_coarse(ncbou) = nbnds
        ENDIF
      ENDDO
      
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Triangulate remaining nodes
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                
      
      xa = 0d0
      ya = 0d0      
            
      
      ALLOCATE(tri(3,3*nn_coarse))
      ALLOCATE(ect_coarse(3,3*nn_coarse))
      CALL delaunay_triangulation(nn_coarse,xy_coarse(1,:),xy_coarse(2,:),ncbou,cbseg_coarse,cbnds_coarse,ntri,tri,xa,ya) 

      ALLOCATE(et_coarse(ntri))
      DO el = 1,ntri
        et_coarse(el) = 1
      ENDDO      
            
    
      
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Create coarse element connectity table
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      PRINT*, "Creating coarse grid data..." 
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
   search: DO bou = 1,nbou_coarse
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
      DO bou = 1,nbou_coarse
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
        
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Adjust node locations to improve mesh quality
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
      CALL elements_per_node(ne_coarse,nn_coarse,nverts,et_coarse,ect_coarse,nepn,mnepn,epn) 
      CALL find_edge_pairs(ne_coarse,nverts,et_coarse,ect_coarse,nepn,epn,ned,ged2el,ged2nn,ged2led) 
      CALL find_adjacent_nodes(nn_coarse,mnepn,ned,ged2nn,nadjnds,adjnds)
      
      ALLOCATE(hnd(nn_coarse))
      DO n1 = 1,nn_coarse
        len_avg = 0d0
        DO j = 1,nadjnds(n1)
          n2 = adjnds(j,n1) 
          x1 = xy_coarse(1,n1)
          x2 = xy_coarse(1,n2)
          y1 = xy_coarse(2,n1)
          y2 = xy_coarse(2,n2)
          
          len_avg = len_avg + sqrt((x1-x2)**2+(y1-y2)**2)
        ENDDO
        hnd(n1) = len_avg/real(nadjnds(n1),rp)
      ENDDO
      
!       srchdp = 20
!       IF (nn_coarse < srchdp) THEN
!         srchdp = nn_coarse
!       ENDIF      
!       
!       tree_xy => kdtree2_create(xy_coarse(1:2,1:nn_coarse), rearrange=.true., sort=.true.)
!       ALLOCATE(closest(MAX(nn_coarse,srchdp)))
      
      
      bou_node = 0
      
      DO bou = 1,nope                    
        nbnds = obseg_coarse(bou)   
        DO i = 1,nbnds   
          nd = obnds_coarse(i,bou)    
          bou_node(nd) = 1
        ENDDO
      ENDDO      
      
      DO bou = 1,nbou_coarse                    
        nbnds = fbseg_coarse(1,bou)   
        DO i = 1,nbnds   
          nd = fbnds_coarse(i,bou)    
          bou_node(nd) = 1
        ENDDO
      ENDDO   
      
      DO bou = 1,ncbou                    
        nbnds = cbseg_coarse(bou)   
        DO i = 1,nbnds   
          nd = cbnds_coarse(i,bou)    
          bou_node(nd) = 1
        ENDDO
      ENDDO         
      
      ALLOCATE(frc(nn_coarse,2))
      dt = 0.2d0
      
      DO it = 1,20
 nodes: DO n1 = 1,nn_coarse
          
          frc(n1,1) = 0d0
          frc(n1,2) = 0d0
          
          IF (bou_node(n1) == 1) THEN
            CYCLE nodes
          ENDIF                   
          
          DO i = 1,nadjnds(n1)
            n2 = adjnds(i,n1)
            
            px = xy_coarse(1,n1)-xy_coarse(1,n2)
            py = xy_coarse(2,n1)-xy_coarse(2,n2)
            
            lxy = sqrt(px**2+py**2)
            
!             xy(1) = .5d0*(xy_coarse(1,n1)+xy_coarse(1,n2))            
!             xy(2) = .5d0*(xy_coarse(2,n1)+xy_coarse(2,n2))
!             
!             CALL kdtree2_n_nearest(tp=tree_xy,qv=xy,nn=1,results=closest)
!             hpt = hnd(closest(1)%idx)

            hpt = .5d0*(hnd(n1)+hnd(n2))
            l0 = 1.2d0*hpt
            
            frc(n1,1) = frc(n1,1) + max(l0-lxy,0d0)*px/lxy
            frc(n1,2) = frc(n1,2) + max(l0-lxy,0d0)*py/lxy
            
          ENDDO
          
        ENDDO nodes
        
        
        DO n1 = 1,nn_coarse
        
          xy_coarse(1,n1) = xy_coarse(1,n1) + dt*frc(n1,1)
          xy_coarse(2,n1) = xy_coarse(2,n1) + dt*frc(n1,2)
        
        ENDDO
      ENDDO
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Write output
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      

      CALL print_grid_info(grid_file_out,grid_name,ne_coarse,nn_coarse) 
      PRINT "(A,F8.5)", "element reduction factor: ", real(ne,rp)/real(ne_coarse,rp)
      PRINT "(A,F8.5)", "node reduction factor: ", real(nn,rp)/real(nn_coarse,rp)
      PRINT "(A)", ""
      CALL write_header(grid_file_out,grid_name,ne_coarse,nn_coarse)
      CALL write_coords(nn_coarse,xy_coarse,depth_coarse)
      CALL write_connectivity(ne_coarse,ect_coarse,et_coarse,nverts)
      CALL write_open_boundaries(nope,neta_coarse,obseg_coarse,obnds_coarse)
      CALL write_flow_boundaries(nbou_coarse,nvel_coarse,fbseg_coarse,fbnds_coarse)

      END PROGRAM consolidate