      PROGRAM consolidate

      USE globals, ONLY: rp,ne,nn,nverts,ect,xy,depth,el_type, &
                         nope,neta,obseg,obnds,nvel,nbou,fbseg,fbnds,grid_name, &
                         nepn,mnepn,epn,ned,ged2el,ged2nn,ged2led
      USE grid_file_mod
      USE edge_connectivity_mod      
      USE triangulation, ONLY: delaunay_triangulation      
      USE kdtree2_module      

      IMPLICIT NONE
      
      CHARACTER(200) :: fine_file_in 
      CHARACTER(200) :: fine_file_out      
      CHARACTER(200) :: coarse_file_out
      CHARACTER(100) :: grid_name_coarse
      
      INTEGER :: i,j,k,l,m
      INTEGER :: ed,nd,bou,el,el1,el2
      INTEGER :: n1,n2,n3
      INTEGER :: nodes_kept
      INTEGER :: nbnds
      INTEGER :: bou_type
      INTEGER :: ndcnt
      INTEGER, ALLOCATABLE, DIMENSION(:) :: nadjnds
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: adjnds   
      INTEGER, ALLOCATABLE, DIMENSION(:) :: keep_node
      INTEGER, ALLOCATABLE, DIMENSION(:) :: nd_fine2coarse
      INTEGER, ALLOCATABLE, DIMENSION(:) :: nd_coarse2fine      
      INTEGER, ALLOCATABLE, DIMENSION(:) :: bou_node
      INTEGER, ALLOCATABLE, DIMENSION(:) :: bou_node_coarse 
      INTEGER, ALLOCATABLE, DIMENSION(:) :: node_flag      
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
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: orient      
      INTEGER :: ncbou
      INTEGER :: ncnds
      INTEGER, ALLOCATABLE, DIMENSION(:) :: cbnd_flag
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
      INTEGER, PARAMETER :: cbou_type = 30
      
      TYPE(kdtree2), POINTER :: tree_xy      
      TYPE(kdtree2_result), ALLOCATABLE, DIMENSION(:) :: closest     
      INTEGER :: srchdp
      REAL(rp) :: len_avg
      REAL(rp), DIMENSION(:), ALLOCATABLE :: hnd
      INTEGER :: it,retri
      REAL(rp) :: px,py
      REAL(rp) :: lxy
      REAL(rp) :: hpt,l0
      REAL(rp) :: dt
      REAL(rp), DIMENSION(:,:), ALLOCATABLE ::  frc    
      
      
!       fine_file_in = "/home/sbrus/data-drive/galveston_SL18/grid_dev/v26_cart/galveston_SL18_cart.grd"
!       coarse_file_out = "/home/sbrus/data-drive/galveston_SL18/grid_dev/v26_cart/coarse.grd"
!       fine_file_out = "/home/sbrus/data-drive/galveston_SL18/grid_dev/v26_cart/galveston_SL18_cart_no30.grd"      

!       fine_file_in = "/home/sbrus/data-drive/galveston_SL18/grid_dev/v26_cart/coarse.grd"
!       coarse_file_out = "/home/sbrus/data-drive/galveston_SL18/grid_dev/v26_cart/coarse_x2.grd"
!       fine_file_out = "/home/sbrus/data-drive/galveston_SL18/grid_dev/v26_cart/coarse_no30.grd" 

!       fine_file_in = "./galveston_SL18_cart.grd"
!       coarse_file_out = "./coarse.grd"
      
      fine_file_in = "./coarse_x2_reduced_constrained.grd"
      coarse_file_out = "./coarse_x2_reduced_x2.grd"      

         
      nverts(1) = 3
      nverts(2) = 4
      nverts(3) = 3
      nverts(4) = 4         
         

      CALL read_header(0,fine_file_in,grid_name,ne,nn)        
      CALL read_coords(nn,xy,depth)
      CALL read_connectivity(ne,ect,el_type)                 
      CALL read_open_boundaries(nope,neta,obseg,obnds)            
      CALL read_flow_boundaries(nbou,nvel,fbseg,fbnds)
      CALL print_grid_info(fine_file_in,grid_name,ne,nn)          
      
      CALL elements_per_node(ne,nn,nverts,el_type,ect,nepn,mnepn,epn)       
      CALL find_edge_pairs(ne,nverts,el_type,ect,nepn,epn,ned,ged2el,ged2nn,ged2led)   
!       CALL boundary_orientation(nbou,fbseg,fbnds,ne,ect,xy,nepn,epn,orient)                    
      CALL find_adjacent_nodes(nn,mnepn,ned,ged2nn,nadjnds,adjnds)
     
      
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Eliminate node from input (fine) mesh
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      

      
      ALLOCATE(keep_node(nn))
      ALLOCATE(bou_node(nn))
      ALLOCATE(fbseg_orient(nbou))
      ALLOCATE(cbnd_flag(nn))
      ALLOCATE(node_flag(nn))
      keep_node = -1
      nodes_kept = 0    
      bou_node = 0    
      node_flag = 0
      
      OPEN(unit=123,file="node.list")
      READ(123,*) ndcnt
      DO i = 1,ndcnt
        READ(123,*) nd
        node_flag(nd) = 1
      ENDDO
      CLOSE(123)  
      

      
      
      
      ! Find feature constraint nodes
     
      cbnd_flag = 0
      DO bou = 1,nbou
        nbnds = fbseg(1,bou)
        bou_type = fbseg(2,bou)
        
        IF (bou_type /= cbou_type) THEN
          CYCLE
        ENDIF        
        
        DO j = 1,nbnds
          nd = fbnds(j,bou)
          cbnd_flag(nd) = 1
        ENDDO        
      ENDDO       
            
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
        IF (bou_type /= cbou_type) THEN
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
        
        IF (bou_type == cbou_type) THEN          ! skip feature constraint nodestrings
          CYCLE fbnd
        ENDIF                
        
        nbou_coarse = nbou_coarse + 1
        
        IF (nbnds <= 9) THEN              ! keep small island boundaries
          DO j = 1,nbnds
            nd = fbnds(j,bou)
            IF (keep_node(nd) == -1) THEN
              keep_node(nd) = 1
              nodes_kept = nodes_kept + 1
              PRINT*, "keeping node on small island", nd
            ENDIF
          ENDDO
        ELSE
        
          nd = fbnds(1,bou)               ! keep first node
          keep_node(nd) = 1
          DO i = 1,nadjnds(nd)
            IF (cbnd_flag(adjnds(i,nd)) /= 1) THEN          
              keep_node(adjnds(i,nd)) = 0
            ENDIF
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
                  IF (cbnd_flag(adjnds(i,nd)) /= 1) THEN          
                    keep_node(adjnds(i,nd)) = 0
                  ENDIF
                ENDIF
              ENDDO
              nodes_kept = nodes_kept + 1        
            ENDIF
          ENDDO
        
          nd = fbnds(nbnds,bou)           ! keep last node
          keep_node(nd) = 1
          DO i = 1,nadjnds(nd)
            IF (cbnd_flag(adjnds(i,nd)) /= 1) THEN          
              keep_node(adjnds(i,nd)) = 0
            ENDIF
          ENDDO 
          nodes_kept = nodes_kept + 1   
          
        ENDIF
              
      ENDDO fbnd
      
      nvel_coarse = nodes_kept - neta_coarse
      
      
      ! Interior nodes
      
      DO nd = 1,nn
        IF (keep_node(nd) == -1 .and. cbnd_flag(nd) /= 1) THEN
          keep_node(nd) = 1
          DO i = 1,nadjnds(nd)
            IF (keep_node(adjnds(i,nd)) /= 1 .and. cbnd_flag(adjnds(i,nd)) /= 1) THEN
              keep_node(adjnds(i,nd)) = 0          
            ENDIF
          ENDDO
          nodes_kept = nodes_kept + 1                  
        ENDIF
      ENDDO               
      

      ! Feature constraint nodestrings
      
      ncbou = 0
 cbnd:DO bou = 1,nbou
        nbnds = fbseg(1,bou)
        bou_type = fbseg(2,bou)
        
        IF (bou_type /= cbou_type) THEN
          CYCLE
        ENDIF
        
        ncbou = ncbou + 1
        
        
        nd = fbnds(1,bou)               ! keep first node
        keep_node(nd) = 1
        DO i = 1,nadjnds(nd)
          IF (cbnd_flag(adjnds(i,nd)) == 1) THEN
            IF (adjnds(i,nd) == fbnds(j+1,bou)) THEN
              keep_node(adjnds(i,nd)) = 0
            ENDIF
          ELSE IF (bou_node(adjnds(i,nd)) /= 1) THEN
            keep_node(adjnds(i,nd)) = 0
          ENDIF
        ENDDO 
        nodes_kept = nodes_kept + 1           

        DO j = 2,nbnds-1
          nd = fbnds(j,bou)
     
          IF (keep_node(nd) == -1) THEN   
            keep_node(nd) = 1
            DO i = 1,nadjnds(nd)
              IF (cbnd_flag(adjnds(i,nd)) == 1 ) THEN
                IF (adjnds(i,nd) == fbnds(j-1,bou) .or. adjnds(i,nd) == fbnds(j+1,bou)) THEN              
                  keep_node(adjnds(i,nd)) = 0 
                ENDIF
              ELSE IF (bou_node(adjnds(i,nd)) /= 1) THEN
                keep_node(adjnds(i,nd)) = 0               
              ENDIF
            ENDDO
            nodes_kept = nodes_kept + 1
          ENDIF          
        ENDDO
        
        nd = fbnds(nbnds,bou)               ! keep last node
        keep_node(nd) = 1
        DO i = 1,nadjnds(nd)
          IF (cbnd_flag(adjnds(i,nd)) == 1) THEN
            IF (adjnds(i,nd) == fbnds(j-1,bou)) THEN
              keep_node(adjnds(i,nd)) = 0
            ENDIF
          ELSE IF (bou_node(adjnds(i,nd)) /= 1) THEN
            keep_node(adjnds(i,nd)) = 0
          ENDIF
        ENDDO 
        nodes_kept = nodes_kept + 1           
        
      ENDDO cbnd
      
      
      
      DO i = 1,nn
        IF (node_flag(i) == 0 .and. keep_node(i) == 0) THEN                
          keep_node(i) = 1
          nodes_kept = nodes_kept + 1
        ENDIF
      ENDDO
      nvel_coarse = nvel
            
      
      
      i = 0
      DO nd = 1,nn
        IF (keep_node(nd) < 0) THEN
          i = i+1
        ENDIF
      ENDDO
      
      IF (i > 0) THEN
        PRINT*, i, " unassigned nodes"
        STOP
      ENDIF
            
  
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Keep track of node infromation for coarse mesh
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!           
      
      ALLOCATE(xy_coarse(2,nodes_kept))
      ALLOCATE(depth_coarse(nodes_kept))
      ALLOCATE(bou_node_coarse(nodes_kept))
      ALLOCATE(nd_fine2coarse(nn),nd_coarse2fine(nn))
      
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
          nd_coarse2fine(nn_coarse) = nd
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
      ALLOCATE(cbseg_coarse(ncbou))
      ALLOCATE(cbnds_coarse(nvel_coarse,ncbou))
      
      nbou_coarse = 0
      ncbou = 0 
      DO bou = 1,nbou
        bou_type = fbseg(2,bou)
        nbnds = 0        
        IF (bou_type /= cbou_type) THEN
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
      CALL delaunay_triangulation(nn_coarse,xy_coarse(1,:),xy_coarse(2,:),nbou_coarse,fbseg_coarse,fbnds_coarse, &
                                                                          ne_coarse,ect_coarse, &
                                                                          ncbou,cbseg_coarse,cbnds_coarse) 

      ALLOCATE(et_coarse(3*nn_coarse))
      DO el = 1,3*nn_coarse
        et_coarse(el) = 1
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
      
      DO retri = 1,4
      DO it = 1,10
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
      
      CALL delaunay_triangulation(nn_coarse,xy_coarse(1,:),xy_coarse(2,:),nbou_coarse,fbseg_coarse,fbnds_coarse, &
                                                                          ne_coarse,ect_coarse, &
                                                                          ncbou,cbseg_coarse,cbnds_coarse) 
                                                                          
      CALL elements_per_node(ne_coarse,nn_coarse,nverts,et_coarse,ect_coarse,nepn,mnepn,epn) 
      CALL find_edge_pairs(ne_coarse,nverts,et_coarse,ect_coarse,nepn,epn,ned,ged2el,ged2nn,ged2led) 
      CALL find_adjacent_nodes(nn_coarse,mnepn,ned,ged2nn,nadjnds,adjnds)
      
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
      
      ENDDO
      
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Write output
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      
      ncnds = 0                                        ! enclude feature constraints in coarse mesh so they are there for 
      DO bou = 1,ncbou                                 ! subsequent refinements
        nbnds = cbseg_coarse(bou)
        fbseg_coarse(1,nbou_coarse+bou) = nbnds
        fbseg_coarse(2,nbou_coarse+bou) = cbou_type
        DO i = 1,nbnds
          fbnds_coarse(i,nbou_coarse+bou) = cbnds_coarse(i,bou)
          ncnds = ncnds + 1
        ENDDO
      ENDDO
      nbou_coarse = nbou_coarse + ncbou
      nvel_coarse = nvel_coarse + ncnds
             
             
             
             
             
      grid_name_coarse = grid_name // " coarse"
      
      CALL print_grid_info(coarse_file_out,grid_name_coarse,ne_coarse,nn_coarse) 
      PRINT "(A,F8.5)", "element reduction factor: ", real(ne,rp)/real(ne_coarse,rp)
      PRINT "(A,F8.5)", "node reduction factor: ", real(nn,rp)/real(nn_coarse,rp)
      PRINT "(A)", ""
      CALL write_header(coarse_file_out,grid_name_coarse,ne_coarse,nn_coarse)
      CALL write_coords(nn_coarse,xy_coarse,depth_coarse)
      CALL write_connectivity(ne_coarse,ect_coarse,et_coarse,nverts)
      CALL write_open_boundaries(nope,neta_coarse,obseg_coarse,obnds_coarse)
      CALL write_flow_boundaries(nbou_coarse,nvel_coarse,fbseg_coarse,fbnds_coarse)
      
!       ! Write out fine grid without type 30 boundaries if present 
!       IF (ncbou > 0) THEN
!         CALL write_header(fine_file_out,grid_name,ne,nn)
!         CALL write_coords(nn,xy,depth)
!         CALL write_connectivity(ne,ect,el_type,nverts)
!         CALL write_open_boundaries(nope,neta,obseg,obnds)
!         
!         bou = nbou
!         nbou = 0
!         nvel = 0        
!         DO i = 1,bou
!         
!           bou_type = fbseg(2,i)
!           nbnds = fbseg(1,i)          
!           
!           IF (bou_type == cbou_type) THEN
!             CYCLE
!           ENDIF
!           
!           nbou = nbou + 1
!           
! !           IF (bou_type == 1 .OR. bou_type == 11 .OR. bou_type == 21) THEN
! !             nbnds = nbnds - 1
! !           ENDIF
!           
!           DO j = 1,nbnds
!             nvel = nvel + 1
!             fbnds(j,nbou) = fbnds(j,i)
!           ENDDO
!           fbseg(1,nbou) = nbnds
!           fbseg(2,nbou) = bou_type
!                     
!         ENDDO
!         
!         CALL write_flow_boundaries(nbou,nvel,fbseg,fbnds)        
!       ENDIF

      END PROGRAM consolidate