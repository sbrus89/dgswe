      PROGRAM refine_grid

      USE globals, ONLY: rp,grid_name,ne,nn,nverts,ect,xy,depth,el_type, &
                         nope,neta,obseg,obnds,nvel,nbou,fbseg,fbnds, &
                         ned,mnepn,epn,nepn,ed_type,recv_edge, &
                         ged2nn,ged2el,ged2led, &
                         nied,iedn,nobed,obedn,nfbed,fbedn,nnfbed,nfbedn,nfbednn                         
                         
      USE grid_file_mod, ONLY: read_header,read_coords,read_connectivity, &
                               read_open_boundaries,read_flow_boundaries,print_grid_info, &                               
                               write_header,write_coords,write_connectivity, &
                               write_open_boundaries,write_flow_boundaries
                               
      USE edge_connectivity_mod, ONLY: elements_per_node,find_edge_pairs, &
                                       find_interior_edges,find_open_edges,find_flow_edges, &
                                       find_element_edges,print_connect_info
                                       
      USE version, ONLY: version_information                                       

      IMPLICIT NONE
      
      INTEGER :: myrank = 0
      INTEGER :: nred = 0
      REAL(rp) :: h0 = 0d0
      
      INTEGER :: ed,i,n,nd,bou
      INTEGER :: ged,el,led,et,nv
      INTEGER :: n1,n2
      INTEGER :: ged1,ged2
      INTEGER :: fed,nfed
      INTEGER :: nnr,ner,netar,nvelr
      INTEGER :: segtype
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: el2ged
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: ect_refine
      INTEGER, DIMENSION(:), ALLOCATABLE :: el_type_refine
      INTEGER, DIMENSION(:), ALLOCATABLE :: el2nn
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: obnds_refine
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: fbnds_refine
      INTEGER, DIMENSION(:), ALLOCATABLE :: obnodes
      INTEGER, DIMENSION(:), ALLOCATABLE :: fbnodes    
      REAL(rp), DIMENSION(:,:), ALLOCATABLE :: xy_refine
      REAL(rp), DIMENSION(:), ALLOCATABLE :: depth_refine
      
      CHARACTER(100) :: grid_file
      CHARACTER(40) :: sha1      
      
      nverts(1) = 3
      nverts(2) = 4
      nverts(3) = 3
      nverts(4) = 4      
      
      
      grid_file = '/home/sbrus/data-drive/galveston_spline_flux_fix/grids/galveston_quad.grd'
      
      
      
      
      
      CALL read_header(myrank,grid_file,grid_name,ne,nn)  
      
      CALL read_coords(nn,xy,depth,h0)     

      CALL read_connectivity(ne,ect,el_type)              
      
      CALL read_open_boundaries(nope,neta,obseg,obnds)      
      
      CALL read_flow_boundaries(nbou,nvel,fbseg,fbnds)
      
      CALL print_grid_info(grid_file,grid_name,ne,nn)      


      
      
      CALL elements_per_node(ne,nn,nverts,el_type,ect,nepn,mnepn,epn) 
      
      CALL find_edge_pairs(ne,nverts,el_type,ect,nepn,epn,ned,ged2el,ged2nn,ged2led)
      
      CALL find_interior_edges(ned,ged2el,nied,iedn,ed_type,recv_edge)
      
      CALL find_open_edges(nope,obseg,obnds,ged2nn,nobed,obedn,ed_type,recv_edge)      
      
      CALL find_flow_edges(nbou,fbseg,fbnds,ged2nn,nnfbed,nfbedn,nfbednn,nfbed,fbedn,recv_edge,ed_type)  
      
      CALL find_element_edges(ne,ned,ged2el,ged2led,el2ged)

      CALL print_connect_info(mnepn,ned,nied,nobed,nfbed,nnfbed,nred)      
      
      
            
            
      ALLOCATE(xy_refine(2,nn+ned+ne))   
      ALLOCATE(depth_refine(nn+ned+ne))
      ALLOCATE(el2nn(ne))      
      DO i = 1,nn
        xy_refine(1,i)  = xy(1,i)
        xy_refine(2,i)  = xy(2,i)        
        depth_refine(i) = depth(i)
      ENDDO
      
      nnr = nn
      DO ged = 1,ned
                    
        n1 = ged2nn(1,ged)
        n2 = ged2nn(2,ged)
        
        n = nn + ged
        nnr = nnr + 1
        
        xy_refine(1,n)  = .5d0*(xy(1,n1)  + xy(1,n2))
        xy_refine(2,n)  = .5d0*(xy(2,n1)  + xy(2,n2))        
        depth_refine(n) = .5d0*(depth(n1) + depth(n2))
      ENDDO
      
      DO el = 1,ne
      
        et = el_type(el)              
        
        IF (et == 2) THEN
!           n = nn + ned + el         
          nnr = nnr+1        
          el2nn(el) = nnr
          
          xy_refine(1,nnr)  = 0d0
          xy_refine(2,nnr)  = 0d0
          depth_refine(nnr) = 0d0
          DO i = 1,nverts(et)
            xy_refine(1,nnr)  = xy_refine(1,nnr)  + .25d0*xy(1,ect(i,el)) 
            xy_refine(2,nnr)  = xy_refine(2,nnr)  + .25d0*xy(2,ect(i,el)) 
            depth_refine(nnr) = depth_refine(nnr) + .25d0*depth(ect(i,el))
          ENDDO          
        ENDIF
        
      ENDDO
      
      
      
      
      
      
      
      ALLOCATE(ect_refine(4,4*ne),el_type_refine(4*ne))
      
      ner = 0
      DO el = 1,ne
      
        et = el_type(el)
        nv = nverts(et)
        IF (et == 1) THEN
        
          DO n = 1,nv         
            ner = ner + 1
            
            ged1 = el2ged(el,mod(n+1,nv)+1)
            ged2 = el2ged(el,mod(n+0,nv)+1)
                        
            ect_refine(1,ner) = ect(n,el)
            ect_refine(2,ner) = nn + ged1
            ect_refine(3,ner) = nn + ged2
            
            el_type_refine(ner) = 1
          ENDDO
          
          ner = ner + 1
          DO n = 1,nv
            ect_refine(n,ner) = nn + el2ged(el,n)
          ENDDO
          el_type_refine(ner) = 1
          
        ELSE IF (et == 2) THEN
        
          nd = el2nn(el)
        
          DO n = 1,nv          
            ner = ner + 1
            
            ged1 = el2ged(el,mod(n+2,nv)+1)
            ged2 = el2ged(el,mod(n+1,nv)+1)
            
            ect_refine(1,ner) = ect(n,el)
            ect_refine(2,ner) = nn + ged1
            ect_refine(3,ner) = nd
            ect_refine(4,ner) = nn + ged2
            
            el_type_refine(ner) = 2
          ENDDO       
          
        ENDIF
              
      ENDDO
      
      
      
      
      
      
      ALLOCATE(obnodes(2*neta),obnds_refine(2*neta,nope))
            
      ed = 0
      netar = 0 
      DO bou = 1,nope
        DO nd = 1,obseg(bou)-1
          netar = netar + 1        
          obnds_refine(2*nd-1,bou) = obnds(nd,bou)
          
          netar = netar + 1          
          ed = ed + 1                   
          ged = obedn(ed)
          obnds_refine(2*nd,bou) = nn + ged
        ENDDO         

        nd = obseg(bou)
        obseg(bou) = 2*obseg(bou)-1
        
        netar = netar + 1        
        obnds_refine(2*nd-1,bou) = obnds(nd,bou)
      ENDDO
      
      
      ALLOCATE(fbnodes(2*nvel),fbnds_refine(2*nvel,nbou))
            
      fed = 0
      nfed = 0
      nvelr = 0
      DO bou = 1,nbou
        segtype = fbseg(2,bou)
        DO nd = 1,fbseg(1,bou)-1
        
          nvelr = nvelr + 1
          fbnds_refine(2*nd-1,bou) = fbnds(nd,bou)
          
          IF( segtype == 0 .OR. segtype == 10 .OR. segtype == 20  .OR. &   ! land boundaries
              segtype == 1 .OR. segtype == 11 .OR. segtype == 21 ) THEN    ! island boundaries
              
              nfed = nfed + 1          
              ged = nfbedn(nfed)
          ENDIF
          
          IF ( segtype == 2 .OR. segtype == 12 .OR. segtype == 22 ) THEN             
               
               fed = fed + 1
               ged = fbedn(fed)
          ENDIF          
          
          nvelr = nvelr + 1
          fbnds_refine(2*nd,bou) = nn + ged
        ENDDO

        nd = fbseg(1,bou)
        fbseg(1,bou) = 2*fbseg(1,bou)-1
        IF( segtype == 1 .OR. segtype == 11 .OR. segtype == 21 ) THEN    ! island boundaries
          fbseg(1,bou) = fbseg(1,bou)-1
        ENDIF
        
        nvelr = nvelr + 1        
        fbnds_refine(2*nd-1,bou) = fbnds(nd,bou)
      ENDDO
      
      
      
      CALL write_header('fort.14_refine',grid_name,ner,nnr)
      
      CALL write_coords(nnr,xy_refine,depth_refine)
      
      CALL write_connectivity(ner,ect_refine,el_type_refine,nverts)
      
      CALL write_open_boundaries(nope,netar,obseg,obnds_refine)
      
      CALL write_flow_boundaries(nbou,nvelr,fbseg,fbnds_refine)
      
      
      OPEN(UNIT=14, FILE="fort.14_refine", POSITION="APPEND")
      WRITE(14,"(A)") "4x refinement of base grid "//grid_file
      WRITE(14,"(A)") "base grid SHA "//sha1(grid_file,"./")
      WRITE(14,"(A)") ""
      CALL version_information(14)
      
      WRITE(14,"(A)") "-----------------------------------------------------------------------"           
      CLOSE(14)

      
      
      END PROGRAM refine_grid