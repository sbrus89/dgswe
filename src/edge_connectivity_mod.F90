      MODULE edge_connectivity_mod
      
      USE quit, ONLY: abort
      
      IMPLICIT NONE
      
      INTEGER :: el
      INTEGER :: nd
      INTEGER :: alloc_status

      CONTAINS
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       

      
      SUBROUTINE elements_per_node(ne,nn,nverts,el_type,vct,nepn,mnepn,epn)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: ne
      INTEGER, INTENT(IN) :: nn
      INTEGER, DIMENSION(:), INTENT(IN) :: nverts
      INTEGER, DIMENSION(:), INTENT(IN) :: el_type
      INTEGER, DIMENSION(:,:), INTENT(IN) :: vct
      INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: nepn
      INTEGER, INTENT(OUT) :: mnepn
      INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: epn
      
      
      INTEGER :: nvert  
      INTEGER :: n           
      
      

      ! count the number of elements per node
      
      ALLOCATE(nepn(nn), STAT=alloc_status)
      IF (alloc_status /= 0) THEN
        PRINT*, "Allocation error"
        CALL abort()
      ENDIF         

      nepn(:) = 0
      DO el = 1,ne
        nvert = nverts(el_type(el))
        DO nd = 1,nvert
          n = vct(nd,el)
          nepn(n) = nepn(n) + 1
        ENDDO        
      ENDDO

      
      ! find the maximum number of elements per node

      mnepn = 0
      DO nd = 1,nn
        IF(nepn(nd) > mnepn) THEN
          mnepn = nepn(nd)
        ENDIF
      ENDDO
      
      

      ! find the elements associated with each node
      
      ALLOCATE(epn(mnepn,nn), STAT=alloc_status)
      IF (alloc_status /= 0) THEN
        PRINT*, "Allocation error"
        CALL abort()
      ENDIF         

      nepn(:) = 0
      DO el = 1,ne
        nvert = nverts(el_type(el))
        DO nd = 1,nvert
          n = vct(nd,el)
          nepn(n) = nepn(n) + 1
          epn(nepn(n),n) = el
        ENDDO
      ENDDO      
      
      RETURN
      END SUBROUTINE elements_per_node
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 


      SUBROUTINE find_edge_pairs(ne,nverts,el_type,vct,nepn,epn,ned,ged2el,ged2nn,ged2led)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: ne
      INTEGER, DIMENSION(:), INTENT(IN) :: nverts
      INTEGER, DIMENSION(:), INTENT(IN) :: el_type
      INTEGER, DIMENSION(:,:), INTENT(IN) :: vct
      INTEGER, DIMENSION(:), INTENT(IN) :: nepn
      INTEGER, DIMENSION(:,:), INTENT(IN) :: epn
      INTEGER, INTENT(OUT) :: ned      
      INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: ged2el
      INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: ged2nn
      INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: ged2led
      
      
      INTEGER :: el1,el2
      INTEGER :: led1,led2
      INTEGER :: nvert1,nvert2
      INTEGER :: n1ed1,n2ed1
      INTEGER :: n1ed2,n2ed2
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ged2nn_temp, ged2el_temp, ged2led_temp
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: edflag
      
      
      ALLOCATE(ged2nn_temp(2,4*ne),ged2el_temp(2,4*ne),ged2led_temp(2,4*ne), STAT = alloc_status)      
      IF (alloc_status /= 0) THEN
        PRINT*, "Allocation error"
        CALL abort()
      ENDIF             
      
      ALLOCATE(edflag(4,ne), STAT = alloc_status)
      IF (alloc_status /= 0) THEN
        PRINT*, "Allocation error"
        CALL abort()
      ENDIF       
      
      
      ned = 0
      edflag = 0
      ged2nn_temp = 0
      ged2led_temp = 0
      ged2el_temp = 0
      
   elem1: DO el1 = 1,ne ! loop through trial elements
   
            nvert1 = nverts(el_type(el1))

 local_ed1: DO led1 = 1,nvert1 ! loop through trial edges

              IF(edflag(led1,el1) == 1) THEN ! skip if edge has already been flagged
                CYCLE local_ed1
              ENDIF

              ned = ned + 1 ! increment edge number

              n1ed1 = vct(mod(led1+0,nvert1)+1,el1) ! find nodes on trial edge
              n2ed1 = vct(mod(led1+1,nvert1)+1,el1)
              
              ged2nn_temp(1,ned) = n1ed1 ! set nodes on global edge # ned
              ged2nn_temp(2,ned) = n2ed1

              ged2led_temp(1,ned) = led1 ! set local edge number of first element sharing the edge
              ged2el_temp(1,ned) = el1   ! set the first element that shares the edge

              edflag(led1,el1) = 1  ! flag the edge so it is not repeated

       elem2: DO el = 1,nepn(n1ed1) ! loop through test elements (only those that contain node n1ed1)

                el2 = epn(el,n1ed1) ! choose a test element that contains node n1ed1 
                
                IF(el2 == el1) THEN ! skip if the test element is the same as the trial element
                  CYCLE elem2
                ENDIF
                
                nvert2 = nverts(el_type(el2))

     local_ed2: DO led2 = 1,nvert2 ! loop through local test edge numbers
                  
                  n1ed2 = vct(MOD(led2+0,nvert2)+1,el2) ! find nodes on test edge
                  n2ed2 = vct(MOD(led2+1,nvert2)+1,el2)

                  IF(((n1ed1 == n1ed2) .AND. (n2ed1 == n2ed2)) .OR. & ! check if nodes on trial edge matches test edge
                     ((n1ed1 == n2ed2) .AND. (n2ed1 == n1ed2))) THEN

                     ged2led_temp(2,ned) = led2 ! set local edge number of second element sharing the edge
                     ged2el_temp(2,ned) = el2   ! set the second element that shares the edge
                     
                     edflag(led2,el2) = 1       ! flag the edge so it is not repeated

                     EXIT elem2          
                  ENDIF

                ENDDO local_ed2

              ENDDO elem2

            ENDDO local_ed1
      
          ENDDO elem1

      ALLOCATE(ged2nn(2,ned),ged2el(2,ned),ged2led(2,ned), STAT=alloc_status)
      IF (alloc_status /= 0) THEN
        PRINT*, "Allocation error"
        CALL abort()
      ENDIF          

      ged2nn(:,1:ned) = ged2nn_temp(:,1:ned)
      ged2el(:,1:ned) = ged2el_temp(:,1:ned)
      ged2led(:,1:ned) = ged2led_temp(:,1:ned)      
      
      RETURN
      END SUBROUTINE find_edge_pairs

      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 


      SUBROUTINE find_interior_edges(ned,ged2el,nbed,bedn,nied,iedn,ed_type,recv_edge)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: ned
      INTEGER, DIMENSION(:,:), INTENT(IN) :: ged2el
      INTEGER, INTENT(OUT) :: nbed
      INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: bedn
      INTEGER, INTENT(OUT) :: nied
      INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: iedn
      INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: ed_type
      INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: recv_edge
                  
      INTEGER :: ged
      INTEGER :: el1,el2
      INTEGER, DIMENSION(:), ALLOCATABLE :: ied_temp
      
      ALLOCATE(ied_temp(ned), STAT = alloc_status)
      IF (alloc_status /= 0) THEN
        PRINT*, "Allocation error"
        CALL abort()
      ENDIF       
      
      ALLOCATE(bedn(ned),recv_edge(ned),ed_type(ned), STAT = alloc_status)
      IF (alloc_status /= 0) THEN
        PRINT*, "Allocation error"
        CALL abort()
      ENDIF      
      
      recv_edge = 1 
      ed_type = 0

      nied = 0
      nbed = 0
      ied_temp(:) = 0
      DO ged = 1,ned
        el1 = ged2el(1,ged)
        el2 = ged2el(2,ged)
        IF ((el1 /= 0) .AND. (el2 /= 0)) THEN
          nied = nied + 1
          ied_temp(nied) = ged
          recv_edge(ged) = 0   
          ed_type(ged) = 0
        ELSE 
          nbed = nbed + 1
          bedn(nbed) = ged
        ENDIF        
      ENDDO      
      
      ALLOCATE(iedn(nied), STAT=alloc_status)
      IF (alloc_status /= 0) THEN
        PRINT*, "Allocation error"
        CALL abort()
      ENDIF        
      
      iedn(1:nied) = ied_temp(1:nied)         
      
      RETURN
      END SUBROUTINE find_interior_edges
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE find_open_edges(nope,obseg,obnds,ged2nn,nbed,bedn,nobed,obedn,ed_type,recv_edge)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: nope
      INTEGER, DIMENSION(:), INTENT(IN) :: obseg
      INTEGER, DIMENSION(:,:), INTENT(IN) :: obnds
      INTEGER, DIMENSION(:,:), INTENT(IN) :: ged2nn
      INTEGER, INTENT(IN) :: nbed
      INTEGER, DIMENSION(:), INTENT(IN) :: bedn
      INTEGER, INTENT(OUT) :: nobed
      INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: obedn
      INTEGER, DIMENSION(:), INTENT(INOUT) :: ed_type
      INTEGER, DIMENSION(:), INTENT(INOUT) :: recv_edge
      
      
      INTEGER :: seg,nd,ed,ged
      INTEGER :: n1bed,n2bed
      INTEGER :: n1ed2,n2ed2
      INTEGER, DIMENSION(:), ALLOCATABLE :: obedn_temp
      
      ALLOCATE(obedn_temp(nbed), STAT = alloc_status)
      IF (alloc_status /= 0) THEN
        PRINT*, "Allocation error"
        CALL abort()
      ENDIF        
      
      nobed = 0 
      obedn_temp = 0
      
      DO seg = 1,nope
        DO nd = 1,obseg(seg)-1
        
          n1bed = obnds(nd,seg)
          n2bed = obnds(nd+1,seg)
          
  edges1: DO ed = 1,nbed
  
            ged = bedn(ed) 
            n1ed2 = ged2nn(1,ged)
            n2ed2 = ged2nn(2,ged)
            
            IF(((n1ed2 == n1bed).AND.(n2ed2 == n2bed)).OR. &
               ((n1ed2 == n2bed).AND.(n2ed2 == n1bed))) THEN
               
              nobed = nobed + 1
              obedn_temp(nobed) = ged
              recv_edge(ged) = 0
              ed_type(ged) = 1
              
              EXIT edges1
            ENDIF
            
          ENDDO edges1
        ENDDO
      ENDDO      
      
      ALLOCATE(obedn(nobed), STAT=alloc_status)
      IF (alloc_status /= 0) THEN
        PRINT*, "Allocation error"
        CALL abort()
      ENDIF          
      
      obedn(1:nobed) = obedn_temp(1:nobed)
      
      RETURN
      END SUBROUTINE find_open_edges

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE find_flow_edges(nbou,fbseg,fbnds,ged2nn,nbed,bedn,nnfbed,nfbedn,nfbednn,nfbed,fbedn,recv_edge,ed_type)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: nbou
      INTEGER, DIMENSION(:,:), INTENT(IN) :: fbseg
      INTEGER, DIMENSION(:,:), INTENT(IN) :: fbnds
      INTEGER, DIMENSION(:,:), INTENT(IN) :: ged2nn
      INTEGER, INTENT(IN) :: nbed
      INTEGER, DIMENSION(:), INTENT(IN) :: bedn
      INTEGER, INTENT(OUT) :: nnfbed
      INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: nfbedn
      INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: nfbednn
      INTEGER, INTENT(OUT) :: nfbed
      INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: fbedn
      INTEGER, DIMENSION(:), INTENT(INOUT) :: recv_edge
      INTEGER, DIMENSION(:), INTENT(INOUT) :: ed_type      
      
      INTEGER :: seg,nd,ed
      INTEGER :: ged,segtype
      INTEGER :: n1bed,n2bed
      INTEGER :: n1ed2,n2ed2
      INTEGER :: found
      INTEGER, DIMENSION(:), ALLOCATABLE :: nfbedn_temp
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: nfbednn_temp
      INTEGER, DIMENSION(:), ALLOCATABLE :: fbedn_temp
      
      
      ALLOCATE(nfbedn_temp(nbed),nfbednn_temp(nbed,2),fbedn_temp(nbed), STAT=alloc_status)
      IF (alloc_status /= 0) THEN
        PRINT*, "Allocation error"
        CALL abort()
      ENDIF             
      
      
      nnfbed = 0 ! # no normal flow edge
      nfbedn_temp = 0
      nfbednn_temp = 0
      
      nfbed = 0 ! # flow specifed edges
      fbedn_temp = 0

      found = 0

      DO seg = 1,nbou
      
        segtype = fbseg(2,seg)
              
        DO nd = 1,fbseg(1,seg)-1
        
          n1bed = fbnds(nd,seg)
          n2bed = fbnds(nd+1,seg)
          found = 0 
          
  edges2: DO ed = 1,nbed
  
            ged = bedn(ed)
            n1ed2 = ged2nn(1,ged)
            n2ed2 = ged2nn(2,ged)
            
            IF(((n1ed2 == n1bed).AND.(n2ed2 == n2bed)).OR. &
               ((n1ed2 == n2bed).AND.(n2ed2 == n1bed))) THEN

              ! no normal flow edges
              IF( segtype == 0 .OR. segtype == 10 .OR. segtype == 20  .OR. &   ! land boundaries
                  segtype == 1 .OR. segtype == 11 .OR. segtype == 21 ) THEN    ! island boundaries
                  
                nnfbed = nnfbed + 1
                nfbedn_temp(nnfbed) = ged
                nfbednn_temp(nnfbed,1) = seg
                nfbednn_temp(nnfbed,2) = n1bed
                recv_edge(ged) = 0
                ed_type(ged) = 10
                found = 1
                
                EXIT edges2               
              ENDIF

              ! specified normal flow edges
              IF ( segtype == 2 .OR. segtype == 12 .OR. segtype == 22 ) THEN
              
                nfbed = nfbed + 1
                fbedn_temp(nfbed) = ged
                recv_edge(ged) = 0
                ed_type(ged) = 12
                found = 1
                
                EXIT edges2
              ENDIF

            ENDIF
          ENDDO edges2
          
          IF (found == 0) THEN
            PRINT*, seg, nd, fbseg(1,seg), segtype          
            PRINT "(A)", "  edge not found"
          ENDIF
          
        ENDDO
      ENDDO
      
      ALLOCATE(nfbedn(nnfbed),nfbednn(nnfbed,2),fbedn(nfbed), STAT=alloc_status)
      
      
      nfbedn(1:nnfbed) = nfbedn_temp(1:nnfbed)
      nfbednn(1:nnfbed,1:2) = nfbednn_temp(1:nnfbed,1:2)
      fbedn(1:nfbed) = fbedn_temp(1:nfbed)

      END SUBROUTINE find_flow_edges
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE find_recieve_edges(ned,recv_edge,nred,redn,el_type)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: ned
      INTEGER, DIMENSION(:), INTENT(IN) :: recv_edge
      INTEGER, INTENT(OUT) :: nred
      INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: redn
      INTEGER, DIMENSION(:), INTENT(INOUT) :: el_type
      
      INTEGER :: ged                  

      nred = 0
      
#ifdef CMPI      
      DO ged = 1,ned
        IF(recv_edge(ged) == 1) THEN
          nred = nred + 1
        ENDIF
      ENDDO
      
      ALLOCATE(redn(nred))
      
      nred = 0
      DO ged = 1,ned
        IF(recv_edge(ged) == 1) THEN
          nred = nred + 1
          redn(nred) = ged
          ed_type(ged) = -1
        ENDIF
      ENDDO
#endif         
      
      RETURN
      END SUBROUTINE find_recieve_edges

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE print_connect_info(mnepn,ned,nied,nobed,nfbed,nnfbed,nred)
      
      IMPLICIT NONE      
      
      INTEGER, INTENT(IN) :: mnepn
      INTEGER, INTENT(IN) :: ned
      INTEGER, INTENT(IN) :: nied
      INTEGER, INTENT(IN) :: nobed
      INTEGER, INTENT(IN) :: nfbed
      INTEGER, INTENT(IN) :: nnfbed
      INTEGER, INTENT(IN) :: nred
      
      
      
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", "       Edge Connectivity Information         "
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", " "        
      PRINT "(A,I7)", '   maximum elements per node:', mnepn
      PRINT "(A)", ' '            
      PRINT "(A,I7)", '   number of total edges:', ned
      PRINT "(A)", ' '  
      PRINT "(A,I7)", '   number of interior edges:', nied
      PRINT "(A)", ' '        
      PRINT "(A,I7)", '   number of open boundary edges:', nobed
      PRINT "(A)", ' '                
      PRINT "(A,I7)", '   number of specified normal boundary edges:', nfbed
      PRINT "(A,I7)", '   number of no normal flow boundary edges:', nnfbed
      PRINT "(A)", ' '        
      PRINT "(A,I7)", '   number of recieve edges:', nred
      PRINT "(A)", ' ' 
      PRINT "(A,I7)", '   number of missing edges:',ned-(nied+nobed+nfbed+nnfbed+nred)
      PRINT "(A)", ' '         
      
      
      RETURN
      END SUBROUTINE print_connect_info

      END MODULE edge_connectivity_mod