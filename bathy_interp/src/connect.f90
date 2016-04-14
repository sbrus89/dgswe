      SUBROUTINE connect(mesh)

      USE globals, ONLY: rp,grid
                         

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

      ALLOCATE(mesh%nepn(nn),STAT = alloc_status)
      IF(alloc_status /= 0) THEN
        PRINT*, 'Allocation error: nepn'
      ENDIF 

      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", "       Edge Connectivity Information         "
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", " "

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! count the number of elements per node
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      PRINT "(A)", 'counting elements per node'
      PRINT "(A)", ' '

      mesh%nepn(:) = 0
      DO el = 1,ne
        nvert = mesh%nverts(mesh%el_type(el))
        DO nd = 1,nvert
          n1 = mesh%ect(nd,el)
          mesh%nepn(n1) = mesh%nepn(n1) + 1
        ENDDO
        
      ENDDO

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! find the maximum number of elements per node
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      PRINT "(A)", 'finding maximum elements per node'

      mnepn = maxval(mesh%nepn)
      mesh%mnepn = mnepn

      PRINT "(A,I7)", '   maximum elements per node:', mnepn
      PRINT "(A)", ' '

      ALLOCATE(mesh%epn(mnepn,nn),STAT = alloc_status)
      IF(alloc_status /= 0) THEN
        PRINT*, 'Allocation error: epn'
      ENDIF 

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! find the elements associated with each node
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      PRINT "(A)", 'finding elements associated with each node'
      PRINT "(A)", ' '

      mesh%nepn(:) = 0
      DO el = 1,ne
        nvert = mesh%nverts(mesh%el_type(el))
        DO nd = 1,nvert
          n1 = mesh%ect(nd,el)
          mesh%nepn(n1) = mesh%nepn(n1) + 1
          mesh%epn(mesh%nepn(n1),n1) = el
        ENDDO

      ENDDO

      ALLOCATE(edflag(4,ne),STAT = alloc_status)
      IF(alloc_status /= 0) THEN
        PRINT*, 'Allocation error: edflag'
      ENDIF 

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! find edge pairs
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      PRINT "(A)", 'finding edge pairs'

          ned = 0
          edflag(:,:) = 0
          ged2nn_temp(:,:) = 0
          ged2led_temp(:,:) = 0
          ged2el_temp(:,:) = 0
   elem1: DO el1 = 1,ne ! loop through trial elements
   
            nvert1 = mesh%nverts(mesh%el_type(el1))

 local_ed1: DO led1 = 1,nvert1 ! loop through trial edges

              IF(edflag(led1,el1) == 1) THEN ! skip if edge has already been flagged
                CYCLE local_ed1
              ENDIF

              ned = ned + 1 ! increment edge number

              n1ed1 = mesh%ect(mod(led1+0,nvert1)+1,el1) ! find nodes on trial edge
              n2ed1 = mesh%ect(mod(led1+1,nvert1)+1,el1)
              
              ged2nn_temp(1,ned) = n1ed1 ! set nodes on global edge # ned
              ged2nn_temp(2,ned) = n2ed1

              ged2led_temp(1,ned) = led1 ! set local edge number of first element sharing the edge
              ged2el_temp(1,ned) = el1 ! set the first element that shares the edge

              edflag(led1,el1) = 1 ! flag the edge so it is not repeated

       elem2: DO el = 1,mesh%nepn(n1ed1) ! loop through test elements (only those that contain node n1ed1)

                el2 = mesh%epn(el,n1ed1) ! choose a test element that contains node n1ed1 
                
                IF(el2 == el1) THEN ! skip if the test element is the same as the trial element
                  CYCLE elem2
                ENDIF
                
                nvert2 = mesh%nverts(mesh%el_type(el2))

     local_ed2: DO led2 = 1,nvert2 ! loop through local test edge numbers
                  
                  n1ed2 = mesh%ect(MOD(led2+0,nvert2)+1,el2) ! find nodes on test edge
                  n2ed2 = mesh%ect(MOD(led2+1,nvert2)+1,el2)

                  IF(((n1ed1 == n1ed2) .AND. (n2ed1 == n2ed2)) .OR. & ! check if nodes on trial edge matches test edge
                     ((n1ed1 == n2ed2) .AND. (n2ed1 == n1ed2))) THEN

                     ged2led_temp(2,ned) = led2 ! set local edge number of second element sharing the edge
                     ged2el_temp(2,ned) = el2 ! set the second element that shares the edge
                     edflag(led2,el2) = 1 ! flag the edge so it is not repeated

                     EXIT elem2
          
                  ENDIF

                ENDDO local_ed2

              ENDDO elem2

            ENDDO local_ed1
      
          ENDDO elem1

      PRINT "(A,I7)", '   number of total edges:', ned
      PRINT "(A)", ' '

      ALLOCATE(mesh%ged2nn(2,ned),mesh%ged2el(2,ned),mesh%ged2led(2,ned),mesh%bel2bed(ne,5),STAT = alloc_status)
      IF(alloc_status /= 0) THEN
        PRINT*, 'Allocation error: ged2nn,ged2el,ged2led'
      ENDIF 

      mesh%ned = ned
      mesh%ged2nn(:,1:ned) = ged2nn_temp(:,1:ned)
      mesh%ged2el(:,1:ned) = ged2el_temp(:,1:ned)
      mesh%ged2led(:,1:ned) = ged2led_temp(:,1:ned)

      
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! find interior edges
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      PRINT "(A)", 'finding interior edges'
      PRINT*, " "

      ALLOCATE(mesh%bed_flag(ned),mesh%nepe(ne),mesh%el2el(ne,4))            
      mesh%bed_flag(:) = 0
      mesh%bel2bed(:,:) = 0
      mesh%nepe(:) = 0
      
      mesh%nied = 0
      DO ged = 1,ned
        el1 = mesh%ged2el(1,ged)
        el2 = mesh%ged2el(2,ged)
        IF ((el1 /= 0) .AND. (el2 /= 0)) THEN
        
          mesh%nied = mesh%nied + 1
          
          mesh%nepe(el1) = mesh%nepe(el1) + 1
          mesh%nepe(el2) = mesh%nepe(el2) + 1          
          
          mesh%el2el(el1,mesh%nepe(el1)) = el2 
          mesh%el2el(el2,mesh%nepe(el2)) = el1
          
        ELSE
        
          mesh%bed_flag(ged) = 1
          
          IF ( el1 /= 0 ) THEN
            mesh%bel2bed(el1,1) = mesh%bel2bed(el1,1) + 1
            ed = mesh%bel2bed(el1,1)
            mesh%bel2bed(el1,ed+1) = ged
          ELSE IF (el2 /= 0) THEN
            mesh%bel2bed(el1,1) = mesh%bel2bed(el1,1) + 1
            ed = mesh%bel2bed(el1,1)          
            mesh%bel2bed(el2,ed+1) = ged
          ENDIF
          
        ENDIF
      ENDDO  
      
      
      
      
      ALLOCATE(nfbnd_temp(ned))
      mesh%nnfbed = 0
      
      DO seg = 1,mesh%nbou
      
        segtype = mesh%fbseg(2,seg)
              
        DO nd = 1,mesh%fbseg(1,seg)-1
          n1bed = mesh%fbnds(nd,seg)
          n2bed = mesh%fbnds(nd+1,seg)
          found = 0 
  edges2: DO ged = 1,ned
            n1ed2 = mesh%ged2nn(1,ged)
            n2ed2 = mesh%ged2nn(2,ged)
            IF(((n1ed2 == n1bed).AND.(n2ed2 == n2bed)).OR. &
               ((n1ed2 == n2bed).AND.(n2ed2 == n1bed))) THEN

              ! no normal flow edges
              IF( segtype == 0 .OR. segtype == 10 .OR. segtype == 20  .OR. &   ! land boundaries
                  segtype == 1 .OR. segtype == 11 .OR. segtype == 21 ) THEN    ! island boundaries'
                mesh%nnfbed = mesh%nnfbed + 1
                nfbnd_temp(mesh%nnfbed) = ged
                mesh%bed_flag(ged) = 20                
                found = 1
                EXIT edges2               
              ENDIF

              ! specified normal flow edges
              IF ( segtype == 2 .OR. segtype == 12 .OR. segtype == 22 ) THEN
                found = 1
                EXIT edges2
              ENDIF

            ENDIF
          ENDDO edges2
          IF (found == 0) THEN
            PRINT "(A)", "  edge not found"
          ELSE 
!             PRINT "(A,I5)", "  edge found", nd
          ENDIF
        ENDDO
      ENDDO      
      
      ALLOCATE(mesh%nfbedn(mesh%nnfbed))
      
      mesh%nfbedn(1:mesh%nnfbed) = nfbnd_temp(1:mesh%nnfbed)
      
      
      PRINT "(A)", "---------------------------------------------"
      PRINT*, ""



      RETURN
      END SUBROUTINE connect