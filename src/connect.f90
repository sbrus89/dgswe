      SUBROUTINE connect()

      USE globals, ONLY: pres,nn,ne,ned,ect, &
                         ged2nn,ged2el,ged2led, &
                         nied,iedn,nobed,obedn,nfbed,fbedn,nnfbed,nfbedn, &
                         nope,neta,obseg,obnds,nbou,fbseg,nvel,fbnds
                         

      IMPLICIT NONE
      INTEGER :: el,el1,el2,led1,led2,i,seg,m,ged,ed1,ed2,nd
      INTEGER :: n1,n2,n3
      INTEGER :: n1ed1,n2ed1,n1ed2,n2ed2
      INTEGER :: n1bed,n2bed
      INTEGER :: led(2,3) ! local edge nodes
      INTEGER :: segtype
      INTEGER :: alloc_status
      REAL(pres) :: x1,x2,x3,y1,y2,y3
      
      INTEGER, ALLOCATABLE, DIMENSION(:) :: bnd_temp,nfbnd_temp,fbnd_temp ! temporary arrary for boundary edges
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ged2nn_temp, ged2el_temp, ged2led_temp ! temporary arrays for edge connectivity
      INTEGER, ALLOCATABLE, DIMENSION(:) :: nepn ! number of elements per node
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: epn ! elements per node
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: edflag
      INTEGER :: mnepn ! maximum number of elements per node

      ALLOCATE(bnd_temp(3*ne),nfbnd_temp(3*ne),fbnd_temp(3*ne),STAT = alloc_status)
      IF(alloc_status /= 0) THEN
        PRINT*, 'Allocation error: bnd_temp,nfbnd_temp,fbnd_temp'
      ENDIF 
      ALLOCATE(ged2nn_temp(2,3*ne),ged2el_temp(2,3*ne),ged2led_temp(2,3*ne),STAT = alloc_status)
      IF(alloc_status /= 0) THEN
        PRINT*, 'Allocation error: ged2nn_temp,ged2el_temp,ged2led_temp'
      ENDIF 

      ALLOCATE(nepn(nn),STAT = alloc_status)
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

      nepn(:) = 0
      DO el = 1,ne
        n1 = ect(1,el)
        n2 = ect(2,el)
        n3 = ect(3,el)

        nepn(n1) = nepn(n1) + 1 
        nepn(n2) = nepn(n2) + 1
        nepn(n3) = nepn(n3) + 1
      ENDDO

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! find the maximum number of elements per node
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      PRINT "(A)", 'finding maximum elements per node'

      mnepn = 0
      DO i = 1,nn
        IF(nepn(i) > mnepn) THEN
          mnepn = nepn(i)
        ENDIF
      ENDDO

      PRINT "(A,I7)", '   maximum elements per node:', mnepn
      PRINT "(A)", ' '

      ALLOCATE(epn(mnepn,nn),STAT = alloc_status)
      IF(alloc_status /= 0) THEN
        PRINT*, 'Allocation error: epn'
      ENDIF 

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! find the elements associated with each node
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      PRINT "(A)", 'finding elements associated with each node'
      PRINT "(A)", ' '

      nepn(:) = 0
      DO el = 1,ne
        n1 = ect(1,el)
        n2 = ect(2,el)
        n3 = ect(3,el)

        nepn(n1) = nepn(n1) + 1 
        epn(nepn(n1),n1) = el

        nepn(n2) = nepn(n2) + 1
        epn(nepn(n2),n2) = el

        nepn(n3) = nepn(n3) + 1
        epn(nepn(n3),n3) = el
      ENDDO

      ALLOCATE(edflag(3,ne),STAT = alloc_status)
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

            n1 = ect(1,el1) ! element nodes
            n2 = ect(2,el1)
            n3 = ect(3,el1)

            led(1,1) = n2 ! nodes on local edge 1
            led(2,1) = n3
            led(1,2) = n3 ! nodes on local edge 2
            led(2,2) = n1
            led(1,3) = n1 ! nodes on local edge 3
            led(2,3) = n2

 local_ed1: DO led1 = 1,3 ! loop through trial edges

              IF(edflag(led1,el1) == 1) THEN ! skip if edge has already been flagged
                CYCLE local_ed1
              ENDIF

              ned = ned + 1 ! increment edge number

              n1ed1 = led(1,led1) ! find nodes on trial edge
              n2ed1 = led(2,led1)
              
              ged2nn_temp(1,ned) = n1ed1 ! set nodes on global edge # ned
              ged2nn_temp(2,ned) = n2ed1

              ged2led_temp(1,ned) = led1 ! set local edge number of first element sharing the edge
              ged2el_temp(1,ned) = el1 ! set the first element that shares the edge

              edflag(led1,el1) = 1 ! flag the edge so it is not repeated

       elem2: DO el = 1,nepn(n1ed1) ! loop through test elements (only those that contain node n1ed1)

                el2 = epn(el,n1ed1) ! choose a test element that contains node n1ed1 
                
                IF(el2 == el1) THEN ! skip if the test element is the same as the trial element
                  CYCLE elem2
                ENDIF

     local_ed2: DO led2 = 1,3 ! loop through local test edge numbers
                  
                  n1ed2 = ect(MOD(led2+0,3)+1,el2) ! find nodes on test edge
                  n2ed2 = ect(MOD(led2+1,3)+1,el2)

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

      ALLOCATE(ged2nn(2,ned),ged2el(2,ned),ged2led(2,ned),STAT = alloc_status)
      IF(alloc_status /= 0) THEN
        PRINT*, 'Allocation error: ged2nn,ged2el,ged2led'
      ENDIF 

      ged2nn(:,1:ned) = ged2nn_temp(:,1:ned)
      ged2el(:,1:ned) = ged2el_temp(:,1:ned)
      ged2led(:,1:ned) = ged2led_temp(:,1:ned)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! find interior edges
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      PRINT "(A)", 'finding interior edges'

      nied = 0
      bnd_temp(:) = 0
      DO ged = 1,ned
        el1 = ged2el(1,ged)
        el2 = ged2el(2,ged)
        IF ((el1 /= 0) .AND. (el2 /= 0)) THEN
          nied = nied + 1
          bnd_temp(nied) = ged
        ENDIF
      ENDDO

      ALLOCATE(iedn(nied),STAT = alloc_status)
      IF(alloc_status /= 0) THEN
        PRINT*, 'Allocation error: iedn'
      ENDIF 

      iedn(1:nied) = bnd_temp(1:nied)

      PRINT "(A,I7)", '   number of interior edges:', nied
      PRINT "(A)", ' '

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! find open boundary edges
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      PRINT "(A)", 'finding open boundary edges'

      nobed = 0 
      bnd_temp(:) = 0
      DO seg = 1,nope
        DO nd = 1,obseg(seg)-1
          n1bed = obnds(nd,seg)
          n2bed = obnds(nd+1,seg)
          DO ed2 = 1,ned
            n1ed2 = ged2nn(1,ed2)
            n2ed2 = ged2nn(2,ed2)
            IF(((n1ed2 == n1bed).AND.(n2ed2 == n2bed)).OR. &
               ((n1ed2 == n2bed).AND.(n2ed2 == n1bed))) THEN
              nobed = nobed + 1
              bnd_temp(nobed) = ed2
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      ALLOCATE(obedn(nobed),STAT = alloc_status)
      IF(alloc_status /= 0) THEN
        PRINT*, 'Allocation error: obedn'
      ENDIF 

      obedn(1:nobed) = bnd_temp(1:nobed)

      PRINT "(A,I7)", '   number of open boundary edges:', nobed
      PRINT "(A)", ' '

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! find flow boundary edges
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      PRINT "(A)", 'finding flow boundary edges'

      nnfbed = 0 ! # no normal flow edge
      nfbed = 0 ! # flow specifed edges
      fbnd_temp(:) = 0
      nfbnd_temp(:) = 0

      DO seg = 1,nbou
        DO nd = 1,fbseg(1,seg)-1
          n1bed = fbnds(nd,seg)
          n2bed = fbnds(nd+1,seg)
          DO ed2 = 1,ned
            n1ed2 = ged2nn(1,ed2)
            n2ed2 = ged2nn(2,ed2)
            IF(((n1ed2 == n1bed).AND.(n2ed2 == n2bed)).OR. &
               ((n1ed2 == n2bed).AND.(n2ed2 == n1bed))) THEN

              segtype = fbseg(2,seg)

              ! no normal flow edges
              IF( segtype == 0 .OR. segtype == 10 .OR. segtype == 20 ) THEN
                nnfbed = nnfbed + 1
                nfbnd_temp(nnfbed) = ed2
              ENDIF

              ! specified normal flow edges
              IF ( segtype == 2 .OR. segtype == 12 .OR. segtype == 22 ) THEN
                nfbed = nfbed + 1
                fbnd_temp(nfbed) = ed2
              ENDIF

            ENDIF
          ENDDO
        ENDDO
      ENDDO

      ALLOCATE(fbedn(nfbed),nfbedn(nnfbed),STAT=alloc_status)
      IF(alloc_status /= 0) THEN
        PRINT*, 'Allocation error: fbedn'
      ENDIF    

      nfbedn(1:nnfbed) = nfbnd_temp(1:nnfbed)
      fbedn(1:nfbed) = fbnd_temp(1:nfbed)

      PRINT "(A,I7)", '   number of specified normal boundary edges:', nfbed
      PRINT "(A,I7)", '   number of no normal flow boundary edges:', nnfbed
      PRINT "(A)", ' '

      PRINT "(A,I7)", 'number of missing edges:',ned-(nied+nobed+nfbed+nnfbed)
      PRINT "(A)", ' '



      ! write edge connectivity information in similar format to fort.17
      WRITE(17,*) ned
      DO i = 1,ned
        WRITE(17,*) i,ged2nn(1,i),ged2nn(2,i),ged2el(1,i),ged2el(2,i)
      ENDDO

      WRITE(17,*) 'number of interior edges:', nied
      DO i = 1,nied
        WRITE(17,*) i,iedn(i),ged2nn(1,iedn(i)),ged2nn(2,iedn(i))
      ENDDO

      WRITE(17,*) 'number of no normal flow boundary edges:', nnfbed
      DO i = 1,nnfbed
        WRITE(17,*) i,nfbedn(i),ged2nn(1,nfbedn(i)),ged2nn(2,nfbedn(i))
      ENDDO

      WRITE(17,*) 'number of open boundary edges:', nobed
      DO i = 1,nobed
        WRITE(17,*) i,obedn(i),ged2nn(1,obedn(i)),ged2nn(2,obedn(i))
      ENDDO

      WRITE(17,*) 'number of flow specified boundary edges:', nfbed
      DO i = 1,nfbed
        WRITE(17,*) i,fbedn(i),ged2nn(1,fbedn(i)),ged2nn(2,fbedn(i))
      ENDDO

      WRITE(17,*) "global to local edge table"
      DO i = 1,ned
        WRITE(17,*) i,ged2led(1,i),ged2led(2,i)
      ENDDO

      RETURN
      END SUBROUTINE connect