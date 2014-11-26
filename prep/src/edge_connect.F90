      SUBROUTINE edge_connect()

      USE globals, ONLY: pres,nn,ne,ned,vct,el_type,nverts,curved_grid, &
                         mnepn,epn,nepn, &
                         ged2nn,ged2el,ged2led, &
                         nied,iedn,nobed,obedn,nfbed,fbedn,nnfbed,nfbedn,nbed,bedn, &
                         nope,neta,obseg,obnds,nbou,fbseg,nvel,fbnds, &
                         nelnds,nred,redn 
                         
      USE allocation, ONLY: alloc_connect_arrays       
                         

      IMPLICIT NONE
      INTEGER :: el,el1,el2,led1,led2,i,seg,m,ged,ed1,ed2,nd
      INTEGER :: n1,nnds
      INTEGER :: n1ed1,n2ed1,n1ed2,n2ed2
      INTEGER :: n1bed,n2bed
      INTEGER :: segtype
      INTEGER :: alloc_status
      INTEGER :: found
      INTEGER :: nvert1,nvert2,nvert
      REAL(pres) :: x1,x2,x3,y1,y2,y3
      
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ied_temp,bnd_temp,nfbnd_temp,fbnd_temp ! temporary arrary for boundary edges
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ged2nn_temp, ged2el_temp, ged2led_temp ! temporary arrays for edge connectivity
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: edflag
      
      INTEGER, ALLOCATABLE, DIMENSION(:) :: recv_edge      


      ALLOCATE(ied_temp(3*ne),bnd_temp(3*ne),nfbnd_temp(3*ne),fbnd_temp(3*ne),STAT = alloc_status)
      IF(alloc_status /= 0) PRINT*, 'Allocation error: bnd_temp,nfbnd_temp,fbnd_temp'
        
      ALLOCATE(ged2nn_temp(2,3*ne),ged2el_temp(2,3*ne),ged2led_temp(2,3*ne),STAT = alloc_status)
      IF(alloc_status /= 0) PRINT*, 'Allocation error: ged2nn_temp,ged2el_temp,ged2led_temp'       

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! count the number of elements per node
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      nepn(:) = 0
      DO el = 1,ne
        nvert = nverts(el_type(el))
        DO nd = 1,nvert
          n1 = vct(nd,el)
          nepn(n1) = nepn(n1) + 1
        ENDDO
        
      ENDDO

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! find the maximum number of elements per node
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      mnepn = 0
      DO i = 1,nn
        IF(nepn(i) > mnepn) THEN
          mnepn = nepn(i)
        ENDIF
      ENDDO
      
      CALL alloc_connect_arrays(1)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! find the elements associated with each node
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      nepn(:) = 0
      DO el = 1,ne
        nvert = nverts(el_type(el))
        DO nd = 1,nvert
          n1 = vct(nd,el)
          nepn(n1) = nepn(n1) + 1
          epn(nepn(n1),n1) = el
        ENDDO

      ENDDO

      ALLOCATE(edflag(4,ne),STAT = alloc_status)
      IF(alloc_status /= 0) PRINT*, 'Allocation error: edflag'

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! find edge pairs
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          ned = 0
          edflag(:,:) = 0
          ged2nn_temp(:,:) = 0
          ged2led_temp(:,:) = 0
          ged2el_temp(:,:) = 0
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
              ged2el_temp(1,ned) = el1 ! set the first element that shares the edge

              edflag(led1,el1) = 1 ! flag the edge so it is not repeated

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
                     ged2el_temp(2,ned) = el2 ! set the second element that shares the edge
                     edflag(led2,el2) = 1 ! flag the edge so it is not repeated

                     EXIT elem2
          
                  ENDIF

                ENDDO local_ed2

              ENDDO elem2

            ENDDO local_ed1
      
          ENDDO elem1


      CALL alloc_connect_arrays(2)

      ged2nn(:,1:ned) = ged2nn_temp(:,1:ned)
      ged2el(:,1:ned) = ged2el_temp(:,1:ned)
      ged2led(:,1:ned) = ged2led_temp(:,1:ned)
      
      DEALLOCATE(ged2nn_temp,ged2el_temp,ged2led_temp, STAT = alloc_status)
      IF (alloc_status /= 0) PRINT*, "Deallocation error: ged2nn_temp, ged2el_temp, ged2led_temp"

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! find interior edges
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      ALLOCATE(recv_edge(ned))
      recv_edge = 1      

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
        ELSE
          nbed = nbed + 1
          bedn(nbed) = ged
        ENDIF
      ENDDO



      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! find open boundary edges
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      nobed = 0 
      bnd_temp(:) = 0
      DO seg = 1,nope
        DO nd = 1,obseg(seg)-1
          n1bed = obnds(nd,seg)
          n2bed = obnds(nd+1,seg)
  edges1: DO ed2 = 1,nbed
            ged = bedn(ed2)
            n1ed2 = ged2nn(1,ged)
            n2ed2 = ged2nn(2,ged)
            IF(((n1ed2 == n1bed).AND.(n2ed2 == n2bed)).OR. &
               ((n1ed2 == n2bed).AND.(n2ed2 == n1bed))) THEN
              nobed = nobed + 1
              bnd_temp(nobed) = ged
              recv_edge(ged) = 0              
              EXIT edges1
            ENDIF
          ENDDO edges1
        ENDDO
      ENDDO



      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! find flow boundary edges
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      nnfbed = 0 ! # no normal flow edge
      nfbed = 0 ! # flow specifed edges
      fbnd_temp(:) = 0
      nfbnd_temp(:) = 0
      found = 0

      DO seg = 1,nbou
      
        segtype = fbseg(2,seg)
        print*, seg, fbseg(1,seg), segtype
              
        DO nd = 1,fbseg(1,seg)-1
          n1bed = fbnds(nd,seg)
          n2bed = fbnds(nd+1,seg)
          found = 0 
   edges2: DO ed2 = 1,nbed
            ged = bedn(ed2)
            n1ed2 = ged2nn(1,ged)
            n2ed2 = ged2nn(2,ged)
            IF(((n1ed2 == n1bed).AND.(n2ed2 == n2bed)).OR. &
               ((n1ed2 == n2bed).AND.(n2ed2 == n1bed))) THEN

              ! no normal flow edges
              IF( segtype == 0 .OR. segtype == 10 .OR. segtype == 20  .OR. &
                  segtype == 1 .OR. segtype == 11 .OR. segtype == 21 ) THEN
                nnfbed = nnfbed + 1
                nfbnd_temp(nnfbed) = ged
                recv_edge(ged) = 0
                found = 1
                EXIT edges2
              ENDIF

              ! specified normal flow edges
              IF ( segtype == 2 .OR. segtype == 12 .OR. segtype == 22 ) THEN
                nfbed = nfbed + 1
                fbnd_temp(nfbed) = ged
                recv_edge(ged) = 0
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

      CALL alloc_connect_arrays(3)

      iedn(1:nied) = ied_temp(1:nied)      
      obedn(1:nobed) = bnd_temp(1:nobed)      
      nfbedn(1:nnfbed) = nfbnd_temp(1:nnfbed)
      fbedn(1:nfbed) = fbnd_temp(1:nfbed)


      
      nred = 0
          
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
        ENDIF
      ENDDO
     



      

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


! !       write edge connectivity information in similar format to fort.17
!       OPEN(UNIT=17,FILE=dirname(1:lname)//'/'//'fort.17')
!       
!       WRITE(17,*) ned
!       DO i = 1,ned
!         WRITE(17,*) i,ged2nn(1,i),ged2nn(2,i),ged2el(1,i),ged2el(2,i)
!       ENDDO
! 
!       WRITE(17,*) 'number of interior edges:', nied
!       DO i = 1,nied
!         WRITE(17,*) i,iedn(i),ged2nn(1,iedn(i)),ged2nn(2,iedn(i))
!       ENDDO
! 
!       WRITE(17,*) 'number of no normal flow boundary edges:', nnfbed
!       DO i = 1,nnfbed
!         WRITE(17,*) i,nfbedn(i),ged2nn(1,nfbedn(i)),ged2nn(2,nfbedn(i))
!       ENDDO
! 
!       WRITE(17,*) 'number of open boundary edges:', nobed
!       DO i = 1,nobed
!         WRITE(17,*) i,obedn(i),ged2nn(1,obedn(i)),ged2nn(2,obedn(i))
!       ENDDO
! 
!       WRITE(17,*) 'number of flow specified boundary edges:', nfbed
!       DO i = 1,nfbed
!         WRITE(17,*) i,fbedn(i),ged2nn(1,fbedn(i)),ged2nn(2,fbedn(i))
!       ENDDO
! 
!       WRITE(17,*) "global to local edge table"
!       DO i = 1,ned
!         WRITE(17,*) i,ged2led(1,i),ged2led(2,i)
!       ENDDO
!       
!       CLOSE(17)

      RETURN
      END SUBROUTINE edge_connect