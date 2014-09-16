      SUBROUTINE read_grid()
      
      USE globals, ONLY: grid_file,ne,nn,ect,vct,xy,depth,nelnds,elxy,elhb, &
                         nope,neta,obseg,obnds,nvel,nbou,fbseg,fbnds,grid_name, &
                         el_type,ctp,mnelnds,curved_grid,nverts

      IMPLICIT NONE
      INTEGER :: i,j,k,el
      INTEGER :: nbseg
      INTEGER :: btype
      INTEGER :: nvert
      INTEGER :: alloc_status      
      
      nverts(1) = 3
      nverts(2) = 4
      nverts(3) = 3
      nverts(4) = 4      

      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", "             Grid Information                "
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", " "
       

      ! open fort.14 grid file
      OPEN(UNIT = 14, FILE = grid_file)    

      PRINT "(A,A)", "Grid file: ", grid_file                          
                       
      ! read in name of grid
      READ(14,"(A)"), grid_name                                         

      PRINT "(A,A)", "Grid name: ", grid_name

      ! read in number of elements and number of nodes
      READ(14,*), ne, nn                                                

      PRINT "(A,I5)", "Number of elements: ", ne
      PRINT "(A,I5)", "Number of nodes: ", nn
      PRINT*, " "

      ALLOCATE(ect((ctp+1)*(ctp+1),ne),vct(4,ne),xy(2,nn),depth(nn),nelnds(ne),el_type(ne),STAT = alloc_status)  ! element connectivity table, node coordinate array, depth vector      
      IF(alloc_status /= 0) THEN
        PRINT*, 'Allocation error: ect,xy,depth'
      ENDIF     
      ALLOCATE(elxy((ctp+1)*(ctp+1),ne,2),elhb((ctp+1)*(ctp+1),ne))

      curved_grid = 0
      
      ! read in node coordinates and depths
      DO i = 1,nn                                                      
        READ(14,*), j, xy(1,j), xy(2,j), depth(j)
      ENDDO
!       PRINT "(A)", "Node coordinates and depth: "
!       DO i = 1,nn
!         PRINT "(I5,3(F11.3,3x))", i,xy(1,i), xy(2,i), depth(i)
!       ENDDO
!       PRINT*, " "

      ! read in element connectivity
      ect = 0
      DO i = 1,ne
        READ(14,*) el,nelnds(el),(ect(j,el),j = 1,nelnds(el))
        IF (nelnds(el) == 3) THEN
          el_type(el) = 1
        ELSE IF (nelnds(el) == 4) THEN
          el_type(el) = 2
        ELSE IF (nelnds(el) == (ctp+1)*(ctp+2)/2) THEN
          el_type(el) = 3
          curved_grid = 1
        ELSE IF (nelnds(el) == (ctp+1)*(ctp+1)) THEN
          el_type(el) = 4
          curved_grid = 1
        ELSE
          PRINT*, "Element type not supported or ctp not compatible with grid"
          STOP
        ENDIF 
        
        DO j = 1,nelnds(el)
          elxy(j,el,1) = xy(1,ect(j,el))
          elxy(j,el,2) = xy(2,ect(j,el))
          elhb(j,el)   = depth(ect(j,el))
        ENDDO      
        
      ENDDO
      
      PRINT "(A,I5)", "Curved grid = ",curved_grid
      PRINT*, " "
      
      IF (curved_grid == 1) THEN
        DO i = 1,ne
          nvert = nverts(el_type(i))
          DO j = 1,nvert
            vct(j,i) = ect(ctp*(j-1)+1,i)
          ENDDO
        ENDDO        
      ELSE 
        DO i = 1,ne
          nvert = nverts(el_type(i))
          DO j = 1,nvert
            vct(j,i) = ect(j,i)
          ENDDO
        ENDDO
      ENDIF
      
      mnelnds = maxval(nelnds)
      
!       PRINT "(A)", "Element connectivity table: "
!       DO i = 1,ne
!         PRINT "(2(I5,3x),8x,4(I5,3x))", i,nelnds(i),(ect(j,i),j=1,nelnds(i))
!       ENDDO
!       PRINT*, " "

      READ(14,*) nope  ! number of open boundaries                                                 
      READ(14,*) neta  ! number of total elevation specified boundary nodes

      ALLOCATE(obseg(nope),obnds(neta,nope),STAT = alloc_status)  ! number of nodes in each open boundary segment, array for open boundary nodes
      IF(alloc_status /= 0) THEN
        PRINT*, 'Allocation error: obseg,obnds'
      ENDIF

      DO i = 1,nope                                                     
        READ(14,*), nbseg  ! read in # of nodes in segment, boundary type
        obseg(i) = nbseg
        DO j = 1,nbseg
          READ(14,*) obnds(j,i) ! read in open boundary node numbers
        ENDDO
      ENDDO
!       PRINT "(A)", "Open boundary segments:"
!       DO i = 1,nope
!         nbseg = obseg(i)
!         PRINT "(A,I5,A,I5,A)", "Open boundary segment ",i," contains ",nbseg," nodes"
!         DO j = 1,nbseg
!           PRINT "(I5)",obnds(j,i)
!         ENDDO
!       ENDDO
!       PRINT*, " "

      READ(14,*) nbou  ! number of normal flow boundaries
      READ(14,*) nvel  ! total number of normal flow nodes

      ALLOCATE(fbseg(2,nbou),fbnds(nvel,nbou),STAT = alloc_status)  ! array to indicate number of nodes and type of each normal flow boundary, array for open boundary nodes
      IF(alloc_status /= 0) THEN
        PRINT*, 'Allocation error: nfbseg,nfbnds'
      ENDIF

      DO i = 1,nbou
        READ(14,*), nbseg, btype ! read in # of nodes in segment, boundary type
        fbseg(1,i) = nbseg
        fbseg(2,i) = btype
        DO j = 1,nbseg
          READ(14,*), fbnds(j,i)  ! read in normal flow boundary node numbers
        ENDDO
        IF (btype == 1 .OR. btype == 11 .OR. btype == 21) THEN
          IF (fbnds(nbseg,i) /= fbnds(1,i)) THEN
            fbnds(nbseg+1,i) = fbnds(1,i)  ! close island boundaries
            fbseg(1,i) = fbseg(1,i) + 1
          ENDIF
        ENDIF
        
      ENDDO
!       PRINT "(A)", "Normal flow boundary segments: "
!       DO i = 1,nbou
!         nbseg = fbseg(1,i)
!         btype = fbseg(2,i)
!         PRINT "(A,I3,A,I3,A,I5,A)", "Normal flow boundary segment ",i," type ",btype, " contains ",nbseg," nodes"
!         DO j = 1,nbseg
!           PRINT "(I5)", fbnds(j,i)
!         ENDDO
!       ENDDO
!       PRINT*, " "

      CLOSE(14) 

      RETURN
      END SUBROUTINE read_grid