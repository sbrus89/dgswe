      SUBROUTINE read_grid()
      
      USE globals, ONLY: grid_file,ne,nn,ect,xy,depth,nope,neta, &
                         obseg,obnds,nvel,nbou,fbseg,fbnds,grid_name


      IMPLICIT NONE
      INTEGER :: i,j,k,el
      INTEGER :: nbseg
      INTEGER :: btype
      INTEGER :: alloc_status
      LOGICAL :: file_exists
      
      file_exists = .FALSE.
      
!       DO 
!         PRINT*, "Input grid file name"
!         READ(*,"(A)") grid_file
!         PRINT*, " "      
!         
!         INQUIRE(FILE=grid_file,EXIST=file_exists)
!         IF(file_exists) THEN
!           EXIT
!         ELSE
!           PRINT*, "file not found, try again"
!         ENDIF
!       ENDDO      

      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", "             Grid Information                "
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", " "
       
      ! open fort.14 grid file       
      INQUIRE(FILE=grid_file,EXIST=file_exists)
      IF(file_exists) THEN
        OPEN(UNIT = 14, FILE = grid_file)     
      ELSE
        PRINT*, "file not found"
        STOP
      ENDIF       
       
 
      
      PRINT "(A,A)", "Grid file: ", grid_file                          
                       
      ! read in name of grid
      READ(14,"(A)"), grid_name                                         

      PRINT "(A,A)", "Grid name: ", grid_name

      ! read in number of elements and number of nodes
      READ(14,*), ne, nn                                                

      PRINT "(A,I8)", "Number of elements: ", ne
      PRINT "(A,I8)", "Number of nodes: ", nn
      PRINT*, " "

      ALLOCATE(ect(3,ne),xy(2,nn),depth(nn),STAT = alloc_status)  ! element connectivity table, node coordinate array, depth vector      
      IF(alloc_status /= 0) THEN
        PRINT*, 'Allocation error: ect,xy,depth'
      ENDIF        

      ! read in node coordinates and depths
      DO i = 1,nn                                                      
        READ(14,*), j, xy(1,j), xy(2,j), depth(j)
      ENDDO
!       PRINT "(A)", "Node coordinates and depth: "
!       DO i = 1,nn
!         PRINT "(3(F11.3,3x))", xy(1,i), xy(2,i), depth(i)
!       ENDDO
!       PRINT*, " "

      ! read in element connectivity
      DO i = 1,ne
        READ(14,*) el,k,ect(1,el),ect(2,el),ect(3,el)                               
      ENDDO
!       PRINT "(A)", "Element connectivity table: "
!       DO i = 1,ne
!         PRINT "(4(I5,3x))", i,ect(1,i),ect(2,i),ect(3,i) 
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
      DO i = 1,nope
        nbseg = obseg(i)
        PRINT "(A,I5,A,I5,A)", "Open boundary segment ",i," contains ",nbseg," nodes"
!         DO j = 1,nbseg
!           PRINT "(I5)",obnds(j,i)
!         ENDDO
      ENDDO
      PRINT*, " "

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
      ENDDO
!       PRINT "(A)", "Normal flow boundary segments: "
      DO i = 1,nbou
        nbseg = fbseg(1,i)
        btype = fbseg(2,i)
        PRINT "(A,I3,A,I3,A,I5,A)", "Normal flow boundary segment ",i," type ",btype, " contains ",nbseg," nodes"
!         DO j = 1,nbseg
!           PRINT "(I5)", fbnds(j,i)
!         ENDDO
      ENDDO
      PRINT*, " "

      CLOSE(14) 

      RETURN
      END SUBROUTINE read_grid
