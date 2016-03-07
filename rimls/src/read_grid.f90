      SUBROUTINE read_grid(mesh)
      
      USE globals, ONLY: rp,grid,ctp,nverts, &
                         Erad,lambda0,phi0
      USE allocation, ONLY: grid_alloc                         

      IMPLICIT NONE
      INTEGER :: i,j,k,el,n,nd
      INTEGER :: ne,nn,nbseg
      INTEGER :: btype
      INTEGER :: nvert
      INTEGER :: alloc_status   
      INTEGER, ALLOCATABLE, DIMENSION(:) :: vflag     
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: vxy_temp 
      TYPE(grid) :: mesh
      
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", "             Grid Information                "
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", " "
       

      ! open fort.14 grid file
      OPEN(UNIT = 14, FILE = mesh%grid_file)    

      PRINT "(A,A)", "Grid file: ", mesh%grid_file                          
                       
      ! read in name of grid
      READ(14,"(A)"), mesh%grid_name                                         

      PRINT "(A,A)", "Grid name: ", mesh%grid_name

      ! read in number of elements and number of nodes
      READ(14,*), mesh%ne, mesh%nn    
      
      ne = mesh%ne
      nn = mesh%nn

      PRINT "(A,I9)", "Number of elements: ", ne
      PRINT "(A,I9)", "Number of nodes: ", nn
      PRINT*, " "

      CALL grid_alloc(1,mesh)

      mesh%curved_grid = 0
      
      ! read in node coordinates and depths
      DO i = 1,nn                                                      
        READ(14,*), j, mesh%xy(1,j), mesh%xy(2,j), mesh%depth(j)
      ENDDO
!       PRINT "(A)", "Node coordinates and depth: "
!       DO i = 1,nn
!         PRINT "(I5,3(F11.3,3x))", i,xy(1,i), xy(2,i), depth(i)
!       ENDDO
!       PRINT*, " "

      DO i = 1,nn
        mesh%xy(1,i) = Erad*(mesh%xy(1,i)-lambda0)*cos(phi0)
        mesh%xy(2,i) = Erad*mesh%xy(2,i)
      ENDDO

      ! read in element connectivity
      mesh%ect = 0
      DO i = 1,ne
        READ(14,*) el,mesh%nelnds(el),(mesh%ect(j,el),j = 1,mesh%nelnds(el))
        IF (mesh%nelnds(el) == 3) THEN
          mesh%el_type(el) = 1
        ELSE IF (mesh%nelnds(el) == 4) THEN
          mesh%el_type(el) = 2
        ELSE IF (mesh%nelnds(el) == (ctp+1)*(ctp+2)/2) THEN
          mesh%el_type(el) = 3
          mesh%curved_grid = 1
        ELSE IF (mesh%nelnds(el) == (ctp+1)*(ctp+1)) THEN
          mesh%el_type(el) = 4
          mesh%curved_grid = 1
        ELSE
          PRINT*, "Element type not supported or ctp not compatible with grid"
          STOP
        ENDIF 
        
        DO j = 1,mesh%nelnds(el)
          mesh%elxy(j,el,1) = mesh%xy(1,mesh%ect(j,el))
          mesh%elxy(j,el,2) = mesh%xy(2,mesh%ect(j,el))
          mesh%elhb(j,el)   = mesh%depth(mesh%ect(j,el))
        ENDDO      
        
      ENDDO
      
      PRINT "(A,I5)", "Curved grid = ", mesh%curved_grid
      PRINT*, " "
      
      ALLOCATE(vxy_temp(2,nn),vflag(nn))
      vflag = 0
      
      IF (mesh%curved_grid == 1) THEN
        DO i = 1,ne
          nvert = nverts(mesh%el_type(i))
          DO j = 1,nvert
            mesh%vct(j,i) = mesh%ect(ctp*(j-1)+1,i)
          ENDDO
        ENDDO        
      ELSE 
        DO i = 1,ne
          nvert = nverts(mesh%el_type(i))
          DO j = 1,nvert
            mesh%vct(j,i) = mesh%ect(j,i)
          ENDDO
        ENDDO
      ENDIF
      
      n = 0
      DO i = 1,ne
        nvert = nverts(mesh%el_type(i))
        DO j = 1,nvert
          nd = mesh%vct(j,i)
          IF (vflag(nd) == 0) THEN  
            n = n + 1            
            mesh%vxyn(n) = nd
            vxy_temp(1,n) = mesh%xy(1,nd)
            vxy_temp(2,n) = mesh%xy(2,nd)
            vflag(nd) = 1
          ENDIF
        ENDDO
      ENDDO         
      
      ALLOCATE(mesh%vxy(2,n))
      
      DO i = 1,n
        mesh%vxy(1,i) = vxy_temp(1,i)
        mesh%vxy(2,i) = vxy_temp(2,i)
      ENDDO
      
      mesh%mnelnds = maxval(mesh%nelnds)      
      
      DO i = 1,nn  
        mesh%xyhv(1,i,1) = mesh%xy(1,i)
        mesh%xyhv(1,i,2) = mesh%xy(2,i)
        mesh%xyhv(1,i,3) = mesh%depth(i)
      ENDDO
      
!       PRINT "(A)", "Element connectivity table: "
!       DO i = 1,ne
!         PRINT "(2(I5,3x),8x,4(I5,3x))", i,nelnds(i),(ect(j,i),j=1,nelnds(i))
!       ENDDO
!       PRINT*, " "

      READ(14,*) mesh%nope  ! number of open boundaries                                                 
      READ(14,*) mesh%neta  ! number of total elevation specified boundary nodes

      CALL grid_alloc(2,mesh)

      DO i = 1,mesh%nope                                                     
        READ(14,*), nbseg  ! read in # of nodes in segment, boundary type
        mesh%obseg(i) = nbseg
        DO j = 1,nbseg
          READ(14,*) mesh%obnds(j,i) ! read in open boundary node numbers
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

      READ(14,*) mesh%nbou  ! number of normal flow boundaries
      READ(14,*) mesh%nvel  ! total number of normal flow nodes

      CALL grid_alloc(3,mesh)
      
      DO i = 1,mesh%nbou
        READ(14,*), nbseg, btype ! read in # of nodes in segment, boundary type
        mesh%fbseg(1,i) = nbseg
        mesh%fbseg(2,i) = btype
        DO j = 1,nbseg
          READ(14,*), mesh%fbnds(j,i)  ! read in normal flow boundary node numbers
        ENDDO
        IF (btype == 1 .OR. btype == 11 .OR. btype == 21) THEN
          IF (mesh%fbnds(nbseg,i) /= mesh%fbnds(1,i)) THEN
            mesh%fbnds(nbseg+1,i) = mesh%fbnds(1,i)  ! close island boundaries
            mesh%fbseg(1,i) = mesh%fbseg(1,i) + 1
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