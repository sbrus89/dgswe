      MODULE read_grid
      
      USE globals, ONLY: grid,base,eval
      
      CONTAINS
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
      SUBROUTINE read_grids()
      
      USE allocation, ONLY: sizes
      
      IMPLICIT NONE
      
      INTEGER :: i,j
      
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", "           Base Grid Information             "
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", " "      
      
      CALL sizes(base)
      CALL read_nodes(base)      
      
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", "            Eval Grid Information            "
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", " "     
     
      CALL sizes(eval)
      CALL read_nodes(eval)
      
  
      
      RETURN
      END SUBROUTINE read_grids
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE read_nodes(mesh)
      
      USE globals, ONLY: rp,grid
      USE allocation, ONLY: grid_alloc

      IMPLICIT NONE
      INTEGER :: i,j,k,el,n,nd
      INTEGER :: ne,nn,ctp,nvert,nbseg,btype,nope,nbou
      INTEGER :: mnepn,n1,n2,found            
      INTEGER :: curved_grid
      INTEGER, ALLOCATABLE, DIMENSION(:) :: vflag     
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: vxy_temp
      TYPE(grid) :: mesh
  
      ! open fort.14 grid file
      OPEN(UNIT = 14, FILE = mesh%grid_file)      
      
      PRINT "(A,A)", "Grid file: ", mesh%grid_file                          
                       
      ! read in name of grid
      READ(14,"(A)"), mesh%grid_name                                         

      PRINT "(A,A)", "Grid name: ", mesh%grid_name

      ! read in number of elements and number of nodes
      READ(14,*), mesh%ne, mesh%nn                                                

      PRINT "(A,I5)", "Number of elements: ", mesh%ne
      PRINT "(A,I5)", "Number of nodes: ", mesh%nn
      PRINT*, " "
      
      nn = mesh%nn
      ne = mesh%ne
      ctp = mesh%ctp

      CALL grid_alloc(1,mesh)

      curved_grid = 0
      
      ! read in node coordinates and depths
      DO i = 1,nn                                                      
        READ(14,*), j, mesh%xy(1,j), mesh%xy(2,j), mesh%depth(j)
      ENDDO
!       PRINT "(A)", "Node coordinates and depth: "
!       DO i = nn-15,nn !nn
!         PRINT "(I5,3(F11.3,3x))", i,mesh%xy(1,i), mesh%xy(2,i), mesh%depth(i)
!       ENDDO
!       PRINT*, " "

      ! read in element connectivity
      mesh%ect = 0
      DO i = 1,ne
        READ(14,*) el,mesh%nelnds(el),(mesh%ect(j,el),j = 1,mesh%nelnds(el))
        IF (mesh%nelnds(el) == 3) THEN
          mesh%el_type(el) = 1
        ELSE IF (mesh%nelnds(el) == 4) THEN
          mesh%el_type(el) = 2
        ELSE IF (mesh%nelnds(el) == (mesh%ctp+1)*(mesh%ctp+2)/2) THEN
          mesh%el_type(el) = 3
          curved_grid = 1
        ELSE IF (mesh%nelnds(el) == (mesh%ctp+1)*(mesh%ctp+1)) THEN
          mesh%el_type(el) = 4
          curved_grid = 1
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
      
      PRINT "(A,I5)", "Curved grid = ",curved_grid
      PRINT*, " "
      
      
      ! find element verticies
      ALLOCATE(vxy_temp(2,mesh%nn),vflag(mesh%nn))
      vflag = 0
      
      IF(curved_grid == 1) THEN
        DO i = 1,ne
          nvert = mesh%nverts(mesh%el_type(i))
          DO j = 1,nvert
            mesh%vct(j,i) = mesh%ect(ctp*(j-1)+1,i)
          ENDDO
        ENDDO        
      ELSE ! in case grid isn't curved but ctp > 1
        DO i = 1,ne
          nvert = mesh%nverts(mesh%el_type(i))
          DO j = 1,nvert
            mesh%vct(j,i) = mesh%ect(j,i) 
          ENDDO
        ENDDO           
      ENDIF

      
      n = 0
      DO i = 1,ne
        nvert = mesh%nverts(mesh%el_type(i))
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

      DO i = 1,n
        mesh%vxy(1,i) = vxy_temp(1,i)
        mesh%vxy(2,i) = vxy_temp(2,i)
      ENDDO
      
      mesh%mnelnds = maxval(mesh%nelnds)
      
!       PRINT "(A)", "Element connectivity table: "
!       DO i = 1,15
!         PRINT "(2(I5,3x),8x,9(I5,3x))", i,mesh%nelnds(i),(mesh%ect(j,i),j=1,mesh%nelnds(i))
!       ENDDO
!       PRINT*, " "


      ! Find elements associated with each node            
      ALLOCATE(mesh%nepn(mesh%nn))
      
      mesh%nepn(:) = 0
      DO el = 1,mesh%ne
        nvert = mesh%nverts(mesh%el_type(el))
        DO nd = 1,nvert
          n1 = mesh%vct(nd,el)
          mesh%nepn(n1) = mesh%nepn(n1) + 1
        ENDDO
      ENDDO
      
      mnepn = maxval(mesh%nepn)
      
      ALLOCATE(mesh%epn(mnepn,mesh%nn))
      
      mesh%nepn(:) = 0
      DO el = 1,mesh%ne
        nvert = mesh%nverts(mesh%el_type(el))
        DO nd = 1,nvert
          n1 = mesh%vct(nd,el)
          mesh%nepn(n1) = mesh%nepn(n1) + 1
          mesh%epn(mesh%nepn(n1),n1) = el
        ENDDO
      ENDDO 


      ! Read in open boundaries
      READ(14,*) mesh%nope  ! number of open boundaries                                                 
      READ(14,*) mesh%neta  ! number of total elevation specified boundary nodes
      
      nope = mesh%nope
      
      CALL grid_alloc(2,mesh)

      DO i = 1,nope                                                     
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
!       

      ! read in normal flow boundaries
      READ(14,*) mesh%nbou  ! number of normal flow boundaries
      READ(14,*) mesh%nvel  ! total number of normal flow nodes

      nbou = mesh%nbou
       
      CALL grid_alloc(3,mesh)

      DO i = 1,nbou
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
      END SUBROUTINE read_nodes
      
      END MODULE read_grid
