      MODULE nodal_attributes_mod
      
      USE globals, ONLY: rp
      USE kdtree2_module      
      
      IMPLICIT NONE
      
      TYPE :: char_array
        CHARACTER(80) :: line
      END TYPE      
      
      TYPE :: type13
        CHARACTER(:), ALLOCATABLE :: file_name
        CHARACTER(100) :: grid_name
        INTEGER :: nn
        INTEGER :: n
        INTEGER :: nattr
        
        TYPE(char_array), DIMENSION(:), ALLOCATABLE :: name
        TYPE(char_array), DIMENSION(:), ALLOCATABLE :: units
        
        INTEGER, DIMENSION(:), ALLOCATABLE :: nvals
        INTEGER, DIMENSION(:), ALLOCATABLE :: nndns
        INTEGER, DIMENSION(:,:), ALLOCATABLE :: ndns
        INTEGER, DIMENSION(:), ALLOCATABLE :: attr_index        
        REAL(rp), DIMENSION(:), ALLOCATABLE :: default_val
        REAL(rp), DIMENSION(:,:), ALLOCATABLE :: defined_val
      END TYPE


      CONTAINS
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE print_fort13_info(fort13)
      
      IMPLICIT NONE
      
      TYPE(type13), INTENT(IN) :: fort13      
      
      INTEGER :: attr
      INTEGER :: n
      INTEGER :: k
      
      n = 0
      DO attr = 1,fort13%nattr
        PRINT*, TRIM(fort13%name(attr)%line)
        PRINT("(12(F15.5))"), (fort13%default_val(n+k), k = 1,fort13%nvals(attr))
        PRINT*, fort13%nndns(attr)
        PRINT("(12(F15.5))"), (MAXVAL(fort13%defined_val(:,n+k)), k = 1,fort13%nvals(attr))
        PRINT("(12(F15.5))"), (MINVAL(fort13%defined_val(:,n+k)), k = 1,fort13%nvals(attr))        
        n = n + fort13%nvals(attr)
        PRINT*, "" 
      ENDDO      
      
      RETURN
      END SUBROUTINE print_fort13_info

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE read_fort13(fort13)
      
      IMPLICIT NONE
      
      TYPE(type13), INTENT(INOUT) :: fort13
      
      PRINT*, "Reading default values"
      CALL read_default_attribute(fort13%file_name,fort13%grid_name,fort13%nn,fort13%nattr, &
                                  fort13%n,fort13%name,fort13%units,fort13%nvals,fort13%attr_index,fort13%default_val)
      PRINT*, "Applying nodal attributes"
      CALL apply_default_attributes(fort13%nn,fort13%n,fort13%default_val,fort13%defined_val)
      PRINT*, "Reading defined attributes"      
      CALL read_defined_attributes(fort13%nattr,fort13%name,fort13%nvals,fort13%nn,fort13%nndns,fort13%ndns,fort13%defined_val)
      
      RETURN
      END SUBROUTINE read_fort13

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE read_default_attribute(file_name,grid_name,nn,nattr,n,name,units,nvals,attr_index,default_val)
      
      IMPLICIT NONE
      
      CHARACTER(*), INTENT(IN) :: file_name
      CHARACTER(*), INTENT(OUT) :: grid_name
      INTEGER, INTENT(OUT) :: nn
      INTEGER, INTENT(OUT) :: nattr
      INTEGER, INTENT(OUT) :: n
      TYPE(char_array), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: name
      TYPE(char_array), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: units
      INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: nvals
      INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: attr_index      
      REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: default_val      
      
      INTEGER :: attr
      INTEGER :: k
      LOGICAL :: file_exists
      
      INQUIRE(FILE=file_name, EXIST=file_exists)
      IF(file_exists .eqv. .FALSE.) THEN
        PRINT*, "fort.13 file does not exist"
        STOP       
      ENDIF      
      
      OPEN(UNIT=13,FILE=file_name)
      READ(13,"(A100)") grid_name      
      READ(13,*) nn
      READ(13,*) nattr

      n = 0
      ALLOCATE(default_val(30))
      ALLOCATE(name(nattr),units(nattr))
      ALLOCATE(nvals(nattr))
      ALLOCATE(attr_index(nattr))
      DO attr = 1,nattr
      
        READ(13,*) name(attr)%line
        READ(13,*) units(attr)%line
        READ(13,*) nvals(attr)
        
        READ(13,*) (default_val(n+k), k = 1,nvals(attr))
        attr_index(attr) = n+1
        n = n + nvals(attr)
              
      ENDDO
      
      
      RETURN
      END SUBROUTINE read_default_attribute
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE apply_default_attributes(nn,n,default_val,defined_val)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: nn
      INTEGER, INTENT(IN) :: n
      REAL(rp), DIMENSION(:), INTENT(IN) :: default_val
      REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: defined_val
      
      INTEGER :: i,j
      
      IF (NOT(ALLOCATED(defined_val))) THEN      
        ALLOCATE(defined_val(nn,n))
      ENDIF
      
      DO i = 1,n
        DO j = 1,nn
          defined_val(j,i) = default_val(i)
        ENDDO
      ENDDO
      
      RETURN
      END SUBROUTINE apply_default_attributes
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
      
      SUBROUTINE read_defined_attributes(nattr,name,nvals,nn,nndns,ndns,defined_val)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: nattr
      TYPE(char_array), DIMENSION(:), INTENT(IN) :: name
      INTEGER, DIMENSION(:), INTENT(IN) :: nvals
      INTEGER, INTENT(IN) :: nn
      INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: nndns
      INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: ndns      
      REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: defined_val
      
      INTEGER :: attr
      INTEGER :: i,k,nd
      INTEGER :: n
      CHARACTER(80) :: check
      CHARACTER(80) :: units
      
      n = 0
      ALLOCATE(nndns(nattr))
      ALLOCATE(ndns(nn,nattr))
      DO attr = 1,nattr
      
        READ(13,*) check
        IF (check /= name(attr)%line) THEN
          PRINT*, "Matching error in fort.13"
          STOP
        ENDIF
        READ(13,*) nndns(attr)
        
        DO i = 1,nndns(attr)
          READ(13,*) nd, (defined_val(nd,n+k), k = 1,nvals(attr))
          ndns(i,attr) = nd
        ENDDO
        n = n + nvals(attr)
      
      ENDDO
               
      CLOSE(13)
   
      RETURN
      END SUBROUTINE read_defined_attributes
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 
!       SUBROUTINE find_mode()
!       
!       IMPLICIT NONE
!             
!       DO attr = 1,nattr
!         n = 0
!         DO k = 1,nvals(attr)
!           DO nd = 1,nn
!             defined_val(nd,n+k)           
!           ENDDO
!         ENDDO        
!         n = n + nvals(attr)
!       ENDDO
!       
!       RETURN
!       END SUBROUTINE find_mode

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE write_fort13(fort13)
      
      TYPE(type13), INTENT(IN) :: fort13
      
      CALL write_default_attribute(fort13%file_name,fort13%grid_name,fort13%nn,fort13%nattr, &
                                   fort13%name,fort13%units,fort13%nvals,fort13%default_val)
      CALL write_defined_attribute(fort13%nattr,fort13%name,fort13%nvals,fort13%nndns,fort13%ndns,fort13%defined_val)
      
      RETURN
      END SUBROUTINE write_fort13
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE write_default_attribute(file_name,grid_name,nn,nattr,name,units,nvals,default_val)
      
      IMPLICIT NONE
      
      CHARACTER(*), INTENT(IN) :: file_name
      CHARACTER(*), INTENT(IN) :: grid_name
      INTEGER, INTENT(IN) :: nn
      INTEGER, INTENT(IN) :: nattr
      TYPE(char_array), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: name
      TYPE(char_array), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: units
      INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(IN) :: nvals    
      REAL(rp), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: default_val      
      
      INTEGER :: attr
      INTEGER :: k
      INTEGER :: n
      
      PRINT*, "writing header"
      OPEN(UNIT=13,FILE=file_name)
      WRITE(13,"(A80)") grid_name      
      WRITE(13,*) nn
      WRITE(13,*) nattr
      
      PRINT*, "writing values"
      n = 0
      DO attr = 1,nattr
      
        WRITE(13,*) TRIM(name(attr)%line)
        WRITE(13,*) TRIM(units(attr)%line)
        WRITE(13,*) nvals(attr)
        
        WRITE(13,*) (default_val(n+k), k = 1,nvals(attr))
        n = n + nvals(attr)
              
      ENDDO      
      
      RETURN
      END SUBROUTINE write_default_attribute
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE write_defined_attribute(nattr,name,nvals,nndns,ndns,defined_val)
     
      IMPLICIT NONE
     
      INTEGER, INTENT(IN) :: nattr
      TYPE(char_array), DIMENSION(:), INTENT(IN) :: name
      INTEGER, DIMENSION(:), INTENT(IN) :: nvals
      INTEGER, DIMENSION(:), INTENT(IN) :: nndns
      INTEGER, DIMENSION(:,:), INTENT(IN) :: ndns 
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: defined_val
      
      INTEGER :: attr
      INTEGER :: i,k,nd
      INTEGER :: n
      
      n = 0
      DO attr = 1,nattr
      
        WRITE(13,*) TRIM(name(attr)%line)
        WRITE(13,*) nndns(attr)
        
        DO i = 1,nndns(attr)
          nd = ndns(i,attr)
          WRITE(13,*) nd, (defined_val(nd,n+k), k = 1,nvals(attr))
        ENDDO
        n = n + nvals(attr)
      
      ENDDO
               
      CLOSE(13)     
     
     
     RETURN
     END SUBROUTINE write_defined_attribute

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     

      SUBROUTINE determine_non_default_nodes(nn,val,default_val,nndns,ndns,defined_val)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: nn
      REAL(rp), DIMENSION(:), INTENT(IN) :: val
      REAL(rp), INTENT(IN) :: default_val      
      INTEGER, INTENT(OUT) :: nndns
      INTEGER, DIMENSION(:), INTENT(OUT) :: ndns
      REAL(rp), DIMENSION(:), INTENT(OUT) :: defined_val
      
      
      INTEGER :: nd
      
      
      nndns = 0
      DO nd = 1,nn
        defined_val(nd) = val(nd)      
        
        IF (abs(val(nd)-default_val) > 1d-12) THEN          
          nndns = nndns + 1
          ndns(nndns) = nd          
        ENDIF
      ENDDO
      
      RETURN
      END SUBROUTINE determine_non_default_nodes


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE create_neighbor_table(nn_search,xy_search,nn_find,xy_find,neighbor_table)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: nn_search
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: xy_search
      INTEGER, INTENT(IN) :: nn_find
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: xy_find
      INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: neighbor_table
      
      TYPE(kdtree2), POINTER :: tree_xy      
      TYPE(kdtree2_result), ALLOCATABLE, DIMENSION(:) :: closest   
      REAL(rp) :: xy(2)    
      INTEGER :: nd
      
      PRINT*, "Building kd-tree"
      tree_xy => kdtree2_create(xy_search(1:2,1:nn_search), rearrange=.true., sort=.true.)
      ALLOCATE(closest(1))          
      PRINT*, "done"
     
      ALLOCATE(neighbor_table(nn_find))
      DO nd = 1,nn_find
      
        xy(1) = xy_find(1,nd)
        xy(2) = xy_find(2,nd)
        
        CALL kdtree2_n_nearest(tp=tree_xy,qv=xy,nn=1,results=closest) 
        
        neighbor_table(nd) = closest(1)%idx

      ENDDO      
      
      RETURN
      END SUBROUTINE
      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE nearest_neighbor_interpolation(nn,neighbor_table,val,default_val,nndns,ndns,defined_val)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: nn
      INTEGER, DIMENSION(:), INTENT(IN) :: neighbor_table
      REAL(rp), DIMENSION(:), INTENT(IN) :: val
      REAL(rp), INTENT(IN) :: default_val
      INTEGER, INTENT(OUT) :: nndns
      INTEGER, DIMENSION(:), INTENT(OUT) :: ndns
      REAL(rp), DIMENSION(:), INTENT(OUT) :: defined_val
      
      INTEGER :: nd
      INTEGER :: clnd
      
      nndns = 0
      DO nd = 1,nn
      
        clnd = neighbor_table(nd)
        defined_val(nd) = val(clnd)        
        
        IF (abs(val(clnd)-default_val) > 1d-12) THEN
          nndns = nndns + 1
          ndns(nndns) = nd
        ENDIF
          
      ENDDO      
      
      
      
      RETURN
      END SUBROUTINE nearest_neighbor_interpolation


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
     
      END MODULE nodal_attributes_mod