      MODULE nodal_attributes_mod
      
      USE globals, ONLY: rp
      
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
        INTEGER, DIMENSION(:), ALLOCATABLE :: attr_index        
        REAL(rp), DIMENSION(:), ALLOCATABLE :: default_val
        REAL(rp), DIMENSION(:,:), ALLOCATABLE :: defined_val
      END TYPE


      CONTAINS

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
      CALL read_defined_attributes(fort13%nattr,fort13%name,fort13%nvals,fort13%nndns,fort13%defined_val)
      
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
      
      PRINT*, "reading header"
      OPEN(UNIT=13,FILE=file_name)
      READ(13,"(A100)") grid_name      
      READ(13,*) nn
      READ(13,*) nattr
      
      PRINT*, "reading values"
      n = 0
      ALLOCATE(default_val(30))
      ALLOCATE(name(nattr),units(nattr))
      ALLOCATE(nvals(nattr))
      ALLOCATE(attr_index(nattr))
      DO attr = 1,nattr
      
        READ(13,*) name(attr)%line
        READ(13,*) units(attr)%line
        READ(13,*) nvals(attr)
        
        PRINT*, name(attr)%line, units(attr)%line, nvals(attr)
        
        READ(13,*) (default_val(n+k), k = 1,nvals(attr))
        PRINT*, (default_val(n+k), k = 1,nvals(attr))
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
      REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: defined_val
      
      INTEGER :: i,j
      
      ALLOCATE(defined_val(nn,n))
      DO i = 1,n
        DO j = 1,nn
          defined_val(j,i) = default_val(i)
        ENDDO
      ENDDO
      
      RETURN
      END SUBROUTINE apply_default_attributes
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
      
      SUBROUTINE read_defined_attributes(nattr,name,nvals,nndns,defined_val)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: nattr
      TYPE(char_array), DIMENSION(:), INTENT(IN) :: name
      INTEGER, DIMENSION(:), INTENT(IN) :: nvals
      INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: nndns
      REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: defined_val
      
      INTEGER :: attr
      INTEGER :: i,k,nd
      INTEGER :: n
      CHARACTER(80) :: check
      CHARACTER(80) :: units
      
      n = 0
      ALLOCATE(nndns(nattr))
      DO attr = 1,nattr
      
        READ(13,*) check
        PRINT*, check
        PRINT*, name(attr)%line
        IF (check /= name(attr)%line) THEN
          PRINT*, "Matching error in fort.13"
          STOP
        ENDIF
        READ(13,*) nndns(attr)
        PRINT*, nndns(attr)
        
        DO i = 1,nndns(attr)
          READ(13,*) nd, (defined_val(nd,n+k), k = 1,nvals(attr))
        ENDDO
        n = n + nvals(attr)
      
      ENDDO
               
      CLOSE(13)
   
      RETURN
      END SUBROUTINE read_defined_attributes
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 
      SUBROUTINE find_mode()
      
      IMPLICIT NONE
            
      DO attr = 1,nattr
        n = 0
        DO k = 1,nvals(attr)
          DO nd = 1,nn
            defined_val(nd,n+k)           
          ENDDO
        ENDDO        
        n = n + nvals(attr)
      ENDDO
      
      RETURN
      END SUBROUTINE find_mode

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE write_fort13()
      
      
      RETURN
      END SUBROUTINE write_fort13
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE write_default_attribute()
      
      
      RETURN
      END SUBROUTINE write_default_attribute
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     SUBROUTINE write_defined_attribute()
     
     
     RETURN
     END SUBROUTINE write_defined_attribute

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
     
      END MODULE nodal_attributes_mod