      PROGRAM pointer_test
      IMPLICIT NONE
      INTEGER, PARAMETER :: pres = kind(1d0)

      INTEGER :: m,i,j

      TYPE :: ptr_array
        REAL(pres), POINTER :: ptr
      END TYPE ptr_array

      TYPE(ptr_array), ALLOCATABLE, TARGET, DIMENSION(:) :: test

      REAL(pres), TARGET, DIMENSION(5,2) :: H
      REAL(pres), DIMENSION(10) :: H2

      TYPE(ptr_array), POINTER :: ptr2type

      ALLOCATE(test(10))

      m = 0
      DO i = 1,2
        DO j = 1,5
          m = m+1
          H(j,i) = m
          test(m)%ptr => H(j,i)
        ENDDO
      ENDDO

      DO i = 1,5
        PRINT*, (H(i,j),j=1,2)
      ENDDO

      PRINT*, ' ' 

      DO i = 1,m
        PRINT*, test(i)%ptr
      ENDDO

      PRINT*, ' '


      m = 10
      DO i = 1,2
        DO j = 1,5
          m = m+1
          H(j,i) = m
        ENDDO
      ENDDO

      DO i = 1,5
        PRINT*, (H(i,j),j=1,2)
      ENDDO

      PRINT*, ' ' 

      DO i = 1,10
        PRINT*, test(i)%ptr
      ENDDO

      PRINT*, ' '

      DO i = 1,10
        H2(i) = test(i)%ptr*test(i)%ptr
      ENDDO

      DO i = 1,10
        PRINT*, H2(i)
      ENDDO
        
      ptr2type => test

      DO i = 1,10
        PRINT*, ptr2type(i)
      ENDDO
      

      END PROGRAM pointer_test