      ! This module is used to maintain compatibility with the read_dginp module in 
      ! read_input.F90 used in dgswe.  This allows both the dgswe and dgprep to 
      ! use the same input reading subroutines so that features can be added more
      ! easily 
      
      
      MODULE messenger2

      CHARACTER(6), PARAMETER :: dirname = "."
      INTEGER, PARAMETER :: lname = 1
      INTEGER, PARAMETER :: myrank = 0
      INTEGER, PARAMETER :: nthreads = 1


      CONTAINS



      SUBROUTINE abort()

      STOP

      END SUBROUTINE abort



      END MODULE messenger2