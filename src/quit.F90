      MODULE quit 

#ifdef CMPI        
      USE mpi
#endif      

      USE globals, ONLY: rp
      
      IMPLICIT NONE

      CONTAINS 
      


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       SUBROUTINE finish()          

       IMPLICIT NONE
       
       INTEGER :: ierr       
       
#ifdef CMPI       
       CALL MPI_FINALIZE(ierr)
       
       PRINT*, "MPI terminated, status = ", ierr
#endif
       
       STOP

       RETURN
       END SUBROUTINE finish

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       SUBROUTINE abort()

       IMPLICIT NONE
       
       INTEGER :: errorcode
       INTEGER :: ierr
       
       errorcode = 0
        
#ifdef CMPI       
       CALL MPI_ABORT(MPI_COMM_WORLD,errorcode,ierr)
       
       PRINT*, "MPI aborted, status = ", ierr
#endif      

       STOP         

       RETURN
       END SUBROUTINE abort
                 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      END MODULE quit