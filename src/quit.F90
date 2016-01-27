      MODULE quit 

#ifdef CMPI        
      USE mpi
#endif      

      USE globals, ONLY: rp
      
      IMPLICIT NONE

      CONTAINS 
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       SUBROUTINE end_time(t_start,nproc)

       IMPLICIT NONE
       
       INTEGER :: i
       INTEGER :: ierr
       INTEGER :: nproc
       REAL(rp) :: t_start,t_end       
       REAL(rp) :: t_max,t_min,t_avg
       REAL(rp) :: cpu_times(nproc)
       

#ifdef openmp      
      t_end = omp_get_wtime()
#else
      CALL CPU_TIME(t_end)     
#endif      
       
#ifdef CMPI       
      CALL MPI_GATHER(t_end-t_start,1,MPI_DOUBLE_PRECISION,cpu_times,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      
      IF (myrank == 0) THEN
        t_avg = 0d0
        DO i = 1,nproc
          t_avg = t_avg + cpu_times(i)
        ENDDO
        t_avg = t_avg/nproc
      
        t_max = MAXVAL(cpu_times)
        t_min = MINVAL(cpu_times)
      
        PRINT*, ' '      
        PRINT("(A,F25.5,A)"), "Average CPU time = ",t_avg," seconds"
        PRINT("(A,F25.5,A)"), "Minimum CPU time = ",t_min," seconds"
        PRINT("(A,F25.5,A)"), "Maximum CPU time = ",t_max," seconds" 
      ENDIF   
      
#else      
      PRINT*, ' '      
      PRINT("(A,F25.5,A)"), "CPU time = ",t_end-t_start," seconds"      
#endif




       RETURN
       END SUBROUTINE end_time

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       SUBROUTINE finish()          

       IMPLICIT NONE
       
       INTEGER :: ierr       
       
#ifdef CMPI       
       CALL MPI_FINALIZE(ierr)
       
       IF(myrank == 0) THEN
         PRINT*, "MPI terminated, status = ", ierr
       ENDIF
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