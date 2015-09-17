      PROGRAM dgprep
      
      USE globals
      USE allocation, ONLY: sizes
      USE read_dginp
      USE messenger2, ONLY: message_init,nproc

      IMPLICIT NONE
      
      INTEGER :: pe
      REAL :: elsum,elavg,elstd

      PRINT*, "dgprep"
      PRINT*, " "            

      CALL message_init()      

      PRINT*, "Input number of processors"
      READ(*,*) nproc
      PRINT*, " "
      
      CALL read_input()
      
      CALL sizes()
      
      CALL read_grid()
      
      CALL read_forcing()
      
      CALL connect()
      
      CALL metis2(nproc)
      
      CALL decomp2()
      
      CALL write_files()
      
      elsum = 0.0
      DO pe = 1,nproc
        elsum = elsum + nresel(pe)
      ENDDO
      elavg = elsum/nproc
      
      PRINT*,"Average elements per partition: " , elavg
      
      elsum = 0.0
      DO pe = 1,nproc
        elsum = elsum + (nresel(pe)-elavg)**2
      ENDDO
      elstd = sqrt(elsum/nproc)
      
      PRINT*, "Standard deviation", elstd
      
      PRINT*, "Maximum number of elements", MAXVAL(nresel)
      PRINT*, "Minimum number of elements", MINVAL(nresel)
      PRINT*, " "
      
      END PROGRAM dgprep