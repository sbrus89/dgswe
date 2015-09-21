      PROGRAM dgprep
      
      USE globals
      USE allocation, ONLY: sizes, alloc_trans_arrays
      USE read_dginp
      USE messenger2, ONLY: message_init,nproc
      USE vandermonde, ONLY: area_vandermonde
      USE shape_functions, ONLY: shape_functions_qpts

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
      
      CALL edge_qpts()
      
      mnqpta = 0
      nqpta = 0 
      
      CALL alloc_trans_arrays()
      
      CALL area_vandermonde()
      
      CALL shape_functions_qpts()
      
      CALL normals()
      
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