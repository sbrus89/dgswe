      PROGRAM swe_tri

!$    USE omp_lib      
      USE globals
      USE basis
      USE allocation
      USE messenger2
      USE read_dginp

      IMPLICIT NONE
      INTEGER :: it,tskp,cnt
      INTEGER :: i,j,et
      REAL(rp) :: tstep,t_start,t_end
      


      
      ! Initialize MPI and determine rank if appropriate
      CALL message_init()
      
      ! Read in the keyword (or fixed) format input file 
      CALL read_input()

      ! Compute # of dofs and # of nodes straight/curved elements etc.
      CALL sizes()

      ! Read in grid file
      CALL read_grid()

      ! Read in forcing file
      CALL read_forcing()     
   
      ! Find edge connectivity
      CALL connect()

      ! Get quadrature points for area integration
      CALL area_qpts()

      ! Calculate basis function values at edge quadrature points
      CALL edge_qpts()   
      
      ! Calculate basis function and derivative values at area quadrature points
      CALL area_basis()  

      ! Calculate basis function and derivative values at edge quadrature points
      CALL edge_basis()            
      
      ! Compute element area, edge length, edge normals, and bathymetry derivatives
      CALL element_data()           
      
      ! Set up netcdf output files
!       CALL file_setup()

      ! Compute initial condition
      CALL initial()

      ! Boundary forcing interpolation
      CALL interp_forcing()
      
      ! Partition domain for element/edge blocking
      CALL metis2(npart)
      
      ! Decompose domain and prep element/edge blocking
      CALL edge_partition2()
      
      ! Read the local-to-global element and message passing files
      CALL read_message_files()      
      
      ! Set up send/recieve buffers and edge data structures
      CALL message_setup()
      
      ! Initialize the MPI persistent message passing calls
      CALL communication_setup()    

      OPEN(unit=195,file=trim(out_direc) // 'edge_check.d')        

      DO et = 1,2*nel_type
        DO i = 1,mnnds
          WRITE(195,"(100(ES24.17,1x))") (dpsids(i,j,et), j = 1,mnqpta+4*mnqpte)
        ENDDO
      WRITE(195,*) " "
      ENDDO      
      
#ifdef openmp
      t_start = omp_get_wtime()
#else
      CALL CPU_TIME(t_start)
#endif
      

      
      tf = tf*86400d0
      tstep = int(tf/dt)
      tskp = int(tf/(lines*dt))             

      t = 0d0
      cnt = 0
 
      CALL write_output(.true.)

      DO it = 1,tstep

        CALL rk()

         t = t + dt

         cnt = cnt + 1
         IF(cnt == tskp) THEN

           CALL write_output(.false.)             
           cnt = 0

         ENDIF

      ENDDO

#ifdef openmp      
      t_end = omp_get_wtime()
#else
      CALL CPU_TIME(t_end)
#endif

      PRINT*, ' '
      PRINT("(A,F25.5,A)"), "CPU time = ",t_end-t_start," seconds"
      
!       OPEN(UNIT=101, FILE='CPUtime.log', STATUS='OLD', POSITION='APPEND')
!       WRITE(101,*) "grid = ", grid_file
!       WRITE(101,*) "p = ", p
!       WRITE(101,*) "dt, tf = ", dt,tf
! #ifdef rk22
!       WRITE(101,*) "RK 22"
! #endif
! #ifdef rk33
!       WRITE(101,*) "RK 33"
! #endif
!       WRITE(101,*) "nblk, npart = ", nblk, npart
!       WRITE(101,*) "mnpartel,mnparted = ", mnpartel,mnparted   
!       WRITE(101,*) "CPU time = ", t_end-t_start
!       WRITE(101,*) " "
!       CLOSE(101)      
      
      CALL finish()
      
      END PROGRAM swe_tri
