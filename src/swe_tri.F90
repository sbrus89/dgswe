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
      


      coord_sys = 1
      slam0 = -79.0d0*deg2rad
      sphi0 = 35.0d0*deg2rad
      h0 = 0d0
      
      CALL message_init()
      
      CALL read_input()

      CALL sizes()

      ! Read in grid file
      CALL read_grid()

      ! Read in forcing file
      CALL read_forcing()
      
      CALL read_message_files()
   
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
      
!       STOP
      
      ! Set up netcdf output files
!       CALL file_setup()

      ! Compute initial condition
      CALL initial()

      ! Boundary forcing interpolation
      CALL interp_forcing()
      
      CALL metis2()
      
      CALL decomp2()
      
      CALL message_setup()
      
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
      
      IF (myrank == 0) THEN
        PRINT "(A)", "---------------------------------------------"
        PRINT "(A)", "               Time Stepping                 "
        PRINT "(A)", "---------------------------------------------"
        PRINT "(A)", " "

        PRINT "(A,e12.4)", "Time step: ",dt
        PRINT "(A,e12.4)", "Final time: ",tf

        PRINT "(A)", " "
      ENDIF
      
      tf = tf*86400d0
      tstep = int(tf/dt)
      tskp = int(tf/(lines*dt))     
      
      ark(1) = 0d0
      ark(2) = - 567301805773d0 / 1357537059087d0
      ark(3) = -2404267990393d0 / 2016746695238d0
      ark(4) = -3550918686646d0 / 2091501179385d0
      ark(5) = -1275806237668d0 / 842570457699d0
      
      brk(1) = 1432997174477d0 / 9575080441755d0
      brk(2) = 5161836677717d0 / 13612068292357d0
      brk(3) = 1720146321549d0 / 2090206949498d0
      brk(4) = 3134564353537d0 / 4481467310338d0
      brk(5) = 2277821191437d0 / 14882151754819d0
      
      crk(1) = 0d0
      crk(2) = 1432997174477d0 / 9575080441755d0
      crk(3) = 2526269341429d0 / 6820363962896d0
      crk(4) = 2006345519317d0 / 3224310063776d0
      crk(5) = 2802321613138d0 / 2924317926251d0           

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
