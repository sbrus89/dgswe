      PROGRAM swe_tri

!$    USE omp_lib      
      USE globals
      USE basis

      IMPLICIT NONE
      INTEGER :: it,tskp,cnt,myid,mndof
      INTEGER :: i,j
      REAL(pres) :: tstep,t_start,t_end

      OPEN(unit=63,file='../output/solution_H.d')
      OPEN(unit=641,file='../output/solution_Qx.d')
      OPEN(unit=642,file='../output/solution_Qy.d')
      OPEN(unit=17,file='../output/edge_connect.d')

      PRINT*, ' '

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      nrblk = 1

      CALL read_input()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ndof(1) = (p+1)*(p+2)/2
      ndof(2) = (p+1)**2
      ndof(3) = ndof(1)
      ndof(4) = ndof(2)
      
      mndof = maxval(ndof)
      
      nverts(1) = 3
      nverts(2) = 4
      nverts(3) = 3
      nverts(4) = 4
      
      tstep = int(tf/dt)
      tskp = int(tf/(lines*dt)) 
      
#ifdef openmp      
      PRINT*, "Thread numbers: "
!$OMP  parallel private(myid)
      myid = omp_get_thread_num()
      PRINT*, myid
      nrblk = omp_get_num_threads()
!$OMP end parallel
#endif

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

      ! Allocate arrays needed in time-stepping and rhs evaluation
      CALL alloc_arrays()

      CALL ptr_arrays()
      
      ! Set up netcdf output files
!       CALL file_setup()

      ! Compute initial condition, boundary forcing interpolation
      CALL initial()
      
      CALL metis2()
      
      CALL decomp2()
      
      
#ifdef openmp
      t_start = omp_get_wtime()
#else
      CALL CPU_TIME(t_start)
#endif
      
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", "               Time Stepping                 "
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", " "

      PRINT "(A,e12.4)", "Time step: ",dt
      PRINT "(A,e12.4)", "Final time: ",tf

      PRINT "(A)", " "

      t = 0d0
      cnt = 0
      DO it = 1,tstep

        CALL rk()

         t = t + dt

         cnt = cnt + 1
         IF(cnt == tskp) THEN

           PRINT("(A,e15.8)"), 't = ', t

           WRITE(63,"(e24.17)") t
           DO dof = 1,mndof
             WRITE(63,"(16000(e24.17,1x))") (Hwrite(el,dof)%ptr, el = 1,ne)
           ENDDO

           WRITE(641,"(e24.17)") t
           DO dof = 1,mndof
             WRITE(641,"(16000(e24.17,1x))") (Qxwrite(el,dof)%ptr, el = 1,ne)
           ENDDO

           WRITE(642,"(e24.17)") t
           DO dof = 1,mndof
             WRITE(642,"(16000(e24.17,1x))") (Qywrite(el,dof)%ptr, el = 1,ne)
           ENDDO           
             
           cnt = 0

         ENDIF

      ENDDO

#ifdef openmp      
      t_end = omp_get_wtime()
#else
      CALL CPU_TIME(t_end)
#endif

      PRINT*, ' '
      PRINT("(A,F10.5,A)"), "CPU time = ",t_end-t_start," seconds"
      
      OPEN(UNIT=101, FILE='CPUtime.log', STATUS='OLD', POSITION='APPEND')
      WRITE(101,*) "grid = ", grid_file
      WRITE(101,*) "p = ", p
      WRITE(101,*) "dt, tf = ", dt,tf
#ifdef rk22
      WRITE(101,*) "RK 22"
#endif
#ifdef rk33
      WRITE(101,*) "RK 33"
#endif
      WRITE(101,*) "nblk, npart = ", nblk, npart
      WRITE(101,*) "mnpartel,mnparted = ", mnpartel,mnparted   
      WRITE(101,*) "CPU time = ", t_end-t_start
      WRITE(101,*) " "
      CLOSE(101)
      

      PRINT*, ' '
      END PROGRAM swe_tri
