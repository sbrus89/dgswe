      PROGRAM dgprep
      
      USE globals
      USE allocation, ONLY: sizes, alloc_trans_arrays
      USE read_dginp
      USE messenger2, ONLY: message_init,nproc,nied_pe
      USE edge_qpts_mod, ONLY: edge_qpts
      USE grid_file_mod, ONLY: courant,read_stations

      IMPLICIT NONE
      
      INTEGER :: pe,el
      INTEGER :: et,ed,eln
      INTEGER :: ncalc,ncalc_area,ncalc_edge,ncalc_solve
      REAL :: elsum,elavg,elstd
      REAL(rp) :: cfl,u       

      PRINT*, "dgprep"
      PRINT*, " "            

      CALL message_init()      

      PRINT*, "Input number of processors"
      READ(*,*) nproc
      PRINT*, " "
      
      CALL read_input(0,".")
      
      CALL sizes()
      
      CALL read_grid()   
      
      cfl = 1d0
      u = 0d0
      CALL courant(p,ne,u,cfl,el_type,nverts,nnds,elhb,el_size)
      
      CALL read_forcing()              
      
      CALL connect()
      
      IF (sta_opt > 0) THEN
        CALL read_stations(0,stations_file,sta_opt,nsta,xysta)
        CALL find_stations()
      ENDIF       
      
      CALL edge_qpts(0,p,ctp,nel_type,nqpte,mnqpte,wpte,qpte)
      
      mnqpta = 0
      nqpta = 0 
      
      CALL alloc_trans_arrays()
      
      CALL shape_functions_area_qpts()
      
      CALL shape_functions_edge_qpts()
      
      CALL curvilinear()          
      
      CALL normals()
      
      CALL edge_transformation()      
      
      CALL bathymetry_interp_edge_qpts()
      
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
      
      ncalc_area = 0
      ncalc_edge = 0
      ncalc_solve = 0
      DO pe = 1,nproc
        DO el = 1,nresel(pe) 
          eln = el_l2g(el,pe)        
          et = el_type(eln)
          IF (et == 1) THEN
            ncalc_area = ncalc_area + 12
          ELSE
            ncalc_area = ncalc_area + 36
          ENDIF
          
        ENDDO
        
        DO ed = 1,nied_pe(pe)
          ncalc = ncalc + 4
        ENDDO
        
        
        
      ENDDO
      
      END PROGRAM dgprep
