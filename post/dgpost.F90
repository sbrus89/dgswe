      PROGRAM dgpost
      
      USE globals, ONLY:rp
      USE read_write_output
      USE version
      USE read_dginp

      IMPLICIT NONE
      
      INTEGER :: i,pe,el,dof,sta
      INTEGER :: npe,tne,mne,ne,ndof,nout_sol,nsta,mnlsta,nlsta,nout_sta
      INTEGER :: tstep
      INTEGER :: read_stat
      INTEGER :: lname
      CHARACTER(6) :: dirname
      CHARACTER(7) :: pe_direc
      CHARACTER(100) :: file_name
      INTEGER, ALLOCATABLE, DIMENSION(:) :: lel2gel,sta_l2g
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: Z_global,Qx_global,Qy_global,hb_global
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: Zsta_global,Qxsta_global,Qysta_global,hbsta_global      
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: t_sol,t_sta
      LOGICAL :: post=.true.
      
      CALL version_information(unit=6)
      
      CALL read_input(0,".")
      

      OPEN(UNIT=80,FILE='PE0000/fort.80')      
      READ(80,*) npe
      CLOSE(80)
      
      
      
      

      PRINT*, "Globalizing solutions/stations"
      DO pe = 1,npe
      
        dirname = "PE0000"
        lname = 6
        
        WRITE(dirname(3:lname),"(I4.4)") pe-1
        pe_direc = dirname(1:lname)//'/'
        
        OPEN(UNIT=80,FILE=pe_direc//'fort.80',POSITION='rewind')
        
        READ(80,*)
        READ(80,*) tne
        READ(80,*) mne
        READ(80,*) ne
        READ(80,*) ndof
        READ(80,*) nout_sol
        
        IF(.not. ALLOCATED(lel2gel)) THEN
          ALLOCATE(lel2gel(mne))
        ENDIF
        
        DO el = 1,ne
          READ(80,*) i, lel2gel(el)
        ENDDO
        
        CLOSE(80)   
        
        
        

        
        IF(.not. ALLOCATED(Z_global))THEN
          ALLOCATE(t_sol(nout_sol+1))
          ALLOCATE(Z_global(tne,ndof,nout_sol+1))
          ALLOCATE(Qx_global(tne,ndof,nout_sol+1))
          ALLOCATE(Qy_global(tne,ndof,nout_sol+1))
          ALLOCATE(hb_global(tne,ndof,1))
        ENDIF  
        
        CALL read_solution_full(pe_direc,'hb.sol',t_sol,hb_global,"T",lel2gel)
        CALL read_solution_full(pe_direc,'Z.sol',t_sol,Z_global,"T",lel2gel)       
        CALL read_solution_full(pe_direc,'Qx.sol',t_sol,Qx_global,"T",lel2gel)           
        CALL read_solution_full(pe_direc,'Qy.sol',t_sol,Qy_global,"T",lel2gel)        


        
        
        
        OPEN(UNIT=82,FILE=dirname(1:lname)//'/'//'fort.82')        
        READ(82,*) 
        READ(82,*) nsta
        READ(82,*) mnlsta
        READ(82,*) nlsta
        READ(82,*) nout_sta
        
        IF (.not. ALLOCATED(sta_l2g)) THEN
          ALLOCATE(sta_l2g(nsta))
        ENDIF
        
        DO sta = 1,nlsta
          READ(82,*) sta_l2g(sta)
        ENDDO
        CLOSE(82)
        
        
        

        
        IF (.not. ALLOCATED(Zsta_global)) THEN
          ALLOCATE(t_sta(nout_sta+1))
          ALLOCATE(Zsta_global(nsta,1,nout_sta+1))
          ALLOCATE(Qxsta_global(nsta,1,nout_sta+1))
          ALLOCATE(Qysta_global(nsta,1,nout_sta+1))
          ALLOCATE(hbsta_global(nsta,1,1))
        ENDIF
        
        CALL read_solution_full(pe_direc,'hb.sta',t_sta,hbsta_global,"N",sta_l2g)
        CALL read_solution_full(pe_direc,'Z.sta',t_sta,Zsta_global,"N",sta_l2g) 
        CALL read_solution_full(pe_direc,'Qx.sta',t_sta,Qxsta_global,"N",sta_l2g) 
        CALL read_solution_full(pe_direc,'Qy.sta',t_sta,Qysta_global,"N",sta_l2g)      
        
        
        
      ENDDO
      PRINT*, "  done"      
      
      
      
      
      
      
      
      PRINT*, "Writing solutions"
      CALL write_solution_full(out_direc,'hb.sol',ndof,tne,1,t_sol,hb_global,"T",post)      
      CALL write_solution_full(out_direc,'Z.sol',ndof,tne,nout_sol+1,t_sol,Z_global,"T",post)
      CALL write_solution_full(out_direc,'Qx.sol',ndof,tne,nout_sol+1,t_sol,Qx_global,"T",post)      
      CALL write_solution_full(out_direc,'Qy.sol',ndof,tne,nout_sol+1,t_sol,Qy_global,"T",post) 
      PRINT*, "  done"

      
      PRINT*, "Writing stations"
      CALL write_solution_full(out_direc,'hb.sta',nsta,1,1,t_sta,hbsta_global,"N",post)      
      CALL write_solution_full(out_direc,'Z.sta',nsta,1,nout_sta+1,t_sta,Zsta_global,"N",post)
      CALL write_solution_full(out_direc,'Qx.sta',nsta,1,nout_sta+1,t_sta,Qxsta_global,"N",post)
      CALL write_solution_full(out_direc,'Qy.sta',nsta,1,nout_sta+1,t_sta,Qysta_global,"N",post)      
      PRINT*, "  done"

      
      CALL system('cp PE0000/modal2nodal.d .')
      CALL system('cp PE0000/projection.d .')
      
      END PROGRAM dgpost
