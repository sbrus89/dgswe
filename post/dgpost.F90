      PROGRAM dgpost
      
      USE globals, ONLY:rp
      USE read_write_output
      USE version
      USE read_dginp

      IMPLICIT NONE
      
      INTEGER :: i,pe,el,dof,sta
      INTEGER :: npe,tne,mne,ne,ndof,nout_sol,nsta,mnlsta,nlsta,nout_sta
      INTEGER :: tstep
      INTEGER :: nsnap_sol,nsnap_sta
      INTEGER :: read_stat
      INTEGER :: lname
      CHARACTER(6) :: dirname
      CHARACTER(7) :: pe_direc
      CHARACTER(100) :: file_name
      INTEGER, ALLOCATABLE, DIMENSION(:) :: lel2gel,sta_l2g
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: Z_global,Qx_global,Qy_global,hb_global
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: Zsta_global,Qxsta_global,Qysta_global,hbsta_global      
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: t_sol,t_sta
      INTEGER :: ndof_hb      
      LOGICAL :: post=.true.
      
      CALL version_information(unit=6)
      
      CALL read_input(0,".")
      
      ndof_hb = (hbp+1)*(hbp+1)         

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
          ALLOCATE(hb_global(tne,ndof_hb,1))
        ENDIF  
        
        nsnap_sol = 1
        CALL read_solution_full(pe_direc,'hb.sol',"T",t_sol,hb_global,nsnap_sol,lel2gel)
        nsnap_sol = 9999
        CALL read_solution_full(pe_direc,'Z.sol',"T",t_sol,Z_global,nsnap_sol,lel2gel)       
        CALL read_solution_full(pe_direc,'Qx.sol',"T",t_sol,Qx_global,nsnap_sol,lel2gel)           
        CALL read_solution_full(pe_direc,'Qy.sol',"T",t_sol,Qy_global,nsnap_sol,lel2gel)        


        
        
        IF (sta_opt > 0) THEN
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
        
          nsnap_sta = 1
          CALL read_solution_full(pe_direc,'hb.sta',"N",t_sta,hbsta_global,nsnap_sta,sta_l2g)
          nsnap_sta = 9999
          CALL read_solution_full(pe_direc,'Z.sta',"N",t_sta,Zsta_global,nsnap_sta,sta_l2g) 
          CALL read_solution_full(pe_direc,'Qx.sta',"N",t_sta,Qxsta_global,nsnap_sta,sta_l2g) 
          CALL read_solution_full(pe_direc,'Qy.sta',"N",t_sta,Qysta_global,nsnap_sta,sta_l2g)
        ENDIF
        
        
        
      ENDDO
      PRINT*, "  done"      
      
      
      
      
      
      
      
      PRINT*, "Writing solutions"
      CALL write_solution_full(out_direc,'hb.sol',ndof_hb,tne,1,"T",t_sol,hb_global,post)      
      CALL write_solution_full(out_direc,'Z.sol',ndof,tne,nsnap_sol,"T",t_sol,Z_global,post)
      CALL write_solution_full(out_direc,'Qx.sol',ndof,tne,nsnap_sol,"T",t_sol,Qx_global,post)      
      CALL write_solution_full(out_direc,'Qy.sol',ndof,tne,nsnap_sol,"T",t_sol,Qy_global,post) 
      PRINT*, "  done"

      IF (sta_opt > 0) THEN
        PRINT*, "Writing stations"
        CALL write_solution_full(out_direc,'hb.sta',nsta,1,1,"N",t_sta,hbsta_global,post)      
        CALL write_solution_full(out_direc,'Z.sta',nsta,1,nsnap_sta,"N",t_sta,Zsta_global,post)
        CALL write_solution_full(out_direc,'Qx.sta',nsta,1,nsnap_sta,"N",t_sta,Qxsta_global,post)
        CALL write_solution_full(out_direc,'Qy.sta',nsta,1,nsnap_sta,"N",t_sta,Qysta_global,post)      
        PRINT*, "  done"
      ENDIF

      
      CALL system('cp PE0000/modal2nodal.d .')
      CALL system('cp PE0000/projection.d .')
      
      END PROGRAM dgpost
