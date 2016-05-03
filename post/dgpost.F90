      PROGRAM dgpost
      
      USE globals, ONLY:rp
      USE globalize

      IMPLICIT NONE
      
      INTEGER :: i,pe,el,dof,sta
      INTEGER :: npe,tne,mne,ne,ndof,nout_sol,nsta,mnlsta,nlsta,nout_sta
      INTEGER :: tstep
      INTEGER :: read_stat
      INTEGER :: lname
      CHARACTER(6) :: dirname
      CHARACTER(100) :: file_name
      INTEGER, ALLOCATABLE, DIMENSION(:) :: lel2gel,sta_l2g
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: H_global, Qx_global, Qy_global
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: Hsta_global, Qxsta_global, Qysta_global      
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: hb_local, hb_global
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: hbsta_global
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: t_sol,t_sta
      
      CALL version()

      OPEN(UNIT=80,FILE='PE0000/fort.80')      
      READ(80,*) npe
      CLOSE(80)
      
      
      
      

      PRINT*, "Globalizing solutions/stations"
      DO pe = 1,npe
      
        dirname = "PE0000"
        lname = 6
        
        WRITE(dirname(3:lname),"(I4.4)") pe-1
        
        OPEN(UNIT=80,FILE=dirname(1:lname)//'/'//'fort.80',POSITION='rewind')
        
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
        
        
        

        
        IF(.not. ALLOCATED(H_global))THEN
          ALLOCATE(t_sol(nout_sol+2))
          ALLOCATE(H_global(tne,ndof,nout_sol+1),Qx_global(tne,ndof,nout_sol+1),Qy_global(tne,ndof,nout_sol+1),hb_global(tne,ndof))
        ENDIF  
        
        file_name = dirname(1:lname)//'/'//'solution_H.d'
        CALL globalize_solution(file_name,ne,ndof,lel2gel,t_sol,H_global)
        
        file_name = dirname(1:lname)//'/'//'solution_Qx.d'
        CALL globalize_solution(file_name,ne,ndof,lel2gel,t_sol,Qx_global)   
        
        file_name = dirname(1:lname)//'/'//'solution_Qy.d'
        CALL globalize_solution(file_name,ne,ndof,lel2gel,t_sol,Qy_global)        

        

        IF(.not. ALLOCATED(hb_local)) THEN
          ALLOCATE(hb_local(mne,ndof))
        ENDIF        
        
        OPEN(UNIT=65,FILE=dirname(1:lname)//'/'//'hb_modal.d')
        DO dof = 1,4
          READ(65,6364) (hb_local(el,dof), el = 1,ne)
        ENDDO
          
        DO dof = 1,4
          DO el = 1,ne
            hb_global(lel2gel(el),dof) = hb_local(el,dof)
          ENDDO
        ENDDO
          
        CLOSE(65)        

        
        
        
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
        
        
        

        
        IF (.not. ALLOCATED(Hsta_global)) THEN
          ALLOCATE(t_sta(nout_sta+2))
          ALLOCATE(Hsta_global(nsta,nout_sta+1),Qxsta_global(nsta,nout_sta+1),Qysta_global(nsta,nout_sta+1),hbsta_global(nsta))
        ENDIF
        
        file_name = dirname(1:lname)//'/'//'station_H.d'
        CALL globalize_stations(file_name,nlsta,sta_l2g,t_sta,Hsta_global) 
        
        file_name = dirname(1:lname)//'/'//'station_Qx.d'
        CALL globalize_stations(file_name,nlsta,sta_l2g,t_sta,Qxsta_global) 
        
        file_name = dirname(1:lname)//'/'//'station_Qy.d'
        CALL globalize_stations(file_name,nlsta,sta_l2g,t_sta,Qysta_global)      
        
        
        OPEN(UNIT=612,FILE=dirname(1:lname)//'/'//'station_hb.d')        
        READ(612,*)
        READ(612,*) nlsta
        DO sta = 1,nlsta
          READ(612,*) hbsta_global(sta_l2g(sta))
        ENDDO        
        CLOSE(612)
        
      ENDDO
      PRINT*, "  done"      
      
      
      
      
      
      
      
      PRINT*, "Writing solutions"
      file_name(:) = ' '
      
      file_name = 'solution_H.d'
      CALL write_global_sol(file_name,nout_sol,tne,ndof,t_sol,H_global)
      
      file_name = 'solution_Qx.d'
      CALL write_global_sol(file_name,nout_sol,tne,ndof,t_sol,Qx_global)
      
      file_name = 'solution_Qy.d'
      CALL write_global_sol(file_name,nout_sol,tne,ndof,t_sol,Qy_global) 
      PRINT*, "  done"
      
      
      
      
      PRINT*, "Writing stations"
      file_name(:) = ' '  
      
      file_name = 'station_H.d'
      CALL write_global_sta(file_name,nout_sta,nsta,t_sta,Hsta_global)
      
      file_name = 'station_Qx.d'
      CALL write_global_sta(file_name,nout_sta,nsta,t_sta,Qxsta_global)
      
      file_name = 'station_Qy.d'
      CALL write_global_sta(file_name,nout_sta,nsta,t_sta,Qysta_global)      
      PRINT*, "  done"
      
      OPEN(UNIT=61,FILE='station_hb.d')
      WRITE(61,*) 'Globalized Stations'
      WRITE(61,*) nsta
      DO sta = 1,nsta
        WRITE(61,"(e24.17)") hbsta_global(sta)
      ENDDO
      CLOSE(61)      
      
      CALL system('cp PE0000/modal2nodal.d .')
      CALL system('cp PE0000/projection.d .')
      
 6364  FORMAT(160000(e24.17,1x))
      END PROGRAM dgpost
