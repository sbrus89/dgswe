      PROGRAM dgpost
      
      USE netcdf

      IMPLICIT NONE
      
      INTEGER :: pe,el,dof,l
      INTEGER :: npe,tne,mne,ne,ndof,lines
      INTEGER :: tstep
      INTEGER :: read_stat
      INTEGER, PARAMETER :: pres = kind(1d0)
      INTEGER :: lname
      CHARACTER(6) :: dirname
      INTEGER, ALLOCATABLE, DIMENSION(:) :: lel2gel
      REAL(pres), ALLOCATABLE, DIMENSION(:,:,:) :: H_local, Qx_local, Qy_local
      REAL(pres), ALLOCATABLE, DIMENSION(:,:,:) :: H_global, Qx_global, Qy_global
      REAL(pres), ALLOCATABLE, DIMENSION(:) :: t
      REAL(pres), ALLOCATABLE, DIMENSION(:,:) :: ml2,mml
      
      INTEGER :: global_ncid, local_ncid
      INTEGER :: dimID_time,dimID_dof,dimID_el,dimID_ldof
      INTEGER :: global_hid,global_qxid,global_qyid,global_tid
      INTEGER :: local_hid,local_qxid,local_qyid,local_tid
      INTEGER :: global_ml2id,global_mmlid     
      INTEGER :: local_ml2id,local_mmlid        
      INTEGER, DIMENSION(3) :: sol_dim
      INTEGER, DIMENSION(2) :: ml2_dim          
      INTEGER, DIMENSION(2) :: mml_dim
      INTEGER, DIMENSION(3) :: var_start,var_end      

      OPEN(UNIT=19,FILE='PE0000/fort.80')      
      READ(19,*) npe
      READ(19,*) tne
      READ(19,*) mne
      READ(19,*) ne
      READ(19,*) ndof
      READ(19,*) lines      
      CLOSE(19)
      
      PRINT*, "Total number of sub-domains", npe
      PRINT*, "Number of global elements", tne
      
      ALLOCATE(ml2(3,ndof),mml(3,3))
      
      CALL check(NF90_CREATE('global_solution.nc',NF90_CLOBBER,global_ncid))
      
      CALL check(NF90_DEF_DIM(global_ncid,'time',NF90_UNLIMITED,dimID_time))
      CALL check(NF90_DEF_DIM(global_ncid,'dof',ndof,dimID_dof))
      CALL check(NF90_DEF_DIM(global_ncid,'el',tne,dimID_el))
      CALL check(NF90_DEF_DIM(global_ncid,'ldof',3,dimID_ldof)) 
      
      sol_dim = (/ dimID_el, dimID_dof, dimID_time /)
          
      CALL check(NF90_DEF_VAR(global_ncid,'H_global',NF90_DOUBLE,sol_dim,global_hid))
      CALL check(NF90_DEF_VAR(global_ncid,'Qx_global',NF90_DOUBLE,sol_dim,global_qxid))
      CALL check(NF90_DEF_VAR(global_ncid,'Qy_global',NF90_DOUBLE,sol_dim,global_qyid))
      CALL check(NF90_DEF_VAR(global_ncid,'t_global',NF90_DOUBLE,dimID_time,global_tid))
          
      ml2_dim = (/ dimID_ldof, dimID_dof /)
          
      CALL check(NF90_DEF_VAR(global_ncid,'ml2',NF90_DOUBLE,ml2_dim,global_ml2id))
          
      mml_dim = (/ dimID_ldof, dimID_ldof /)
          
      CALL check(NF90_DEF_VAR(global_ncid,'mml',NF90_DOUBLE,mml_dim,global_mmlid))
          
      CALL check(NF90_ENDDEF(global_ncid))      
      
      
      ALLOCATE(lel2gel(mne))
      ALLOCATE(H_local(mne,ndof,lines+1),Qx_local(mne,ndof,lines+1),Qy_local(mne,ndof,lines+1))  
        
      ALLOCATE(t(lines+1))
      ALLOCATE(H_global(tne,ndof,lines+1),Qx_global(tne,ndof,lines+1),Qy_global(tne,ndof,lines+1))
      
      

      PRINT*, "Globalizing solutions"
      DO pe = 1,npe
      
        dirname = "PE0000"
        lname = 6
        
        WRITE(dirname(3:lname),"(I4.4)") pe-1
        
        OPEN(UNIT=19,FILE=dirname(1:lname)//'/'//'fort.80',POSITION='rewind')
        
        READ(19,*)
        READ(19,*) tne
        READ(19,*) mne
        READ(19,*) ne
        READ(19,*) ndof
        READ(19,*) lines
        
        DO el = 1,ne
          READ(19,*) lel2gel(el)
        ENDDO
        
        CLOSE(19)   
        
        
       
       
        
!         var_start = (/1, 1, 1/)
!         var_end = (/ne, ndof, lines/)
        
        CALL check(NF90_OPEN(dirname(1:lname)//'/'//'solution.nc',NF90_NOWRITE,local_ncid))

        CALL check(NF90_INQ_VARID(local_ncid,'H',local_hid))
        CALL check(NF90_INQ_VARID(local_ncid,'Qx',local_qxid))
        CALL check(NF90_INQ_VARID(local_ncid,'Qy',local_qyid))

        CALL check(NF90_GET_VAR(local_ncid,local_hid,H_local(1:ne,1:ndof,1:lines+1)))
                
        DO l = 1,lines+1
          DO dof = 1,ndof
            DO el = 1,ne
              H_global(lel2gel(el),dof,l) = H_local(el,dof,l)
            ENDDO
          ENDDO
        ENDDO
        
        CALL check(NF90_GET_VAR(local_ncid,local_qxid,Qx_local(1:ne,1:ndof,1:lines+1)))        
                 
        DO l = 1,lines+1
          DO dof = 1,ndof
            DO el = 1,ne
              Qx_global(lel2gel(el),dof,l) = Qx_local(el,dof,l)
            ENDDO
          ENDDO
        ENDDO
 
        CALL check(NF90_GET_VAR(local_ncid,local_qyid,Qy_local(1:ne,1:ndof,1:lines+1)))   
 
        DO l = 1,lines+1
          DO dof = 1,ndof
            DO el = 1,ne
              Qy_global(lel2gel(el),dof,l) = Qy_local(el,dof,l)
            ENDDO
          ENDDO
        ENDDO
        
        CALL check(NF90_CLOSE(local_ncid))

        
      ENDDO
      PRINT*, "  done"      
      
      
      
      PRINT*, "Writing solutions"
      
      var_start = (/1, 1, 1/)
      var_end = (/tne, ndof, lines+1/)
        
      CALL check(NF90_PUT_VAR(global_ncid,global_hid,H_global,var_start,var_end))
      CALL check(NF90_PUT_VAR(global_ncid,global_qxid,Qx_global,var_start,var_end))
      CALL check(NF90_PUT_VAR(global_ncid,global_qyid,Qy_global,var_start,var_end))      
      
      
      CALL check(NF90_OPEN('PE0000/solution.nc',NF90_NOWRITE,local_ncid))
      CALL check(NF90_INQ_VARID(local_ncid,'ml2',local_ml2id))
      CALL check(NF90_INQ_VARID(local_ncid,'mml',local_mmlid))  
      CALL check(NF90_INQ_VARID(local_ncid,'t',local_tid))
      CALL check(NF90_GET_VAR(local_ncid,local_ml2id,ml2(:,:)))
      CALL check(NF90_GET_VAR(local_ncid,local_mmlid,mml(:,:)))
      CALL check(NF90_GET_VAR(local_ncid,local_tid,t(:)))
      CALL check(NF90_PUT_VAR(global_ncid,global_ml2id,ml2(:,:)))
      CALL check(NF90_PUT_VAR(global_ncid,global_mmlid,mml(:,:)))
      CALL check(NF90_PUT_VAR(global_ncid,global_tid,t(:)))
      
      CALL check(NF90_CLOSE(local_ncid))
      PRINT*, "  done"
      
          
!      CALL system('cp PE0000/projection.d .')
      
      CALL check(NF90_CLOSE(global_ncid))
      
      END PROGRAM dgpost

      
      
        SUBROUTINE check(status)

          USE netcdf
        
          IMPLICIT NONE
          INTEGER :: status

          IF(status /= NF90_NOERR) THEN
            PRINT("(A,A)"), "fatal error from ", TRIM(NF90_STRERROR(status))  
          ENDIF
               
        END SUBROUTINE check      
