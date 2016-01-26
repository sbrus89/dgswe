      MODULE output
      
      USE netcdf       
       
      IMPLICIT NONE
      INTEGER :: ncid
      INTEGER :: dimID_time,dimID_mndof,dimID_el
      INTEGER :: dimID_ldof,dimID_tdof,dimID_mnnds,dimID_neltype
      INTEGER :: zid,qxid,qyid,tid,hbid   
      INTEGER :: ml2id,mmlid,m2nid
      INTEGER, DIMENSION(3) :: sol_dim
      INTEGER, DIMENSION(2) :: hb_dim        
      INTEGER, DIMENSION(2) :: ml2_dim          
      INTEGER, DIMENSION(2) :: mml_dim
      INTEGER, DIMENSION(3) :: m2n_dim
!       INTEGER, DIMENSION(0) :: scal_dim

      INTEGER :: nsnap
       
      CONTAINS
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      

      SUBROUTINE nc_setup()

      USE globals, ONLY: ne,ndof,mndof,mnnds,ml2,mml,nel_type
      USE messenger2, ONLY: dirname,lname

      IMPLICIT NONE


      CALL check(NF90_CREATE(dirname(1:lname)//'/'//'solution.nc',NF90_CLOBBER,ncid))

      CALL check(NF90_DEF_DIM(ncid,'time',NF90_UNLIMITED,dimID_time))
      CALL check(NF90_DEF_DIM(ncid,'tri_dof',ndof(1),dimID_tdof))
      CALL check(NF90_DEF_DIM(ncid,'mndof',mndof,dimID_mndof))   
      CALL check(NF90_DEF_DIM(ncid,'mnnds',mnnds,dimID_mnnds))            
      CALL check(NF90_DEF_DIM(ncid,'el',ne,dimID_el))
      CALL check(NF90_DEF_DIM(ncid,'linear_dof',3,dimID_ldof))
      CALL check(NF90_DEF_DIM(ncid,'nel_type',nel_type,dimID_neltype))          
          
      sol_dim = (/ dimID_el, dimID_mndof, dimID_time /)
        
      CALL check(NF90_DEF_VAR(ncid,'Z',NF90_DOUBLE,sol_dim,zid))
      CALL check(NF90_DEF_VAR(ncid,'Qx',NF90_DOUBLE,sol_dim,qxid))
      CALL check(NF90_DEF_VAR(ncid,'Qy',NF90_DOUBLE,sol_dim,qyid))
      CALL check(NF90_DEF_VAR(ncid,'t',NF90_DOUBLE,dimID_time,tid))
          
      hb_dim = (/ dimID_el, dimID_mnnds /)
      CALL check(NF90_DEF_VAR(ncid,'hb',NF90_DOUBLE,hb_dim,hbid))          
          
      m2n_dim = (/ dimID_mnnds, dimID_mndof, dimID_neltype/)
      CALL check(NF90_DEF_VAR(ncid,'m2n',NF90_DOUBLE,m2n_dim,m2nid))
          
      ml2_dim = (/ dimID_ldof, dimID_tdof /)        
      CALL check(NF90_DEF_VAR(ncid,'ml2',NF90_DOUBLE,ml2_dim,ml2id))
          
      mml_dim = (/ dimID_ldof, dimID_ldof /)          
      CALL check(NF90_DEF_VAR(ncid,'mml',NF90_DOUBLE,mml_dim,mmlid))
          
      CALL check(NF90_ENDDEF(ncid))
          
          
      CALL check(NF90_PUT_VAR(ncid,ml2id,ml2))
      CALL check(NF90_PUT_VAR(ncid,mmlid,mml))   
          
      nsnap = 1
          
          
      END SUBROUTINE nc_setup
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
        
      SUBROUTINE write_output(init)
      
      USE globals, ONLY: rp,t,mndof,ne, &
                         Hwrite,Zwrite,Qxwrite,Qywrite, &
                         Znc,Qxnc,Qync
      USE messenger2, ONLY: myrank
      USE read_dginp, ONLY: out_direc,grid_file,dt,tf

      IMPLICIT NONE

      INTEGER :: dof,el
      LOGICAL :: init
      
      INTEGER, DIMENSION(3) :: var_start,var_end
      REAL(rp), DIMENSION(1) :: t_tmp
     

      IF (init) THEN
      
        ! Initialize files and write initial condition      
      
        IF (myrank == 0) THEN
          PRINT "(A)", "---------------------------------------------"
          PRINT "(A)", "               Time Stepping                 "
          PRINT "(A)", "---------------------------------------------"
          PRINT "(A)", " "

          PRINT "(A,e12.4)", "Time step: ",dt
          PRINT "(A,e12.4)", "Final time: ",tf

          PRINT "(A)", " "
        ENDIF      

        OPEN(unit=63,file=trim(out_direc) // 'solution_H.d')
        OPEN(unit=641,file=trim(out_direc) // 'solution_Qx.d')
        OPEN(unit=642,file=trim(out_direc) // 'solution_Qy.d')

        WRITE(63,"(A)") grid_file
        WRITE(641,"(A)") grid_file
        WRITE(642,"(A)") grid_file
        
        ! Set up netcdf output files
        CALL nc_setup()       

      ELSE

        IF(myrank == 0) THEN
          PRINT("(A,e15.8)"), 't = ', t
        ENDIF
        
                  
      ENDIF



      WRITE(63,"(e24.17)") t
      DO dof = 1,mndof
        WRITE(63,"(16000(e24.17,1x))") (Zwrite(el,dof)%ptr, el = 1,ne)
      ENDDO

      WRITE(641,"(e24.17)") t
      DO dof = 1,mndof
        WRITE(641,"(16000(e24.17,1x))") (Qxwrite(el,dof)%ptr, el = 1,ne)
      ENDDO

      WRITE(642,"(e24.17)") t
      DO dof = 1,mndof
        WRITE(642,"(16000(e24.17,1x))") (Qywrite(el,dof)%ptr, el = 1,ne)
      ENDDO     
      
      
      
      
      
      DO dof = 1,mndof
        DO el = 1,ne
          Znc(el,dof)  = Zwrite(el,dof)%ptr
          Qxnc(el,dof) = Qxwrite(el,dof)%ptr
          Qync(el,dof) = Qywrite(el,dof)%ptr
        ENDDO
      ENDDO
      
      
      
      var_start = (/ 1, 1, nsnap /)
      var_end = (/ ne, mndof, 1 /)
      
      t_tmp(1) = t
      CALL check(NF90_PUT_VAR(ncid,tid,t_tmp,(/nsnap/),(/1/)))
      CALL check(NF90_PUT_VAR(ncid,zid,Znc,var_start,var_end ))
      CALL check(NF90_PUT_VAR(ncid,qxid,Qxnc,var_start,var_end ))
      CALL check(NF90_PUT_VAR(ncid,qyid,Qync,var_start,var_end ))

      nsnap = nsnap + 1
      

      RETURN
      END SUBROUTINE
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                    
      
      SUBROUTINE close_output()     
      
      IMPLICIT NONE
      
      
      CLOSE(63)
      CLOSE(641)
      CLOSE(642)
      
      CALL check(NF90_CLOSE(ncid))           
      
      RETURN
      END SUBROUTINE        
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
                
      SUBROUTINE check(status)

      IMPLICIT NONE
      INTEGER :: status

      IF(status /= NF90_NOERR) THEN
        PRINT("(A,A)"), "fatal error from ", TRIM(NF90_STRERROR(status))  
      ENDIF
               
      END SUBROUTINE check

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
      END MODULE output
