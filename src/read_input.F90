      MODULE read_dginp
      
      ! Module containing subroutines to read the fort.dg file
      !
      !  - Supports both fixed and keyword formats
      !
      !  - Keyword format is intended to increase flexibility in adding/depreciating features
      !    while maintaining forward compatibility and some degree of backward compatibility 
      !    for the fort.dg file.  
      !
      !      * The keywords are configured in the FORT_DG_SETUP subroutine.
      !        (See subrotine header for details.)
      !
      !      * keyword fort.dg file format rules are:
      !
      !          1) options are assigned in keyword = value format (e.g. fluxtype = 1)
      !          2) one option per line
      !          3) options can be specified in any order
      !          4) line beginning with ! indicates a comment
      !          5) blank lines are skipped
      !          6) comments may follow an assignment (e.g. fluxtype = 1 !comment)
      !          7) unrecognized keyword assignments are skipped      
      !          8) unassigned options that are not required will use 
      !             default values specified in FORT_DG_SETUP
      !          
      !  - Subroutines contained are:
      !
      !      1) READ_FIXED_FORT_DG
      !         
      !         * reads old fixed format fort.dg used in dgswem v11.13/dg-adcirc v22
      !
      !      2) READ_KEYWORD_FORT_DG
      !
      !         * reads keyword format fort.dg described above
      !
      !      3) CHECK_ERRORS
      !
      !         * Handles missing options 
      !         * Terminates if required options are missing
      !         * Warns that default values are used for missing optional options and continues      
      !
      !      4) FORT_DG_SETUP
      !
      !         * Responsible for configuring fort.dg options
      !         * MODIFICATIONS FOR ADDITION/REMOVAL OF FORT.DG OPTIONS SHOULD BE DONE HERE
      
      USE globals, ONLY: pres
      
      TYPE :: key_val
        CHARACTER(15) :: key            ! keyword
        REAL(pres), POINTER :: rptr       ! pointer to real target
        INTEGER, POINTER :: iptr        ! pointer to integer target
        CHARACTER(100), POINTER :: cptr ! pointer to character target      
        
        INTEGER :: vartype              ! target type indicator: 1=integer, 2=real, 3=character
        
        LOGICAL :: required             ! required/optional flag
        
        INTEGER :: flag                 ! successful read flag
      END TYPE key_val
      
      INTEGER, PARAMETER :: maxopt = 100          ! maximum allowable fort.dg options
      TYPE(key_val), DIMENSION(maxopt) :: dginp      
      
      INTEGER :: nopt                             ! number of valid options in dginp structure
      INTEGER, DIMENSION(maxopt) :: dginp_ind    ! indicies of valid options in dginp structure
 
      CONTAINS 
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
      SUBROUTINE read_input()
      
      USE messenger2, ONLY: abort,dirname,lname,myrank
      
      IMPLICIT NONE
      
      INTEGER :: opt_read
      INTEGER :: read_stat
      INTEGER :: skipped
      INTEGER :: read_file
      LOGICAL :: file_exists
      CHARACTER(100) :: temp
      
      ! initialize the dginp option structure
      CALL dginp_setup()      

      INQUIRE(FILE=dirname(1:lname)//'/'//'dgswe.inp', EXIST = file_exists)
      IF(file_exists == .FALSE.) THEN
        PRINT*, "dgswe.inp file does not exist"
        CALL abort()
      ENDIF      
      
      OPEN(UNIT=15,FILE=dirname(1:lname)//'/'//'dgswe.inp') 
      
      opt_read = 0
      skipped = 0 
      read_file = 0

      DO WHILE (opt_read < nopt)
      
        READ(15,"(A100)",IOSTAT=read_stat) temp
        IF(read_stat /= 0) THEN                    ! check for end-of-file
          EXIT
        ENDIF
        
        IF ( INDEX(temp,"!") == 1 .or. LEN(TRIM(temp)) == 0) THEN
            skipped = skipped + 1
        ELSE  
          read_file = 1
          EXIT
        ENDIF
        
      ENDDO
      
      CLOSE(15)
      

      IF (myrank == 0) THEN
        PRINT "(A)", "---------------------------------------------"
        PRINT "(A)", "             Input Information               "
        PRINT "(A)", "---------------------------------------------"
        PRINT "(A)", " "      

      
#ifdef rk22
        PRINT*, "RK22 timestepping"
#elif rk33
        PRINT*, "RK33 timestepping"
#else
        PRINT*, "Forward Euler timestepping"
#endif      
        PRINT*, " "      
      ENDIF        
      

      IF (read_file == 1) THEN
        IF (INDEX(temp,"=") > 0) THEN
          PRINT*, "Reading keyword format input file"
          CALL read_keyword_dginp()
        ELSE
          PRINT*, "Reading fixed format input file"
          CALL read_fixed_dginp()
        ENDIF
      ELSE 
        PRINT*, "ERROR: dgswe.inp does not contain any information"
        CALL abort()
      ENDIF
          
#ifdef openmp      
      IF (npart < nthreads) THEN
        npart = nthreads
      ENDIF  
#endif               
      
      
      END SUBROUTINE read_input
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      

      SUBROUTINE read_fixed_dginp()
      
      USE globals, ONLY: grid_file,forcing_file,p,ctp,hbp,dt,tf,dramp,cf,lines,out_direc,npart
      USE messenger2, ONLY: myrank,dirname,lname,finish,nthreads

      IMPLICIT NONE
      
      INTEGER, PARAMETER :: ninp = 12
      INTEGER :: i
      INTEGER :: exind
      INTEGER :: inp_read,skipped
      CHARACTER(100) :: temp
      CHARACTER(100) :: line
      LOGICAL :: file_exists
      
!       OPEN(UNIT=16, FILE=dirname(1:lname)//'/'//'fort.16')      
      
      OPEN(UNIT=15,FILE=dirname(1:lname)//'/'//'dgswe.inp',POSITION="rewind")      
      
      inp_read = 0
      skipped = 0
      DO WHILE (inp_read < ninp)
      
        READ(15,"(A100)") temp
                    
        IF ( INDEX(temp,"!") == 1 .or. LEN(TRIM(temp)) == 0) THEN
            skipped = skipped + 1
        ELSE
        
          exind = INDEX(temp,"!")
          IF (exind > 0) THEN
            line = ADJUSTL(temp(1:exind-1)) 
          ELSE
            line = ADJUSTL(temp)          
          ENDIF
          
          inp_read = inp_read + 1
          SELECT CASE (inp_read)
            CASE (1)
              grid_file = TRIM(line)
            CASE (2)
              forcing_file = TRIM(line)
            CASE (3)
              READ(line,*) p
            CASE (4)
              READ(line,*) ctp
            CASE (5)
              READ(line,*) hbp            
            CASE (6)
              READ(line,*) dt
            CASE (7)
              READ(line,*) tf
!               tf = tf*86400d0
            CASE (8)
              READ(line,*) dramp
            CASE (9)
              READ(line,*) cf
            CASE (10)
              READ(line,*) lines
            CASE (11)
              out_direc = TRIM(line)
            CASE (12)
              READ(line,*) npart
          END SELECT
            
        ENDIF
      
      
      ENDDO
           
      CLOSE(15)
      
 
      
      IF (myrank == 0) THEN     
      
      ! print inputs
      DO i = 1,maxopt        
        IF (ASSOCIATED(dginp(i)%iptr)) THEN  
          PRINT("(A,A,I8)"), dginp(i)%key," = ",dginp(i)%iptr
        ENDIF
        
        IF (ASSOCIATED(dginp(i)%rptr)) THEN 
          PRINT("(A,A,E21.8)"), dginp(i)%key," = ",dginp(i)%rptr
        ENDIF
        
        IF (ASSOCIATED(dginp(i)%cptr)) THEN 
          PRINT("(A,A,A)"), dginp(i)%key," = ",dginp(i)%cptr 
        ENDIF
      ENDDO      
      
      PRINT*, " "
      
      ENDIF

      RETURN
      END SUBROUTINE read_fixed_dginp
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 
      SUBROUTINE read_keyword_dginp()
      
      USE messenger2, ONLY: myrank,dirname,lname

      IMPLICIT NONE
      
      INTEGER :: i,j,opt
      INTEGER :: read_stat
      INTEGER :: opt_read
      INTEGER :: comment,blank
      INTEGER :: eqind,exind
      LOGICAL :: found      
      CHARACTER(100) :: temp,line,empty
      CHARACTER(15) :: test_opt
      CHARACTER(100) :: test_val     
      
      opt_read = 0
      comment = 0 
      blank = 0
      
      
      OPEN(25,FILE=dirname(1:lname)//'/'//'dgswe.inp',POSITION="rewind")        
     
      
      DO WHILE (opt_read < nopt)
      
        READ(25,"(A100)",IOSTAT=read_stat) temp
        IF(read_stat /= 0) THEN                    ! check for end-of-file
          EXIT
        ENDIF
        
        line = ADJUSTL(temp)
        
        IF(INDEX(line,"!") == 1) THEN              ! lines beginning with ! are skipped
        
          comment = comment + 1
          
        ELSE IF (LEN(TRIM(line)) == 0) THEN        ! blank lines are skipped
        
          blank = blank + 1
          
        ELSE
  
          ! determine keyword and assignment value
          eqind = INDEX(line,"=")
          exind = INDEX(line,"!")
          test_opt = line(1:eqind-1)
          IF (exind > 0) THEN                          ! handle trailing comment 
            test_val = ADJUSTL(line(eqind+1:exind-1))  ! (only necessary if there is no space between value and the !)
          ELSE
            test_val = ADJUSTL(line(eqind+1:))
          ENDIF         
          
          ! Look for a match for the keyword
          found = .false.
    test: DO opt = 1,nopt
    
            i = dginp_ind(opt)    
    
            IF (test_opt == dginp(i)%key) THEN
              
              ! Set variables equal to value from fort.dg through pointer using an internal read
              SELECT CASE (dginp(i)%vartype) 
                CASE(1)
                  READ(test_val,*) dginp(i)%iptr
                  IF (myrank == 0 ) PRINT("(A,A,I8)"), test_opt," = ",dginp(i)%iptr
                CASE(2)
                  READ(test_val,*) dginp(i)%rptr
                  IF (myrank == 0 ) PRINT("(A,A,E21.8)"), test_opt," = ",dginp(i)%rptr                  
                CASE(3)
                  dginp(i)%cptr = TRIM(test_val)
                  IF (myrank == 0 ) PRINT("(A,A,A)"), test_opt," = ",dginp(i)%cptr                  
              END SELECT

              found = .true.          ! flag match
              opt_read = opt_read + 1
              dginp(i)%flag = 1      ! flag option as found
              
              EXIT test
              
            ENDIF
          ENDDO test
                    
          IF (found == .false. .and. eqind > 0) THEN
            ! unmatched lines with an equal sign are either incorrect or no longer supported
            PRINT("(3A)"),"*** WARNING: ",test_opt, " is an incorrect or depreciated value ***"            
          ELSE IF (found == .false.) THEN
            ! unmatched lines without an equal sign are ignored
            PRINT("(A)"), "*** WARNING: non-comment line does not contain a keyword assignment***"           
          ENDIF
          
        ENDIF
      ENDDO 
      
      PRINT*, ""
     
      CALL check_errors(opt_read)
      
      PRINT*, ""
      CLOSE(25)
            
      END SUBROUTINE read_keyword_dginp
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     SUBROUTINE check_errors(opt_read)
     
     IMPLICIT NONE
     
     INTEGER :: i,j,opt
     INTEGER :: opt_read
     INTEGER :: quit
     
     IF(opt_read /= nopt) THEN

       ! check for required options that are unspecifed 
       quit = 0
       DO opt = 1,nopt
         i = dginp_ind(opt)
         IF (dginp(i)%flag == 0 .and. dginp(i)%required == .true.) THEN
           quit = 1   ! flag fatal error
         ENDIF
       ENDDO
        
       IF (quit == 1) THEN
        
          PRINT("(A)"), "*** ERROR: There are missing required options in the fort.dg file ***"  
          PRINT("(A)"), "           The following options must be specified: "      
          j = 0        
          DO opt = 1,nopt
            i = dginp_ind(opt)
            IF (dginp(i)%flag == 0 .and. dginp(i)%required == .true.) THEN
              j = j+1
              PRINT "(A,I3,2A)", "              ",j,") ",dginp(i)%key
            ENDIF
          ENDDO          
          
          PRINT("(A)"), "!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!"
          STOP
          
       ELSE
        
          PRINT("(A)"), "*** WARNING: There are missing optional options in the fort.dg file ***"
          PRINT("(A)"), "             The following default values will be used: "    
          j = 0        
          DO opt = 1,nopt
            i = dginp_ind(opt)
            IF (dginp(i)%flag == 0 .and. dginp(i)%required == .false.) THEN
              
              j = j+1
              SELECT CASE (dginp(i)%vartype) 
                CASE(1)
                  PRINT("(A,I3,A,A,A,I8)"),     "              ",j,") ",dginp(i)%key," = ",dginp(i)%iptr
                CASE(2)
                  PRINT("(A,I3,A,A,A,E21.8)"),  "              ",j,") ",dginp(i)%key," = ",dginp(i)%rptr                  
                CASE(3)
                  PRINT("(A,I3,A,A,A,A)"),      "              ",j,") ",dginp(i)%key," = ",dginp(i)%cptr                  
              END SELECT
              
            ENDIF
          ENDDO 
          
          PRINT("(A)"), '!!!!!! EXECUTION WILL CONTINUE !!!!!!!!'
          
       ENDIF       
                  
     ENDIF        
     
     
     RETURN
     END SUBROUTINE check_errors
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!           
      
      SUBROUTINE dginp_setup()
      
      ! Subroutine that configures the fort.dg options
      !
      !   This subroutine is meant to add flexibility in adding/depreciating  
      !   features while maintaining forward (and some degree of backward) compatibility 
      !
      !   - Options can be added to the fort.dg file by:
      !       1) Specifying a keyword in a unused index (<= maxopt) of the dginp structure
      !       2) Associating the appropriate pointer with the corresponding variable
      !          Note: pointer must agree with the associated variable type 
      !                (iptr=integer, rptr=real, cptr=character)      
      !          Note: the associated variable must be declared using the TARGET attribute
      !       3) Specifying whether the variable is required 
      !       4) Providing a default value
      ! 
      !   - Options can be removed from the fort.dg file by:
      !       1) Commenting out or deleting an existing entry in the dginp structure
      !          Note: re-indexing subsequent entries is not necessary (see dginp(17) below)
      !       
      !       OR
      !
      !       2) Setting the dginp(i)%required variable to .false.
      ! 
      !   - New features should be added as dginp(i)%required = .false. as much as possible 
      !     to maintain backward compatibility, older fort.dg files not containing these 
      !     options will cause provided default values to be used (these should be set so 
      !     the feature is turned off)
      !
      !   - fort.dg files containing new feature options can still be used for previous  
      !     versions of the code because the new options will be ignored
      
      
      USE globals, ONLY: grid_file,forcing_file,p,ctp,hbp,dt,tf,dramp,cf,lines,out_direc,npart
      
      IMPLICIT NONE        
      
      INTEGER :: i
      INTEGER :: ncheck
      CHARACTER(15) :: empty
      
      ! initialize dginp structure
      DO i = 1,maxopt
        NULLIFY(dginp(i)%iptr)
        NULLIFY(dginp(i)%rptr)
        NULLIFY(dginp(i)%cptr)
        
        dginp(i)%key = empty
        dginp(i)%flag = 0        
      ENDDO
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Configure fort.dg options here:
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      !    keywords                         target variables                      requirement                 default values
      dginp(1)%key = "grid_file";       dginp(1)%cptr => grid_file ;      dginp(1)%required = .true.;      dginp(1)%cptr = ""
      dginp(2)%key = "forcing_file";    dginp(2)%cptr => forcing_file;    dginp(2)%required = .true.;      dginp(2)%cptr = ""
      dginp(3)%key = "p";               dginp(3)%iptr => p;               dginp(3)%required = .true.;      dginp(3)%iptr = 1
      dginp(4)%key = "ctp";             dginp(4)%iptr => ctp;             dginp(4)%required = .false.;     dginp(4)%iptr = 1
      dginp(5)%key = "hbp";             dginp(5)%iptr => hbp;             dginp(5)%required = .false.;     dginp(5)%iptr = 1
      dginp(6)%key = "dt";              dginp(6)%rptr => dt;              dginp(6)%required = .true.;      dginp(6)%rptr = 0.5d0
      dginp(7)%key = "tf";              dginp(7)%rptr => tf;              dginp(7)%required = .true.;      dginp(7)%rptr = 1
      dginp(8)%key = "dramp";           dginp(8)%rptr => dramp;           dginp(8)%required = .true.;      dginp(8)%rptr = 0.08d0
      dginp(9)%key = "cf";              dginp(9)%rptr => cf;              dginp(9)%required = .true.;      dginp(9)%rptr = 0.0025d0
      dginp(10)%key = "lines";          dginp(10)%rptr => lines;          dginp(10)%required = .true.;     dginp(10)%rptr = 10d0
      dginp(11)%key = "out_direc";      dginp(11)%cptr => out_direc;      dginp(11)%required = .false.;    dginp(11)%cptr = "./"
      dginp(12)%key = "npart";          dginp(12)%iptr => npart;          dginp(12)%required = .true.;     dginp(12)%iptr = 1

      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! End configuration
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      nopt = 0
      ncheck = 0
      DO i = 1,maxopt
      
        ! find and keep track of populated indicies
        IF (dginp(i)%key .ne. empty) THEN      
          nopt = nopt + 1      
          dginp_ind(nopt) = i
        ENDIF
        
        ! determine target variable type by checking association status
        dginp(i)%vartype = 0    
        
        IF (ASSOCIATED(dginp(i)%iptr)) THEN  ! integer
          ncheck = ncheck + 1   
          dginp(i)%vartype = 1
        ENDIF
        
        IF (ASSOCIATED(dginp(i)%rptr)) THEN ! real
          ncheck = ncheck + 1
          dginp(i)%vartype = 2
        ENDIF
        
        IF (ASSOCIATED(dginp(i)%cptr)) THEN ! character
          ncheck = ncheck + 1        
          dginp(i)%vartype = 3        
        ENDIF
      ENDDO
      
!       PRINT*, "Number of options = ", nopt
!       PRINT*, "Number of pointer associations = ", ncheck
      
      ! ensure user has associated each keyword pointer
      IF (nopt /= ncheck) THEN
        PRINT("(A)"), "*** ERROR: fort.dg option pointer association error ***"
        PRINT("(A)"), "           check keyword configuration in dginp_setup subroutine"
        STOP
      ENDIF
          
      
      RETURN
      END SUBROUTINE dginp_setup
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!            
      
      END MODULE read_dginp