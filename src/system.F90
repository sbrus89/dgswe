      SUBROUTINE EXECUTE_COMMAND_LINE(command)

      IMPLICIT NONE
      
      CHARACTER(*) :: command
      
      CALL SYSTEM(command)

      RETURN
      END SUBROUTINE EXECUTE_COMMAND_LINE
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      FUNCTION sha1(file_name,out_dir)
      
      IMPLICIT NONE
      
      CHARACTER(40) :: sha1
      CHARACTER(*)  :: file_name
      CHARACTER(*)  :: out_dir
      LOGICAL :: file_exists
      
      INQUIRE(FILE=file_name,EXIST=file_exists)
      IF (file_exists) THEN
        CALL EXECUTE_COMMAND_LINE("sha1sum "//TRIM(ADJUSTL(file_name))//" > "//TRIM(ADJUSTL(out_dir))//"sha1.out")    
        OPEN(UNIT=9, FILE=TRIM(ADJUSTL(out_dir))//"sha1.out")
        READ(9,"(A)") sha1
        CLOSE(9, STATUS='DELETE')              
!        PRINT*, sha1
      ELSE
        sha1 = "                                                     "
      ENDIF
      
      RETURN
      END FUNCTION sha1
      