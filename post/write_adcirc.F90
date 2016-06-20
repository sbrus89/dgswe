      MODULE write_adcirc
      
      USE globals, ONLY:rp
      
      INTEGER :: unit_counter = 10
      
      CONTAINS
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
      
      SUBROUTINE write_6364_header(filename,rundes,runid,agrid,ndsetse,np,dtdp,nspoolge,irtype,file_unit)
      
      IMPLICIT NONE
      
      CHARACTER(*), INTENT(IN) :: filename
      CHARACTER(32), INTENT(IN) :: rundes
      CHARACTER(24), INTENT(IN)  :: runid
      CHARACTER(24), INTENT(IN)  :: agrid
      INTEGER, INTENT(IN) :: ndsetse
      INTEGER, INTENT(IN) :: np
      REAL(rp), INTENT(IN)  :: dtdp
      INTEGER, INTENT(IN)  :: nspoolge
      INTEGER, INTENT(IN) :: irtype
      INTEGER, INTENT(OUT) :: file_unit
      
      file_unit = unit_counter
      OPEN(UNIT=file_unit,FILE=filename)
      unit_counter = unit_counter + 1
      
      WRITE(file_unit,"(A32,2x,A24,2x,A24)") rundes,runid,agrid
      WRITE(file_unit,*) ndsetse,np,dtdp*nspoolge,nspoolge,irtype
      
      RETURN
      END SUBROUTINE write_6364_header

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
      
      SUBROUTINE write_6364_snap(file_unit,np,n,time,it,nodal_sol)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: file_unit
      INTEGER, INTENT(IN) :: np
      INTEGER, INTENT(IN) :: n
      REAL(rp), INTENT(IN) :: time
      INTEGER, INTENT(IN) :: it
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: nodal_sol
      
      INTEGER :: k,j
      
      WRITE(file_unit,*) time,it
      DO k = 1,np
        WRITE(file_unit,*) k,(nodal_sol(k,j), j = 1,n)
      ENDDO
      
      RETURN
      END SUBROUTINE write_6364_snap
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       

      END MODULE write_adcirc