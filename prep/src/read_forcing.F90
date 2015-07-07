      SUBROUTINE read_forcing()

      USE globals, ONLY: pres,neta,nvel,nbou,fbseg, &
                         nobfr,obtag,obtag2,obfreq,obnfact,obeq,obamp,obph, &
                         nfbfr,fbtag,fbtag2,fbfreq,fbnfact,fbeq,fbamp,fbph
                         
      USE allocation, ONLY: alloc_forcing_arrays    
      USE read_dginp, ONLY: forcing_file

      IMPLICIT NONE
      INTEGER :: alloc_stat
      INTEGER :: bfr,node,seg,segtype
      LOGICAL :: file_exists, any_nfb
      
      file_exists = .FALSE.
      
!       DO 
!         PRINT*, "Input boundary forcing file"
!         READ(*,"(A)") forcing_file
!         PRINT*, " "    
!         
!         INQUIRE(FILE=forcing_file,EXIST=file_exists)
!         IF(file_exists) THEN
!           EXIT
!         ELSE  
!           PRINT*, "file not found, try again"
!         ENDIF
!       ENDDO         

      
      INQUIRE(FILE=forcing_file,EXIST=file_exists)
      IF(file_exists) THEN
        OPEN(unit=15,file=forcing_file)          
      ELSE  
        PRINT*, "file not found, try again"
        STOP
      ENDIF

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Read in open boundary forcing data
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      READ(15,*) nobfr
!       PRINT*, 'number of oben boundary forcings',nobfr

      CALL alloc_forcing_arrays(1)

      DO bfr = 1,nobfr
        READ(15,*) obtag(bfr)
        READ(15,*) obfreq(bfr),obnfact(bfr),obeq(bfr)
      ENDDO

      DO bfr = 1,nobfr
        READ(15,*) obtag2(bfr) 
        DO node = 1,neta
          READ(15,*) obamp(node,bfr),obph(node,bfr)
        ENDDO
      ENDDO

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Read in flow boundary forcing data
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      any_nfb = .false. ! determine if there are normal flow boundaries
      DO seg = 1,nbou
        segtype = fbseg(2,seg)
        IF(segtype == 2 .OR. segtype == 12 .OR. segtype == 22)THEN
          any_nfb = .true.
        ENDIF        
      ENDDO
      
      nfbfr = 0
      IF (any_nfb) then
        READ(15,*) nfbfr

        CALL alloc_forcing_arrays(2)

        DO bfr = 1,nfbfr
          READ(15,*) fbtag(bfr)
          READ(15,*) fbfreq(bfr),fbnfact(bfr),fbeq(bfr)
        ENDDO

        DO bfr = 1,nfbfr
          READ(15,*) fbtag2(bfr) 
          DO seg = 1,nbou
              segtype = fbseg(2,seg)
              IF(segtype == 2 .OR. segtype == 12 .OR. segtype == 22)THEN
              DO node = 1,fbseg(1,seg)
                READ(15,*) fbamp(node,seg,bfr),fbph(node,seg,bfr)
              ENDDO
            ENDIF
          ENDDO
        ENDDO
        
      ENDIF
      CLOSE(15)
      
      RETURN
      END SUBROUTINE read_forcing
