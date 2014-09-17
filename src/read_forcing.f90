      SUBROUTINE read_forcing()

      USE globals, ONLY: pres,forcing_file,neta,nvel,nbou,fbseg, &
                         nobfr,obtag,obtag2,obfreq,obnfact,obeq,obamp,obph, &
                         nfbfr,fbtag,fbtag2,fbfreq,fbnfact,fbeq,fbamp,fbph, &
                         fbper,obper,pi
                         
      USE allocation, ONLY: alloc_forcing_arrays                         

      IMPLICIT NONE
      INTEGER :: bfr,node,seg,segtype
      REAL(pres) :: deg2rad

      deg2rad = pi/180d0
      
      OPEN(unit=15,file=forcing_file)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Read in open boundary forcing data
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      READ(15,*) nobfr
      PRINT "(A,I5)", 'Number of open boundary forcings',nobfr

      CALL alloc_forcing_arrays(1)

      DO bfr = 1,nobfr
        READ(15,*) obtag(bfr)
        PRINT "(A,A)", "  ",obtag(bfr)
        READ(15,*) obfreq(bfr),obnfact(bfr),obeq(bfr)
        obeq(bfr) = obeq(bfr)*deg2rad
        IF(obfreq(bfr) == 0.) THEN
          obper(bfr) = 1d0
        ELSE
          obper(bfr) = 2d0*pi/obfreq(bfr)
        ENDIF
      ENDDO

      DO bfr = 1,nobfr
        READ(15,*) obtag2(bfr) 
        DO node = 1,neta
          READ(15,*) obamp(node,bfr),obph(node,bfr)
          obph(node,bfr) = obph(node,bfr)*deg2rad
        ENDDO
      ENDDO

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Read in flow boundary forcing data
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      READ(15,*) nfbfr
      PRINT "(A,I5)", 'Number of flow boundary forcings',nfbfr
      PRINT*, ' '

      CALL alloc_forcing_arrays(2)

      DO bfr = 1,nfbfr
        READ(15,*) fbtag(bfr)
        PRINT "(A,A)", "  ",obtag(bfr)        
        READ(15,*) fbfreq(bfr),fbnfact(bfr),fbeq(bfr)
!         PRINT*, fbfreq(bfr),fbnfact(bfr),fbeq(bfr)
        fbeq(bfr) = fbeq(bfr)*deg2rad
        IF(fbfreq(bfr) == 0.) THEN
          fbper(bfr) = 1d0
        ELSE
          fbper(bfr) = 2d0*pi/fbfreq(bfr)        
        ENDIF
      ENDDO


      DO bfr = 1,nfbfr
        READ(15,*) fbtag2(bfr) 
        DO seg = 1,nbou
            segtype = fbseg(2,seg)
            IF(segtype == 2 .OR. segtype == 12 .OR. segtype == 22)THEN
            DO node = 1,fbseg(1,seg)
              READ(15,*) fbamp(node,bfr),fbph(node,bfr)
              fbph(node,bfr) = fbph(node,bfr)*deg2rad
            ENDDO
          ENDIF
        ENDDO
      ENDDO


      RETURN
      END SUBROUTINE read_forcing