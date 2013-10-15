      SUBROUTINE read_forcing()

      USE globals, ONLY: pres,forcing_file,neta,nvel,nbou,fbseg, &
                         nobfr,obtag,obtag2,obfreq,obnfact,obeq,obamp,obph, &
                         nfbfr,fbtag,fbtag2,fbfreq,fbnfact,fbeq,fbamp,fbph, &
                         fbper,obper,pi

      IMPLICIT NONE
      INTEGER :: alloc_stat
      INTEGER :: bfr,node,seg,segtype
      REAL(pres) :: deg2rad

      deg2rad = pi/180d0
      
      OPEN(unit=15,file=forcing_file)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Read in open boundary forcing data
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      READ(15,*) nobfr
!       PRINT*, 'number of oben boundary forcings',nobfr

      ALLOCATE(obtag(nobfr),obfreq(nobfr),obnfact(nobfr),obeq(nobfr),obper(nobfr),STAT = alloc_stat)
      IF(alloc_stat /= 0) THEN
        PRINT*, "Allocation error: obtag, obfreq, obnfact, obeq"
      ENDIF

      DO bfr = 1,nobfr
        READ(15,*) obtag(bfr)
        READ(15,*) obfreq(bfr),obnfact(bfr),obeq(bfr)
        obeq(bfr) = obeq(bfr)*deg2rad
        IF(obfreq(bfr) == 0.) THEN
          obper(bfr) = 1d0
        ELSE
          obper(bfr) = 2d0*pi/obfreq(bfr)
        ENDIF
      ENDDO

      ALLOCATE(obtag2(nobfr),obamp(neta,nobfr),obph(neta,nobfr),STAT = alloc_stat)
      IF(alloc_stat /= 0) THEN
        PRINT*, "Allocation error: obtag2, obamp, obph"
      ENDIF

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
!       PRINT*, 'number of flow boundary forcings',nfbfr
!       PRINT*, ' '

      ALLOCATE(fbtag(nfbfr),fbfreq(nfbfr),fbnfact(nfbfr),fbeq(nfbfr),fbper(nfbfr),STAT = alloc_stat)
      IF(alloc_stat /= 0) THEN
        PRINT*, "Allocation error: fbtag, fbfreq, fbnfact, fbeq"
      ENDIF

      DO bfr = 1,nfbfr
        READ(15,*) fbtag(bfr)
        READ(15,*) fbfreq(bfr),fbnfact(bfr),fbeq(bfr)
!         PRINT*,fbtag(bfr)
!         PRINT*, fbfreq(bfr),fbnfact(bfr),fbeq(bfr)
        fbeq(bfr) = fbeq(bfr)*deg2rad
        IF(fbfreq(bfr) == 0.) THEN
          fbper(bfr) = 1d0
        ELSE
          fbper(bfr) = 2d0*pi/fbfreq(bfr)        
        ENDIF
      ENDDO

      ALLOCATE(fbtag2(nfbfr),fbamp(nvel,nfbfr),fbph(nvel,nfbfr),STAT = alloc_stat)
      IF(alloc_stat /= 0) THEN
        PRINT*, "Allocation error: fbtag2, fbamp, fbph"
      ENDIF

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