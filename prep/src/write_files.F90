      SUBROUTINE write_files()
      
      USE globals, ONLY: nproc,ne,nresel,nresnd,nd_l2g,xy,depth,lect,nope,nbou,fbseg, &
                         lnope,lneta,lobseg,lobnds,lnbou,lnvel,lfbseg,lfbnds, &
                         nobfr,obtag,obtag2,obfreq,obnfact,obeq,lobamp,lobph, &
                         nfbfr,fbtag,fbtag2,fbfreq,fbnfact,fbeq,lfbamp,lfbph,lnbouf, &
                         nsred,ned_sr,pe_sr,el_sr,led_sr,el_l2g, &
                         grid_file,forcing_file,p,ctp,dt,tf,dramp,cf,lines,nblk,npart,mndof

      IMPLICIT NONE

      
      INTEGER :: pe,pes,nd,el,bnd,bfr,seg,ed
      INTEGER :: segtype
      INTEGER :: gnd
      INTEGER :: npes,ned2pes
      INTEGER :: match
      INTEGER :: lname
      CHARACTER(6) :: dirname
      LOGICAL :: file_exists

      INTEGER, ALLOCATABLE, DIMENSION(:) :: sended

      lname = 6
      dirname = "PE0000"
!      CALL SYSTEM("rm -r PE*")
      
      
      ! Create the local PE directories 
      PRINT*, "Creating local PE directories", nproc
      DO pe = 1,nproc
        WRITE(dirname(3:lname),"(I4.4)") pe-1      
        CALL SYSTEM("mkdir -v "//dirname)      
!        CALL ISHELL("mkdir -p "//dirname)
!        CALL EXECUTE_COMMAND_LINE("mkdir -p "//dirname)
      ENDDO
      PRINT*, "Done creating directories"
      
      PRINT*, "PE    ","NRESEL    ","NRESND    ","NED_SR    "
      
      ! Write the fort.14 grid file
      DO pe = 1,nproc
        
        WRITE(dirname(3:lname),"(I4.4)") pe-1      
        OPEN(UNIT=14,FILE=dirname(1:lname)//'/'//'fort.14')

        WRITE(14,"(A)") dirname(1:lname)
        WRITE(14,*) nresel(pe),nresnd(pe)
        PRINT*, pe, nresel(pe),nresnd(pe),ned_sr(pe)
        
        DO nd = 1,nresnd(pe)
          gnd = nd_l2g(nd,pe)
          WRITE(14,*) nd, xy(1,gnd), xy(2,gnd), depth(gnd)
        ENDDO
        
        DO el = 1,nresel(pe)
          WRITE(14,*) el, 3, lect(1,el,pe), lect(2,el,pe), lect(3,el,pe)
        ENDDO
        
        WRITE(14,*) lnope(pe)
        WRITE(14,*) lneta(pe)
        
        DO bnd = 1,nope
          IF(lobseg(bnd,pe) > 0) THEN
            WRITE(14,*) lobseg(bnd,pe)
            DO nd = 1,lobseg(bnd,pe)
              WRITE(14,*) lobnds(nd,bnd,pe)
            ENDDO
          ENDIF
        ENDDO
        
        WRITE(14,*) lnbou(pe)
        WRITE(14,*) lnvel(pe)
        
        DO bnd = 1,nbou
          IF(lfbseg(bnd,pe) > 0) THEN
            WRITE(14,*) lfbseg(bnd,pe),fbseg(2,bnd)
            DO nd = 1,lfbseg(bnd,pe) 
              WRITE(14,*) lfbnds(nd,bnd,pe)
            ENDDO
          ENDIF
        ENDDO
        
        CLOSE(14)
         
      ENDDO
      
      ! Write the fort.15 boundary forcing file
      DO pe = 1,nproc
        WRITE(dirname(3:lname),"(I4.4)") pe-1      
        OPEN(UNIT=15,FILE=dirname(1:lname)//'/'//'fort.15')      
        
        IF(lnope(pe) > 0) THEN
          WRITE(15,*) nobfr
          DO bfr = 1,nobfr
            WRITE(15,*) obtag(bfr)
            WRITE(15,*) obfreq(bfr),obnfact(bfr),obeq(bfr)
          ENDDO
          DO bfr = 1,nobfr
            WRITE(15,*) obtag2(bfr)
            DO nd = 1,lneta(pe)
              WRITE(15,*) lobamp(nd,bfr,pe),lobph(nd,bfr,pe)
            ENDDO
          ENDDO
        ELSE
          WRITE(15,*) 0
        ENDIF
        
        IF(lnbouf(pe) > 0) THEN
          WRITE(15,*) nfbfr
          DO bfr = 1,nfbfr
            WRITE(15,*) fbtag(bfr)
            WRITE(15,*) fbfreq(bfr),fbnfact(bfr),fbeq(bfr)
          ENDDO
          DO bfr = 1,nfbfr
            WRITE(15,*) fbtag2(bfr)
            DO seg = 1,nbou
              segtype = fbseg(2,seg)
              IF(segtype == 2 .OR. segtype == 12 .OR. segtype == 22)THEN
                DO nd = 1,lfbseg(seg,pe)
                  WRITE(15,*) lfbamp(nd,seg,bfr,pe), lfbph(nd,seg,bfr,pe)
                ENDDO
              ENDIF
            ENDDO
          ENDDO
        ENDIF
           
        CLOSE(15)
      ENDDO
      
      ! Write the fort.18 message passing file
      ALLOCATE(sended(nsred))
      DO pe = 1,nproc
        WRITE(dirname(3:lname),"(I4.4)") pe-1      
        OPEN(UNIT=18,FILE=dirname(1:lname)//'/'//'fort.18')  
         
        npes = 0
        DO pes = 1,nproc   
          IF(pe /= pes) THEN
            match = 0
            DO ed = 1,ned_sr(pe)
              IF(pe_sr(ed,pe) == pes) THEN
                match = 1
              ENDIF
            ENDDO
            npes = npes + match
          ENDIF
        ENDDO

        WRITE(18,182) npes

        DO pes = 1,nproc
          IF(pe /= pes) THEN
            ned2pes = 0
            DO ed = 1,ned_sr(pe)
              IF(pe_sr(ed,pe) == pes) THEN
                ned2pes = ned2pes + 1
                sended(ned2pes) = ed
              ENDIF
            ENDDO
            IF(ned2pes > 0) THEN
              WRITE(18,181) pes-1,ned2pes
              WRITE(18,180) (el_sr(sended(ed),pe), ed = 1,ned2pes)
              WRITE(18,180) (led_sr(sended(ed),pe), ed = 1,ned2pes)
            ENDIF
          ENDIF
        ENDDO
        
        
        CLOSE(18)
      ENDDO
      
      
      ! Write the local to global element table
      DO pe = 1,nproc
        WRITE(dirname(3:lname),"(I4.4)") pe-1      
        OPEN(UNIT=80,FILE=dirname(1:lname)//'/'//'fort.80') 
        
        WRITE(80,*) nproc
        WRITE(80,*) ne
        WRITE(80,*) MAXVAL(nresel)
        WRITE(80,*) nresel(pe)
        WRITE(80,*) mndof
        WRITE(80,*) lines
        
        DO el = 1,nresel(pe)
          WRITE(80,*) el_l2g(el,pe) 
        ENDDO
        
      ENDDO  
      
      ! Write the local input file
      DO pe = 1,nproc
        WRITE(dirname(3:lname),"(I4.4)") pe-1         
        OPEN(UNIT=10,FILE=dirname(1:lname)//'/'//'dgswe.inp')
        
        WRITE(10,"(A)") dirname(1:lname)//'/'//"fort.14"
        WRITE(10,"(A)") dirname(1:lname)//'/'//"fort.15"
        WRITE(10,"(I5)") p
        WRITE(10,"(I5)") ctp
        WRITE(10,"(E24.17)") dt
        WRITE(10,"(E24.17)") tf
        WRITE(10,"(E24.17)") dramp
        WRITE(10,"(E24.17)") cf
        WRITE(10,"(I7)") int(lines)
        WRITE(10,"(A)") './' // dirname(1:lname) // '/'
        WRITE(10,"(I5)") npart
        
        CLOSE(10)
      ENDDO

      
 180  FORMAT(8X,9I8)
 181  FORMAT(8X,2I8) 
 182  FORMAT(8X,I8)
 
      RETURN
      END SUBROUTINE write_files
