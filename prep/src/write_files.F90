      SUBROUTINE write_files()
      
      USE globals, ONLY: ne,nn,nresel,nresnd,nd_l2g,xy,depth,lect,lnelnds,nope,nbou,fbseg, &
                         lnope,lneta,lobseg,lobnds,lnbou,lnvel,lfbseg,lfbnds,lbndxy, &
                         nobfr,obtag,obtag2,obfreq,obnfact,obeq,lobamp,lobph, &
                         nfbfr,fbtag,fbtag2,fbfreq,fbnfact,fbeq,lfbamp,lfbph,lnbouf, &
                         nfbsfr,fbstag,fbstag2,fbsbgn,fbsend,fbssig,lfbsamp, &
                         nsred,el_l2g,el_g2l,mndof, &
                         el_type,elhb,nnds,order, &
                         nsta,nlsta,sta_l2g,xysta
                         
      USE read_dginp, ONLY: write_local,hbp,ctp,curve_file,bathy_file,cb_file_exists,hb_file_exists, &
                            sol_opt,sta_opt,sol_snap,sta_snap,tf,dt
      USE messenger2, ONLY: nproc,nx_sr,ny_sr,nqpte_sr,hb_sr,detJe_sr, &
                            ned_sr,pe_sr,el_sr,led_sr
                            
      USE output, ONLY: time_snaps                            

      IMPLICIT NONE

      
      INTEGER :: i,k,pe,pes,nd,el,bnd,bfr,seg,ed,pt,sed,sta
      INTEGER :: n,segtype,et
      INTEGER :: gnd,gel,lnd,lsta,gsta
      INTEGER :: nnd
      INTEGER :: npes,ned2pes
      INTEGER :: match
      INTEGER :: lname
      INTEGER :: tskp_sol,tskp_sta
      INTEGER :: nout_sol,nout_sta
      CHARACTER(6) :: dirname
      LOGICAL :: file_exists

      INTEGER, ALLOCATABLE, DIMENSION(:) :: sended
      
      PRINT "(A)", ""
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", "           Subdomain Information             "
      PRINT "(A)", "---------------------------------------------"      

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
 grd: DO pe = 1,nproc
        
        WRITE(dirname(3:lname),"(I4.4)") pe-1      
        OPEN(UNIT=14,FILE=dirname(1:lname)//'/'//'fort.14')

        WRITE(14,"(A)") dirname(1:lname)
        WRITE(14,"(2(I8,1x))") nresel(pe),nresnd(pe)
        
        PRINT*, pe, nresel(pe),nresnd(pe),ned_sr(pe)
        
        DO nd = 1,nresnd(pe)
          gnd = nd_l2g(nd,pe)
          WRITE(14,"(I8,1X,3(E24.17,1X))") nd, xy(1,gnd), xy(2,gnd), depth(gnd)
        ENDDO
        
        DO el = 1,nresel(pe)
          WRITE(14,"(15(I8,1x))") el, lnelnds(el,pe), (lect(nd,el,pe), nd = 1,lnelnds(el,pe))
        ENDDO
        
        WRITE(14,"(I8,19x,A)") lnope(pe), "! number of open boundaries"
        WRITE(14,"(I8,19x,A)") lneta(pe), "! number of total elevation specified boundary nodes"
        
        i = 0
        DO bnd = 1,nope
          IF(lobseg(bnd,pe) > 0) THEN
            i = i + 1
            WRITE(14,"(I8,10x,A,1x,I8)") lobseg(bnd,pe), "! number of nodes in open boundary", i
            DO nd = 1,lobseg(bnd,pe)
              WRITE(14,"(I8)") lobnds(nd,bnd,pe)
            ENDDO
          ENDIF
        ENDDO
        
        WRITE(14,"(I8,19x,A)") lnbou(pe), "! number of normal flow boundaries"
        WRITE(14,"(I8,19x,A)") lnvel(pe), "! total number of normal flow nodes"
        
        i = 0
        DO bnd = 1,lnbou(pe)
          IF(lfbseg(1,bnd,pe) > 0) THEN
            i = i + 1
            WRITE(14,"(I8,1x,I8,10x,A,1x,I8)") lfbseg(1,bnd,pe),lfbseg(2,bnd,pe), "! number of nodes in normal flow boundary", i
            DO nd = 1,lfbseg(1,bnd,pe) 
              WRITE(14,"(I8)") lfbnds(nd,bnd,pe)
            ENDDO
          ENDIF
        ENDDO
        
        CLOSE(14)
        
        OPEN(UNIT=141,FILE=dirname(1:lname)//'/'//'bnd.14')   
        
        WRITE(141,*) lnope(pe)
        WRITE(141,*) lneta(pe)
        
        DO bnd = 1,nope
          IF(lobseg(bnd,pe) > 0) THEN
            WRITE(141,*) lobseg(bnd,pe)
            DO nd = 1,lobseg(bnd,pe)
              WRITE(141,*) nd_l2g(lobnds(nd,bnd,pe),pe)
            ENDDO
          ENDIF
        ENDDO
        
        WRITE(141,*) lnbou(pe)
        WRITE(141,*) lnvel(pe)
        
        DO bnd = 1,lnbou(pe)
          IF(lfbseg(1,bnd,pe) > 0) THEN
            WRITE(141,*) lfbseg(1,bnd,pe),lfbseg(2,bnd,pe)
            DO nd = 1,lfbseg(1,bnd,pe) 
              WRITE(141,*) nd_l2g(lfbnds(nd,bnd,pe),pe)
            ENDDO
          ENDIF
        ENDDO        
        
        CLOSE(141)
         
      ENDDO grd
      
      ! Write the fort.15 boundary forcing file
      DO pe = 1,nproc
        WRITE(dirname(3:lname),"(I4.4)") pe-1      
        OPEN(UNIT=15,FILE=dirname(1:lname)//'/'//'fort.15')      
        
        IF(lnope(pe) > 0) THEN
          WRITE(15,*) nobfr
          DO bfr = 1,nobfr
            WRITE(15,*) obtag(bfr)
            WRITE(15,"(3(D24.17,1x))") obfreq(bfr),obnfact(bfr),obeq(bfr)
          ENDDO
          DO bfr = 1,nobfr
            WRITE(15,*) obtag2(bfr)
            DO seg = 1,lnope(pe)
              DO nd = 1,lobseg(seg,pe)
                WRITE(15,"(2(D24.17,1x))") lobamp(nd,seg,bfr,pe),lobph(nd,seg,bfr,pe)
              ENDDO
            ENDDO
          ENDDO
        ELSE
          WRITE(15,*) 0
        ENDIF
        
        IF(lnbouf(pe) > 0) THEN
          WRITE(15,*) nfbfr
          DO bfr = 1,nfbfr
            WRITE(15,*) fbtag(bfr)
            WRITE(15,"(3(D24.17,1x))") fbfreq(bfr),fbnfact(bfr),fbeq(bfr)
          ENDDO
          DO bfr = 1,nfbfr
            WRITE(15,*) fbtag2(bfr)
            DO seg = 1,lnbou(pe)
              segtype = lfbseg(2,seg,pe)
              IF(segtype == 2 .OR. segtype == 12 .OR. segtype == 22)THEN
                DO nd = 1,lfbseg(1,seg,pe)
                  WRITE(15,"(2(D24.17,1x))") lfbamp(nd,seg,bfr,pe), lfbph(nd,seg,bfr,pe)
                ENDDO
              ENDIF
            ENDDO
          ENDDO
          
          WRITE(15,*) nfbsfr
          DO bfr = 1,nfbsfr
            WRITE(15,*) fbstag(bfr)
            WRITE(15,"(3(D24.17,1x))") fbsbgn(bfr),fbsend(bfr),fbssig(bfr)
          ENDDO
          DO bfr = 1,nfbsfr
            WRITE(15,*) fbstag2(bfr)
            DO seg = 1,lnbou(pe)
              segtype = lfbseg(2,seg,pe)
              IF(segtype == 2 .OR. segtype == 12 .OR. segtype == 22)THEN
                DO nd = 1,lfbseg(1,seg,pe)
                  WRITE(15,"(2(D24.17,1x))") lfbsamp(nd,seg,bfr,pe)
                ENDDO
              ENDIF
            ENDDO
          ENDDO          
        ENDIF
           
        CLOSE(15)
      ENDDO
      
      ! Write the local stations file
      DO pe = 1,nproc
        WRITE(dirname(3:lname),"(I4.4)") pe-1      
        OPEN(UNIT=18,FILE=dirname(1:lname)//'/'//'fort.sta') 
        
        WRITE(18,*) nlsta(pe)
        DO lsta = 1,nlsta(pe)
          gsta = sta_l2g(lsta,pe)
          WRITE(18,"(2(e24.17,1x))") xysta(1,gsta), xysta(2,gsta)
        ENDDO
        
        CLOSE(18)
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
              DO ed = 1,ned2pes
                sed = sended(ed)
                WRITE(18,182) nqpte_sr(sed,pe)
                DO pt = 1,nqpte_sr(sed,pe)
                  WRITE(18,183) nx_sr(pt,sed,pe), ny_sr(pt,sed,pe), hb_sr(pt,sed,pe), detJe_sr(pt,sed,pe)
                ENDDO
              ENDDO  
            ENDIF
          ENDIF
        ENDDO
        
        
        CLOSE(18)
      ENDDO 
      
      
      
      
      CALL time_snaps(sol_opt,sol_snap,tf,dt,tskp_sol,nout_sol)
      
      ! Write the local to global element table
      DO pe = 1,nproc
        WRITE(dirname(3:lname),"(I4.4)") pe-1      
        OPEN(UNIT=80,FILE=dirname(1:lname)//'/'//'fort.80') 
        
        WRITE(80,*) nproc
        WRITE(80,*) ne
        WRITE(80,*) MAXVAL(nresel)
        WRITE(80,*) nresel(pe)
        WRITE(80,*) mndof
        WRITE(80,*) nout_sol
        
        DO el = 1,nresel(pe)
          WRITE(80,*) el, el_l2g(el,pe) 
        ENDDO
        CLOSE(80)
      ENDDO  
      
      
      
      ! Write the local to global node table
      DO pe = 1,nproc
        WRITE(dirname(3:lname),"(I4.4)") pe-1      
        OPEN(UNIT=81,FILE=dirname(1:lname)//'/'//'fort.81') 
        
        WRITE(81,*) nproc
        WRITE(81,*) nn
        WRITE(81,*) MAXVAL(nresnd)
        WRITE(81,*) nresnd(pe)
        WRITE(81,*) mndof
        WRITE(81,*) nout_sol
        
        DO el = 1,nresnd(pe)
          WRITE(81,*) nd_l2g(el,pe) 
        ENDDO
        CLOSE(81)
      ENDDO        
      
      
      CALL time_snaps(sta_opt,sta_snap,tf,dt,tskp_sta,nout_sta)      
      
      ! Write the local to global stations table
      DO pe = 1,nproc
        WRITE(dirname(3:lname),"(I4.4)") pe-1      
        OPEN(UNIT=82,FILE=dirname(1:lname)//'/'//'fort.82') 
        
        WRITE(82,*) nproc
        WRITE(82,*) nsta
        WRITE(82,*) MAXVAL(nlsta)
        WRITE(82,*) nlsta(pe)
        WRITE(82,*) nout_sta
        
        DO sta = 1,nlsta(pe)
          WRITE(82,*) sta_l2g(sta,pe) 
        ENDDO
        CLOSE(82)
      ENDDO             
      
      
     ! Write the high order bathymetry file
     IF (hb_file_exists) THEN
     DO pe = 1,nproc
       
       WRITE(dirname(3:lname),"(I4.4)") pe-1      
       OPEN(UNIT=14,FILE=dirname(1:lname)//'/'//'fort.hb')
       WRITE(14,"(A)") "decomposed from: " // bathy_file   
       WRITE(14,"(A)") "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"        
       WRITE(14,"(2(I7))") nresel(pe),hbp
       DO el = 1,nresel(pe)
         gel = el_l2g(el,pe)
         et = el_type(gel)
         
         nnd = nnds(order(et+4))
         
         WRITE(14,"(2(I7),1x,60(D24.17,1x))") el, nnd, (elhb(nd,gel), nd = 1,nnd)
       ENDDO
       
       CLOSE(14)
     ENDDO
     ENDIF
     
     
     ! Write the curved boundary edge file
     IF (cb_file_exists) THEN
     DO pe = 1,nproc
     
       WRITE(dirname(3:lname),"(I4.4)") pe-1      
       OPEN(UNIT=14,FILE=dirname(1:lname)//'/'//'fort.cb')     
        WRITE(14,"(A)") "decomposed from: " // curve_file
        WRITE(14,"(A)") "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"                 
        WRITE(14,"(I8,19x,A)") lnbou(pe), "! number of normal flow boundaries"
        WRITE(14,"(I8,I8,19x,A)") lnvel(pe), ctp, "! total number of normal flow nodes"
        
        i = 0
        DO bnd = 1,lnbou(pe)
          IF(lfbseg(1,bnd,pe) > 0) THEN
            i = i + 1
            segtype = lfbseg(2,bnd,pe) 
            n = lfbseg(1,bnd,pe)
            IF (segtype == 0 .OR. segtype == 10 .OR. segtype == 20 .OR. &
                segtype == 1 .OR. segtype == 11 .OR. segtype == 21) THEN
              WRITE(14,"(I8,1X,I8,10x,A,1x,I8)") n,segtype, "! number of nodes in normal flow boundary", i
              DO nd = 1,n-1
                lnd = lfbnds(nd,bnd,pe)
                gnd = nd_l2g(lnd,pe)
                WRITE(14,"(I8,1X,10(D24.17,1X))") lnd, xy(1,gnd), (lbndxy(1,k,nd,bnd,pe), k=1,ctp-1)
                WRITE(14,"(I8,1X,10(D24.17,1X))") lnd, xy(2,gnd), (lbndxy(2,k,nd,bnd,pe), k=1,ctp-1)
              ENDDO
              lnd = lfbnds(n,bnd,pe)
              gnd = nd_l2g(lnd,pe)              
              WRITE(14,"(I8,1X,D24.17,1X)") lnd, xy(1,gnd)
              WRITE(14,"(I8,1X,D24.17,1X)") lnd, xy(2,gnd)            
              
            ELSE
            
              WRITE(14,"(I8,1x,I8,10x,A,1x,I8)") 0,segtype, "! number of nodes in normal flow boundary", i
            
            ENDIF
          ENDIF
        ENDDO
        
        CLOSE(14)
       
     ENDDO
     ENDIF
     
     
     OPEN(UNIT=89, FILE='fort.89')
     DO el = 1,ne
       WRITE(89,*) el,el_g2l(1,el),el_g2l(2,el)
     ENDDO
     CLOSE(89)
            
      
      ! Write the local input file
      DO pe = 1,nproc
!         WRITE(dirname(3:lname),"(I4.4)") pe-1         
!         OPEN(UNIT=10,FILE=dirname(1:lname)//'/'//'dgswe.inp')
!         
!         WRITE(10,"(A)") dirname(1:lname)//'/'//"fort.14"
!         WRITE(10,"(A)") dirname(1:lname)//'/'//"fort.15"
!         WRITE(10,"(I5)") p
!         WRITE(10,"(I5)") ctp
!         WRITE(10,"(I5)") hbp
!         WRITE(10,"(E24.17)") dt
!         WRITE(10,"(E24.17)") tf
!         WRITE(10,"(E24.17)") dramp
!         WRITE(10,"(E24.17)") cf
!         WRITE(10,"(I7)") int(lines)
!         WRITE(10,"(A)") './' // dirname(1:lname) // '/'
!         WRITE(10,"(I5)") npart
!         
!         CLOSE(10)

        CALL write_local(pe)
      ENDDO

      
 180  FORMAT(8X,9I8)
 181  FORMAT(8X,2I8) 
 182  FORMAT(8X,I8)
 183  FORMAT(8X,4(D24.17,1X))
 
      RETURN
      END SUBROUTINE write_files
