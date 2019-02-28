      SUBROUTINE decomp2()
      
      USE globals, ONLY: rp,nn,ne,ned,part,ect,nverts,lect,lnelnds,mnnds, &
                         ged2el,ged2led,ged2nn,el_type, &
                         nsred,sredn, &
                         nbed,bedn, &
                         nresel,el_g2l,el_l2g, &
                         mned_sr, &
                         nresnd,nd_g2l,nd_l2g, &
                         nope,neta,obseg,obnds, &
                         nbou,nvel,fbseg,fbnds, &
                         lnope,lneta,lobseg,lobnds, &
                         lnbou,lnvel,lfbseg,lfbnds, &
                         nobfr,obamp,obph, &
                         lobamp,lobph, &
                         nfbfr,nfbsfr,fbamp,fbph,fbsamp, &
                         lfbamp,lfbph,lnbouf,lfbsamp, &
                         nx_pt,ny_pt,detJe, &
                         hbqpted,nqpte,mnqpte,elhb, &
                         lbndxy,bndxy, &
                         nsta,elsta,xysta, &
                         nlsta,sta_l2g, &
                         elxy                         
 
      USE read_dginp, ONLY: hbp,ctp                           
                         
      USE messenger2, ONLY: nproc,nqpte_sr,hb_sr,nx_sr,ny_sr,detJe_sr, &
                            ned_sr,pe_sr,el_sr,led_sr,nied_pe                       

      IMPLICIT NONE

      INTEGER :: el,pe,nd,ed,i,j,k,bnd,eln,bfr,pt,sta,bou
      INTEGER :: ged,lnd,nv,n
      INTEGER :: lled,gled,ln1,ln2,gn1,gn2,n1,n2,nd1,nd2
      INTEGER :: el1,el2
      INTEGER :: pe1,pe2
      INTEGER :: led1,led2
      INTEGER :: et1,et2,et
      INTEGER :: gp_in,gp_ex
      INTEGER :: ndcnt,bfnd,nlbnds
      INTEGER :: mnepe
      INTEGER :: mnobnds,mnfbnds
      INTEGER :: segtype
      INTEGER :: bnd_flag,stop_flag,el_flag
      INTEGER :: lbou,mlbou
      
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ndflag 
      INTEGER, ALLOCATABLE, DIMENSION(:) :: sta_flag
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: bou_l2g
      
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: elhb_old
      
       
      ALLOCATE(elhb_old(mnnds,ne))
      elhb_old = elhb
      
      ALLOCATE(nresel(nproc))
      ALLOCATE(el_g2l(2,ne))
      
      nresel = 0
      DO el = 1,ne
        pe = part(el) 
        nresel(pe) = nresel(pe) + 1 ! count resident elements in subdomain pe
      ENDDO
      
      mnepe = MAXVAL(nresel)          
      ALLOCATE(el_l2g(mnepe,nproc))

      
      nresel = 0
      DO el = 1,ne
        pe = part(el) 
        nresel(pe) = nresel(pe) + 1 ! count resident elements in subdomain pe
        
        el_g2l(1,el) = pe ! global element el is in subdomain pe
        el_g2l(2,el) = nresel(pe) ! global element el has local element number nresel(pe) in subdomain pe
        el_l2g(nresel(pe),pe) = el ! local element nresel(pe) on subdomain pe has global element number el
      ENDDO
      
         
      
      ! Find the commincation edges and 
      ! maximunm number of send/recv edges for allocation
      ALLOCATE(ned_sr(nproc))
      ALLOCATE(sredn(ned))
      ALLOCATE(nied_pe(nproc))

      ned_sr = 0
      nied_pe = 0
      nsred = 0
      DO ed = 1,ned 
        el1 = ged2el(1,ed) ! find elements on this edge
        el2 = ged2el(2,ed)
        
        IF(el2 /= 0) THEN ! don't want edges on the boundary of the global domain
          pe1 = el_g2l(1,el1) ! find subdomain numbers for each element
          pe2 = el_g2l(1,el2)
    
          IF(pe1 == pe2) THEN
            ! resident element, do nothing
            nied_pe(pe1) = nied_pe(pe1) + 1 ! keep track of interior edges
            nied_pe(pe2) = nied_pe(pe2) + 1
          ELSE 
            ned_sr(pe1) = ned_sr(pe1) + 1 ! count send/recieve edges for each subdomain
            ned_sr(pe2) = ned_sr(pe2) + 1
            
            nsred = nsred + 1 ! count communication edges
            sredn(nsred) = ed ! keep track of communication edges
          ENDIF
        ENDIF
        
      ENDDO
      
      mned_sr = MAXVAL(ned_sr) ! find the maximum number of communication edges
      
      ALLOCATE(pe_sr(mned_sr,nproc),el_sr(mned_sr,nproc),led_sr(mned_sr,nproc))
      ALLOCATE(nx_sr(mnqpte,mned_sr,nproc),ny_sr(mnqpte,mned_sr,nproc))
      ALLOCATE(hb_sr(mnqpte,mned_sr,nproc))
      ALLOCATE(detJe_sr(mnqpte,mned_sr,nproc))
      ALLOCATE(nqpte_sr(mned_sr,nproc))

      
      

      
      ! Loop through communication edges and keep track of what local edge number from what element gets 
      ! sent/recieved to/from which subdomains
      ned_sr = 0
      DO ed = 1,nsred
      
        ged = sredn(ed)
        
        el1 = ged2el(1,ged) ! find elements on this edge
        el2 = ged2el(2,ged)
        
        led1 = ged2led(1,ged) ! find local edge numbers on this edge
        led2 = ged2led(2,ged) 

        pe1 = el_g2l(1,el1) ! find subdomain number for each element
        pe2 = el_g2l(1,el2)
          
        ned_sr(pe1) = ned_sr(pe1) + 1 ! count send/recieve edges for each subdomain
        ned_sr(pe2) = ned_sr(pe2) + 1
          
        pe_sr(ned_sr(pe1),pe1) = pe2 ! pe1 sends/recieves this edge to/from pe2
        pe_sr(ned_sr(pe2),pe2) = pe1 ! pe2 sends/recieves this edge to/from pe1
        
        el_sr(ned_sr(pe1),pe1) = el_g2l(2,el1) ! pe1 sends/recieves edge from el1 to/from pe2
        el_sr(ned_sr(pe2),pe2) = el_g2l(2,el2) ! pe2 sends/recieves edge from el2 to/from pe1
      
        led_sr(ned_sr(pe1),pe1) = led1 ! pe1 sends local edge led1 from el1 to pe2
        led_sr(ned_sr(pe2),pe2) = led2 ! pe2 sends local edge led2 from el2 to pe1
        
!         et1 = el_type(el1)
!         et2 = el_type(el2)
!         
!         et = max(et1,et2)
        
!         nqpte_sr(ned_sr(pe1),pe1) = nqpte(et)
!         nqpte_sr(ned_sr(pe2),pe2) = nqpte(et)

        nqpte_sr(ned_sr(pe1),pe1) = nqpte(1)  ! assume all internal edges are straight
        nqpte_sr(ned_sr(pe2),pe2) = nqpte(1)
        
        DO pt = 1,nqpte(1) 
        
          gp_in = pt
          gp_ex = nqpte(1) + 1 - pt
          
          nx_sr(pt,ned_sr(pe1),pe1) = nx_pt(ged,gp_in)
          nx_sr(pt,ned_sr(pe2),pe2) = -nx_pt(ged,gp_ex)  
        
          ny_sr(pt,ned_sr(pe1),pe1) = ny_pt(ged,gp_in)
          ny_sr(pt,ned_sr(pe2),pe2) = -ny_pt(ged,gp_ex)   
          
          hb_sr(pt,ned_sr(pe1),pe1) = hbqpted(ged,gp_in)
          hb_sr(pt,ned_sr(pe2),pe2) = hbqpted(ged,gp_ex)
          
          detJe_sr(pt,ned_sr(pe1),pe1) = detJe(ged,gp_in)
          detJe_sr(pt,ned_sr(pe2),pe2) = detJe(ged,gp_ex)
        ENDDO
        
      ENDDO
      
      DEALLOCATE(nx_pt,ny_pt)
      DEALLOCATE(hbqpted)
      DEALLOCATE(detJe)
      DEALLOCATE(elxy)

      
      
      
      ! Determine local element connectivity table and
      ! find the local node numbers
      ! find local global domain boundary nodes
      ALLOCATE(nresnd(nproc))
      ALLOCATE(ndflag(nn))
      ALLOCATE(lect(4,mnepe,nproc))
      ALLOCATE(lnelnds(mnepe,nproc))
      ALLOCATE(nd_g2l(nn))
      ALLOCATE(nd_l2g(4*mnepe,nproc))
      
      mnobnds = MAXVAL(obseg)
      ALLOCATE(lobseg(nope,nproc))
      ALLOCATE(lobnds(mnobnds,nope,nproc))
      ALLOCATE(lneta(nproc))
      ALLOCATE(lnope(nproc))
      ALLOCATE(lobamp(mnobnds,nope,nobfr,nproc))
      ALLOCATE(lobph(mnobnds,nope,nobfr,nproc))
      
     
      mnfbnds = 0                                  ! Find the max number of flow boundary nodes
      DO pe = 1,nproc                              ! in any subdomain.
        n = 0
        ndflag = 0        
        DO el = 1,nresel(pe)
          eln = el_l2g(el,pe)
          et = el_type(eln)
          DO j = 1,nverts(et)
            nd = ect(j,eln)           
              
            IF (ndflag(nd) == 0) THEN
              ndflag(nd) = 1
              
   nd_search: DO bou = 1,nbou
                DO i = 1,fbseg(1,bou)
                  IF (nd == fbnds(i,bou)) THEN
                    n = n+1
                    EXIT nd_search
                  ENDIF
                ENDDO
              ENDDO nd_search
              
            ENDIF
              
          ENDDO
        ENDDO
        
        IF (n > mnfbnds) THEN
          mnfbnds = n
        ENDIF
      ENDDO
      
!       mnfbnds = MAXVAL(fbseg(1,:))
      mlbou = nbou+2                               ! this is a guess for the upper bound, nbou is not enough because boundaries can get split up.      
      
      ALLOCATE(lfbseg(2,mlbou,nproc))
      ALLOCATE(lfbnds(mnfbnds,mlbou,nproc))
      ALLOCATE(lbndxy(2,ctp+1,mnfbnds,mlbou,nproc))      
      ALLOCATE(lnvel(nproc))
      ALLOCATE(lnbou(nproc))
      ALLOCATE(lfbamp(mnfbnds,mlbou,nfbfr,nproc))
      ALLOCATE(lfbph(mnfbnds,mlbou,nfbfr,nproc))
      ALLOCATE(lfbsamp(mnfbnds,mlbou,nfbfr,nproc))      
      ALLOCATE(lnbouf(nproc))
      ALLOCATE(bou_l2g(mlbou,nproc))
      
      lobseg = 0
      lfbseg = 0
  pes:DO pe = 1,nproc

        ndflag = 0
        nresnd(pe) = 0
        nd_g2l = 0 ! some nodes will be included in 2 subdomains so need to start a new global to local table each time
!         DO el = 1,nresel(pe) 
!           eln = el_l2g(el,pe)
!           lnelnds(el,pe) = nelnds(eln)
!           DO j = 1,nelnds(eln) ! loop through each node of all elements on subdomain pe
!             nd = ect(j,eln) ! find global node number
!             IF(ndflag(nd) == 0) THEN  ! decide if it's been counted already
!               nresnd(pe) = nresnd(pe) + 1 ! count as a resident node
!               nd_l2g(nresnd(pe),pe) = nd  ! local node nresnd(pe) on subdomain pe is global node nd
!               nd_g2l(nd) = nresnd(pe) ! global node nd is local node nresnd(pe)
!               ndflag(nd) = 1 ! flag the node so it's not counted again
!               
!               lect(j,el,pe) = nresnd(pe) ! fill in the local element connectivity table
!             ELSE
!               lect(j,el,pe) = nd_g2l(nd) ! if the node has already been assigned a local number, fill in the local connectivity table
!             ENDIF
!           ENDDO
!         ENDDO
        
        
        DO el = 1,nresel(pe) 
          eln = el_l2g(el,pe)
          et = el_type(eln)
          lnelnds(el,pe) = nverts(et)
          DO j = 1,nverts(et) ! loop through each node of all elements on subdomain pe
            nd = ect(j,eln)   ! find global node number
            IF(ndflag(nd) == 0) THEN  ! decide if it's been counted already
              nresnd(pe) = nresnd(pe) + 1 ! count as a resident node
              nd_l2g(nresnd(pe),pe) = nd  ! local node nresnd(pe) on subdomain pe is global node nd
              nd_g2l(nd) = nresnd(pe)     ! global node nd is local node nresnd(pe)
              ndflag(nd) = 1              ! flag the node so it's not counted again
              
              lect(j,el,pe) = nresnd(pe) ! fill in the local element connectivity table
            ELSE
              lect(j,el,pe) = nd_g2l(nd) ! if the node has already been assigned a local number, fill in the local connectivity table
            ENDIF
          ENDDO
        ENDDO        
        
        DO el = 1,nresel(pe)
          eln = el_l2g(el,pe)
          et = el_type(eln)
          nv = nverts(et)
          DO lled = 1,nv
            ln1 = lect(mod(lled+0,nv)+1,el,pe)
            ln2 = lect(mod(lled+1,nv)+1,el,pe)
            
            gn1 = nd_l2g(ln1,pe)
            gn2 = nd_l2g(ln2,pe)
            
            DO gled = 1,nv
              n1 = ect(mod(gled+0,nv)+1,eln)
              n2 = ect(mod(gled+1,nv)+1,eln)              
              
              IF (((n1 == gn1) .AND. (n2 == gn2)) .OR. &
                  ((n1 == gn2) .AND. (n2 == gn1))) THEN
                IF (gled /= lled) THEN
                  PRINT*, "LOCAL EDGE CHANGED", eln,gled,lled
                  DO pt = 1,hbp
                    i = mod(lled,nv)*hbp + pt 
                    j = mod(gled,nv)*hbp + pt 
                    elhb(i,eln) = elhb_old(j,eln)
                  ENDDO
                  
                ELSE  
!                   PRINT*, "LOCAL EDGE SAME"
                ENDIF
                
              ENDIF
            ENDDO
            
          ENDDO          
        ENDDO
        
        
        
        
        lneta(pe) = 0 ! look for open boundary nodes in this subdomain
        lnope(pe) = 0
        ndcnt = 0
        DO bnd = 1,nope
          DO j = 1,obseg(bnd)
            ndcnt = ndcnt + 1
            nd = obnds(j,bnd)
            lnd = nd_g2l(nd)
            IF(lnd == 0) THEN
              ! skip, node is not in this subdomain
            ELSE
              lneta(pe) = lneta(pe) + 1
              lobseg(bnd,pe) = lobseg(bnd,pe) + 1
              nlbnds = lobseg(bnd,pe)
              lobnds(nlbnds,bnd,pe) = lnd 
              
              IF(lobseg(bnd,pe) == 1) THEN
                lnope(pe) = lnope(pe) + 1
                lbou = lnope(pe)
              ENDIF              
              
              DO bfr = 1,nobfr             
                lobamp(nlbnds,lbou,bfr,pe) = obamp(j,bnd,bfr)
                lobph(nlbnds,lbou,bfr,pe) = obph(j,bnd,bfr)              
              ENDDO
              
            ENDIF
          ENDDO
        ENDDO
        
        
        
        
        lnvel(pe) = 0     
        lnbou(pe) = 0
        lnbouf(pe) = 0
        lbou = 1
  bnds: DO bnd = 1,nbou            ! look for flow boundary nodes in this subdomain            
          segtype = fbseg(2,bnd)   !   - this has to take into account that global island boundaries my be broken up into several local   
          bfnd = 0                 !     boundary segments, meaning there may be more than one local boundary per global boundary          
          bnd_flag = 0             !   - when this happens the island boundaries need to be changed to land boundaries because they are no
     nds: DO j = 1,fbseg(1,bnd)    !     longer closed.
            nd1 = fbnds(j,bnd)
            lnd = nd_g2l(nd1)
            
           
            el_flag = 1                                ! Find if element edge is in this subdomain, 
            IF (j <= fbseg(1,bnd)-1) THEN              ! to determine if this breaks up a boundary. 
              nd2 = fbnds(j+1,bnd)                     ! (Assume el_flag=1 to start to account for when j==fbseg(1,bnd))
                                                       ! This is an issue when both nodes belong to this subdomain, but the              
       edges: DO ed = 1,nbed                           ! edge between them belongs to an element from another.
                ged = bedn(ed)
                n1 = ged2nn(1,ged)
                n2 = ged2nn(2,ged)
                
                IF(((nd1 == n1).AND.(nd2 == n2)).OR. &
                   ((nd1 == n2).AND.(nd2 == n1))) THEN
                   
                  el = ged2el(1,ged)                  
                  EXIT edges
                ENDIF
              ENDDO edges
              
              IF (el_g2l(1,el) /= pe) THEN
                el_flag = 0
              ENDIF
            ENDIF
                 
            
            IF(lnd == 0) THEN       
                        
              ! node is not in this subdomain
              
              IF (bnd_flag == 1) THEN
                bnd_flag = 0 ! local boundary is broken

                IF (lfbseg(1,lbou,pe) == 1) THEN ! boundary might break after one node
                  lnbou(pe) = lnbou(pe) - 1      ! if so, don't count it
                  lfbseg(1,lbou,pe) = 0          ! (that node will be included at the end of the global boundary)
                ENDIF
              ENDIF
              
            ELSE
            
              IF (bnd_flag == 0) THEN ! if the boundary was broken, start a new local boundary
              
                IF (j == fbseg(1,bnd)) THEN ! unless it's the "last" node in the global boundary
                  EXIT nds                  ! (island boundaries are closed, i.e. they begin and end with same node,
                ENDIF                       !  so this is really the first node has already been included earlier)
                                        
                lnbou(pe) = lnbou(pe) + 1   ! start a new local boundary
                lbou = lnbou(pe)
                bou_l2g(lbou,pe) = bnd      ! keep track of what global boundary each local segement is from
                
                lfbseg(2,lbou,pe) = segtype
                bnd_flag = 1                ! boundary is no longer broken               
              ENDIF
             
              lnvel(pe) = lnvel(pe) + 1                   ! increment total # of boundary nodes
              lfbseg(1,lbou,pe) = lfbseg(1,lbou,pe) + 1   ! inrement # of segment boundary nodes
              nlbnds = lfbseg(1,lbou,pe) 
              lfbnds(nlbnds,lbou,pe) = lnd                ! keep track of local boundary node numbers
              
              IF (el_flag == 0) THEN
                bnd_flag = 0                     ! break boundary if element edge is not in domain
                IF (lfbseg(1,lbou,pe) == 1) THEN ! boundary might break after one node
                  lnbou(pe) = lnbou(pe) - 1      ! if so, don't count it
                  lfbseg(1,lbou,pe) = 0          ! (that node will be included at the end of the global boundary)
                ENDIF                
              ENDIF
              
              DO k = 1,ctp-1
                lbndxy(1,k,nlbnds,lbou,pe) = bndxy(1,k,j,bnd)
                lbndxy(2,k,nlbnds,lbou,pe) = bndxy(2,k,j,bnd)
              ENDDO
              
              IF(segtype == 2 .OR. segtype == 12 .OR. segtype == 22) THEN
                lnbouf(pe) = lnbouf(pe) + 1   ! count local forced flow boundaries
                DO bfr = 1,nfbfr              ! keep track of forcings
                  lfbamp(nlbnds,lbou,bfr,pe) = fbamp(j,bnd,bfr)
                  lfbph(nlbnds,lbou,bfr,pe) = fbph(j,bnd,bfr) 
                ENDDO
                DO bfr = 1,nfbsfr
                  lfbsamp(nlbnds,lbou,bfr,pe) = fbsamp(j,bnd,bfr)                  
                ENDDO
              ENDIF              
              
            ENDIF
          ENDDO nds                
          
          
        ENDDO bnds
        
        
        
        
        
        
        
        
        
        
        
        DO lbou = 1,lnbou(pe) ! check local flow boundaries      
          bnd = bou_l2g(lbou,pe) ! global boundary that corresponds to local boundary
          
          IF(lfbseg(1,lbou,pe) > 1) THEN
          
            segtype = lfbseg(2,lbou,pe)             
            IF(segtype == 1 .OR. segtype == 11 .OR. segtype == 21) THEN 
              IF(lfbseg(1,lbou,pe) /= fbseg(1,bnd)) THEN  ! if the entire island boundary is not contained in the subdomain
                lfbseg(2,lbou,pe) = 10                    ! change it to a land boundary              
                PRINT("(A,I7,A,I7,A,I7,A,I7,A)"), "Island boundary: ", lbou," / ",lnbou(pe), " on PE: ", pe-1, &
                                                  " changed to land (Global boundary ", bnd, ")"
              ENDIF         
            ENDIF
            
          ELSE IF(lfbseg(1,lbou,pe) == 1) THEN
            PRINT("(A,I7,A)"), "Error: boundary ",lbou," only has one node"          
          ENDIF  
        ENDDO       
        
        
      ENDDO pes
      
      
      
      
      
      
      ALLOCATE(nlsta(nproc))
      ALLOCATE(sta_l2g(nsta,nproc))
      ALLOCATE(sta_flag(nsta))
      
      nlsta = 0
      sta_flag = 0
      DO pe = 1,nproc      
        DO el = 1,nresel(pe) 
          eln = el_l2g(el,pe)
          
          DO sta = 1,nsta   
            IF (elsta(sta) == eln) THEN
            
              nlsta(pe) = nlsta(pe) + 1
              sta_l2g(nlsta(pe),pe) = sta
              sta_flag(sta) = sta_flag(sta) + 1
              
            ENDIF
          ENDDO 
          
        ENDDO      
      ENDDO
      
      stop_flag = 0
      DO sta = 1,nsta
        IF (sta_flag(sta) /= 1) THEN
          PRINT*, "station: ",sta," not found placed in pe",sta_flag(sta), elsta(sta)
          stop_flag = 1
        ENDIF
      ENDDO
      
!       IF (stop_flag) STOP
      
      
      RETURN 
      END SUBROUTINE decomp2
