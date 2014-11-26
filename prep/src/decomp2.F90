      SUBROUTINE decomp2()
      
      USE globals, ONLY: nn,ne,ned,part,nproc,ect,nelnds,mnelnds,xy,lect,lnelnds, &
                         ged2el,ged2led,&
                         nsred,sredn, &
                         nresel,el_g2l,el_l2g, &
                         mned_sr,ned_sr,pe_sr,el_sr,led_sr, &
                         nresnd,nd_g2l,nd_l2g, &
                         nope,neta,obseg,obnds, &
                         nbou,nvel,fbseg,fbnds, &
                         lnope,lneta,lobseg,lobnds, &
                         lnbou,lnvel,lfbseg,lfbnds, &
                         nobfr,obamp,obph,lobamp,lobph, &
                         nfbfr,fbamp,fbph,lfbamp,lfbph,lnbouf

      IMPLICIT NONE

      INTEGER :: el,pe,nd,ed,j,bnd,eln
      INTEGER :: ged,lnd
      INTEGER :: el1,el2
      INTEGER :: pe1,pe2
      INTEGER :: led1,led2
      INTEGER :: ndcnt,bfnd
      INTEGER :: mnepe
      INTEGER :: mnobnds,mnfbnds
      INTEGER :: segtype
      INTEGER :: bnd_flag
      INTEGER :: lbnd
      
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ndflag 
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: bou_l2g
      
      ALLOCATE(nresel(nproc))
      ALLOCATE(el_g2l(2,ne))
      ALLOCATE(el_l2g(ne,nproc))
      
      nresel = 0
      DO el = 1,ne
        pe = part(el) 
        nresel(pe) = nresel(pe) + 1 ! count resident elements in subdomain pe
        
        el_g2l(1,el) = pe ! global element el is in subdomain pe
        el_g2l(2,el) = nresel(pe) ! global element el has local element number nresel(pe) in subdomain pe
        el_l2g(nresel(pe),pe) = el ! local element nresel(pe) on subdomain pe has global element number el
      ENDDO
      
      mnepe = MAXVAL(nresel)     
      
      
      ! Find the commincation edges and 
      ! maximunm number of send/recv edges for allocation
      ALLOCATE(ned_sr(nproc))
      ALLOCATE(sredn(ned))
      
      ned_sr = 0
      nsred = 0
      DO ed = 1,ned 
        el1 = ged2el(1,ed) ! find elements on this edge
        el2 = ged2el(2,ed)
        
        IF(el2 /= 0) THEN ! don't want edges on the boundary of the global domain
          pe1 = el_g2l(1,el1) ! find subdomain numbers for each element
          pe2 = el_g2l(1,el2)
    
          IF(pe1 == pe2) THEN
            ! resident element, do nothing
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

      ENDDO
      
      
      
      ! Determine local element connectivity table and
      ! find the local node numbers
      ! find local global domain boundary nodes
      ALLOCATE(nresnd(nproc))
      ALLOCATE(ndflag(nn))
      ALLOCATE(lect(mnelnds,mnepe,nproc))
      ALLOCATE(lnelnds(mnepe,nproc))
      ALLOCATE(nd_g2l(nn))
      ALLOCATE(nd_l2g(mnelnds*mnepe,nproc))
      
      mnobnds = MAXVAL(obseg)
      ALLOCATE(lobseg(nope,nproc))
      ALLOCATE(lobnds(mnobnds,nope,nproc))
      ALLOCATE(lneta(nproc))
      ALLOCATE(lnope(nproc))
      ALLOCATE(lobamp(neta,nobfr,nproc))
      ALLOCATE(lobph(neta,nobfr,nproc))
      
      mnfbnds = MAXVAL(fbseg(1,:))
      ALLOCATE(lfbseg(2,nbou,nproc))
      ALLOCATE(lfbnds(mnfbnds,nbou,nproc))
      ALLOCATE(lnvel(nproc))
      ALLOCATE(lnbou(nproc))
      ALLOCATE(lfbamp(mnfbnds,nbou,nfbfr,nproc))
      ALLOCATE(lfbph(mnfbnds,nbou,nfbfr,nproc))
      ALLOCATE(lnbouf(nproc))
      ALLOCATE(bou_l2g(nbou,nproc))
      
      lobseg = 0
      lfbseg = 0
      DO pe = 1,nproc
        ndflag = 0
        nresnd(pe) = 0
        nd_g2l = 0 ! some nodes will be included in 2 subdomains so need to start a new global to local table each time
        DO el = 1,nresel(pe) 
          eln = el_l2g(el,pe)
          lnelnds(el,pe) = nelnds(eln)
          DO j = 1,nelnds(eln) ! loop through each node of all elements on subdomain pe
            nd = ect(j,eln) ! find global node number
            IF(ndflag(nd) == 0) THEN  ! decide if it's been counted already
              nresnd(pe) = nresnd(pe) + 1 ! count as a resident node
              nd_l2g(nresnd(pe),pe) = nd  ! local node nresnd(pe) on subdomain pe is global node nd
              nd_g2l(nd) = nresnd(pe) ! global node nd is local node nresnd(pe)
              ndflag(nd) = 1 ! flag the node so it's not counted again
              
              lect(j,el,pe) = nresnd(pe) ! fill in the local element connectivity table
            ELSE
              lect(j,el,pe) = nd_g2l(nd) ! if the node has already been assigned a local number, fill in the local connectivity table
            ENDIF
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
              lobnds(lobseg(bnd,pe),bnd,pe) = lnd 
              
              lobamp(lneta(pe),:,pe) = obamp(ndcnt,:)
              lobph(lneta(pe),:,pe) = obph(ndcnt,:)
            ENDIF
          ENDDO
          IF(lobseg(bnd,pe) > 0) THEN
            lnope(pe) = lnope(pe) + 1
          ENDIF
        ENDDO
        
        lnvel(pe) = 0 ! look for flow boundary nodes in this subdomain
        lnbou(pe) = 0
        lnbouf(pe) = 0
        
        DO bnd = 1,nbou
          segtype = fbseg(2,bnd)          
          bfnd = 0
          bnd_flag = 0 
          DO j = 1,fbseg(1,bnd)
            nd = fbnds(j,bnd)
            lnd = nd_g2l(nd)
            IF(lnd == 0) THEN
            
              ! node is not in this subdomain
              IF (bnd_flag == 1) THEN
                bnd_flag = 0 ! local boundary is broken
              ENDIF
              
            ELSE
            
              IF (bnd_flag == 0) THEN 
                ! start new local boundary
                lnbou(pe) = lnbou(pe) + 1 
                lbnd = lnbou(pe)
                bou_l2g(lbnd,pe) = bnd
                
                lfbseg(2,lbnd,pe) = segtype
                bnd_flag = 1
              ENDIF
             
              lnvel(pe) = lnvel(pe) + 1
              lfbseg(1,lbnd,pe) = lfbseg(1,lbnd,pe) + 1
              lfbnds(lfbseg(1,lbnd,pe),lbnd,pe) = lnd              
              
            ENDIF
          ENDDO                    

        ENDDO
        
        DO lbnd = 1,lnbou(pe)
          bnd = bou_l2g(lbnd,pe)
          IF(lfbseg(1,lbnd,pe) > 1) THEN ! count local flow boundaries
            segtype = lfbseg(2,lbnd,pe)
            
            IF(segtype == 2 .OR. segtype == 12 .OR. segtype == 22) THEN
              lnbouf(pe) = lnbouf(pe) + 1 ! count local forced flow boundaries
              DO bfnd = 1,lfbseg(1,lbnd,pe)
                j = nd_l2g(bfnd,pe)
                lfbamp(bfnd,lbnd,:,pe) = fbamp(j,bnd,:)
                lfbph(bfnd,lbnd,:,pe) = fbph(j,bnd,:) 
              ENDDO
            ENDIF

            IF(segtype == 1 .OR. segtype == 11 .OR. segtype == 21) THEN 
              IF(lfbseg(1,lbnd,pe) /= fbseg(1,bnd)) THEN  ! if the entire island boundary is not contained in the subdomain
                lfbseg(2,lbnd,pe) = 10                    ! change it to a land boundary              
                PRINT("(A,I7,A,I7,A)"), "Island boundary: ", lnbou(pe), " on PE: ", pe-1, " changed to land"
              ENDIF         
            ENDIF
            
          ELSE IF(lfbseg(1,lbnd,pe) == 1) THEN
            PRINT*, "Error: boundary only has one node"           
          ENDIF     
        ENDDO
        
      ENDDO
      
      
      RETURN 
      END SUBROUTINE decomp2
