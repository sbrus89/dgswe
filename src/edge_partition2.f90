      SUBROUTINE metis2()

      USE globals, ONLY: nn,ne,ned,ged2el,npart,part

      IMPLICIT NONE
      
      INTEGER :: ed,nd,j,el
      INTEGER :: el1,el2,nd2
      INTEGER :: tot
      INTEGER :: tmp
      
      INTEGER, ALLOCATABLE, DIMENSION(:) :: nadjels
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: adjels
      
      INTEGER :: numflag,wgtflag
      INTEGER :: edgecut
      INTEGER :: nparts
      INTEGER, ALLOCATABLE, DIMENSION(:) :: xadj
      INTEGER, ALLOCATABLE, DIMENSION(:) :: adjncy
      INTEGER, ALLOCATABLE, DIMENSION(:) :: vwgt
      INTEGER, ALLOCATABLE, DIMENSION(:) :: adjwgt
      INTEGER :: options(5)
      
      edgecut = 0
      
      ! Find the nodes adjacent each node
      ALLOCATE(nadjels(ne))
      ALLOCATE(adjels(4,ne))
      
      nadjels = 0
      DO ed = 1,ned
        el1 = ged2el(1,ed) ! find the node numbers on each edge
        el2 = ged2el(2,ed)
        
        IF(el2 /= 0) THEN
          nadjels(el1) = nadjels(el1) + 1 ! count the nodes adjacent to node n1
          nadjels(el2) = nadjels(el2) + 1 ! count the nodes adjacent to node n2
        
          adjels(nadjels(el1),el1) = el2 ! node n2 is adjacent to node n1
          adjels(nadjels(el2),el2) = el1 ! node n1 is adjacent to node n2       
        ENDIF
      ENDDO
      
      ! sort the node numbers adjacent to each node 
      ! (I'm not sure why this matters, but it seems to reduce the edge cut count)
      DO el = 1,ne
        DO ed = 1,nadjels(el)
          DO j = ed,nadjels(el)
            IF(adjels(j,el) < adjels(ed,el)) THEN
              tmp = adjels(j,el)
              adjels(j,el) = adjels(ed,el)
              adjels(ed,el) = tmp
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      
!       DO el = 1,ne
!         PRINT("(I7,5x,I7,5x,10(I5))"), el, nadjels(el), (adjels(ed,el), ed = 1,nadjels(el))
!       ENDDO      
      
      ! Compute the node weights
      ! (ADCPREP uses the number of adjacent nodes)
      ALLOCATE(vwgt(ne))
      
      DO el = 1,ne
        vwgt(el) = nadjels(el)
!         PRINT*, nd,vwgt(nd)
      ENDDO
      
      ! Create the adjacency arrays used by metis
      ! Also compute the edge weight 
      ! (ADCPREP uses the sum of the adjacent nodes for each node on the edge)
      ALLOCATE(xadj(ne+1))
      ALLOCATE(adjncy(2*ned),adjwgt(2*ned))
      
      xadj(1) = 1
      tot = 0
      DO el = 1,ne
        DO ed = 1,nadjels(el)
          tot = tot + 1
          el2 = adjels(ed,el)
          adjncy(tot) = el2 
          adjwgt(tot) = nadjels(el) + nadjels(el2)
!           PRINT*, tot, adjncy(tot), adjwgt(tot)
        ENDDO
        xadj(el+1) = tot+1
      ENDDO      
      
      numflag = 1 ! use fortran numbering
      wgtflag = 3 ! use edge and node weights
      
      options(1) = 1 ! use non-default options
      options(2) = 3 ! use sorted heavy edge matching type
      options(3) = 1 ! use multilevel recursive bisection algorithm during initial partitioning 
      options(4) = 3 ! minimize connectivity among the subdomains 
      options(5) = 0 ! always set to zero
      
      nparts = npart ! number of partitions
      
      ALLOCATE(part(ne))
      
      IF (nparts > 1) THEN
      
        CALL METIS_PartGraphKway(ne,xadj,adjncy,vwgt,adjwgt,wgtflag, &
                                 numflag,nparts,options,edgecut,part)
                               
      ELSE IF (nparts == 1) THEN
      
        DO el = 1,ne
          part(el) = 1
        ENDDO
      
      ENDIF
      
      PRINT*,"edgecut = ", edgecut
      
      RETURN 
      END SUBROUTINE metis2
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      SUBROUTINE decomp2()
      
      USE globals, ONLY: nn,ne,ndof,nqpta,part,npart,nel_type,el_type,nelnds, &
                         ged2el,&
                         npartel,nparted,npartet, &
                         gel2part,gel2lel,lel2gel, &
                         ael2gel,gel2ael, &
                         nied,iedn, &
                         H,Qx,Qy,Hinit,Qxinit,Qyinit, &
                         dpdx,dpdy,dpdx_init,dpdy_init, &
                         dhbdx,dhbdy,dhbdx_init,dhbdy_init, &
                         phia_int,phia_int_init, &
                         mmi,mmi_init, &
                         Hwrite,Qxwrite,Qywrite, &
                         nblk,elblk,edblk,nfblk,nrblk,rnfblk, &
                         mnpartel,mnparted   


      IMPLICIT NONE

      INTEGER :: el,pe,ed,dof,blk,et
      INTEGER :: ged
      INTEGER :: el_in,el_ex
      INTEGER :: pe_in,pe_ex
      INTEGER :: elcnt,edcnt
      INTEGER :: ted,tel
      INTEGER :: mndof,mnqpta
      
      INTEGER, ALLOCATABLE, DIMENSION(:) :: edflag      
      
      ALLOCATE(npartel(npart))
      ALLOCATE(gel2part(ne))
      ALLOCATE(gel2lel(ne))
      ALLOCATE(lel2gel(ne,npart))
      ALLOCATE(npartet(nel_type,npart))
      
      npartel = 0
      npartet = 0
      DO et = 1,nel_type
        DO el = 1,ne
          IF (el_type(el) == et) THEN
            pe = part(el) 
            npartel(pe) = npartel(pe) + 1 ! count elements in partition pe
            npartet(et,pe) = npartet(et,pe) + 1  
            
            gel2part(el) = pe ! global element el is in parition pe
            gel2lel(el) = npartel(pe) ! global element el has local element number npartel(pe) in partition pe
            lel2gel(npartel(pe),pe) = el ! local element npartel(pe) on partition pe has global element number el
          ENDIF
        ENDDO
      ENDDO
        
        
      ALLOCATE(ael2gel(ne),gel2ael(ne))
      
      mndof = maxval(ndof)
      mnqpta = maxval(nqpta)
      
      elcnt = 0 
      DO pe = 1,npart
!       PRINT*, " "
        DO el = 1,npartel(pe)
          elcnt = elcnt + 1
        
          H(elcnt,1:mndof) = Hinit(lel2gel(el,pe),1:mndof)
          Qx(elcnt,1:mndof) = Qxinit(lel2gel(el,pe),1:mndof)
          Qy(elcnt,1:mndof) = Qyinit(lel2gel(el,pe),1:mndof)  
          
          dpdx(elcnt,1:mndof*mnqpta) = dpdx_init(lel2gel(el,pe),1:mndof*mnqpta)
          dpdy(elcnt,1:mndof*mnqpta) = dpdy_init(lel2gel(el,pe),1:mndof*mnqpta)
          
          phia_int(elcnt,1:mndof*mnqpta) = phia_int_init(lel2gel(el,pe),1:mndof*mnqpta)
          
          dhbdx(elcnt,1:mnqpta) = dhbdx_init(lel2gel(el,pe),1:mnqpta)
          dhbdy(elcnt,1:mnqpta) = dhbdy_init(lel2gel(el,pe),1:mnqpta)
          
          mmi(elcnt,1:mndof*mndof) = mmi_init(lel2gel(el,pe),1:mndof*mndof)
          
          ael2gel(elcnt) = lel2gel(el,pe)
          gel2ael(lel2gel(el,pe)) = elcnt
          
!           PRINT*, elcnt,ael2gel(elcnt),nelnds(ael2gel(elcnt))
        ENDDO
      ENDDO
      
      
      ALLOCATE(Hwrite(ne,mndof),Qxwrite(ne,mndof),Qywrite(ne,mndof))
      
      DO el = 1,ne
        DO dof = 1,mndof
          Hwrite(el,dof)%ptr => H(gel2ael(el),dof)
          Qxwrite(el,dof)%ptr => Qx(gel2ael(el),dof)
          Qywrite(el,dof)%ptr => Qy(gel2ael(el),dof)
        ENDDO
      ENDDO
      
      
      ALLOCATE(edflag(nied))
      ALLOCATE(nparted(npart+1))
      
      edflag = 0
      edcnt = 0
      nparted = 0
      DO pe = 1,npart
        DO ed = 1,nied
        
          ged = iedn(ed)
          el_in = ged2el(1,ged)
          el_ex = ged2el(2,ged)
          
          pe_in = gel2part(el_in) 
          pe_ex = gel2part(el_ex)
          
          IF (pe_in == pe .and. pe_ex == pe) THEN
            
            edflag(ed) = edflag(ed) + 1
            edcnt = edcnt + 1
            
            nparted(pe) = nparted(pe) + 1
            
            CALL point_to_el(edcnt,ged,el_in,el_ex)
            
          ENDIF
        
        ENDDO
      ENDDO
      
      DO ed = 1,nied
        IF (edflag(ed) == 0) THEN
          edflag(ed) = edflag(ed) + 1
          edcnt = edcnt + 1
          
          ged = iedn(ed)
          el_in = ged2el(1,ged)
          el_ex = ged2el(2,ged)
          
          nparted(npart+1) = nparted(npart+1) + 1
          
          CALL point_to_el(edcnt,ged,el_in,el_ex)          
          
        ENDIF
      ENDDO
      
      
      
      
      ! check for errors
      DO ed = 1,nied
        IF(edflag(ed) /= 1) THEN
          PRINT*, "Edge partition error: edflag /= 1"
          STOP
        ENDIF
      ENDDO    
      
      ted = 0
      DO pe = 1,npart+1
        ted = ted + nparted(pe)
      ENDDO
      
      IF(ted /= nied .or. edcnt /= nied) THEN
        PRINT*, "Error: sum of partition interior edges /= total interior edges"
        PRINT*, "ted = ", ted
        PRINT*, "nied = ", nied
        PRINT*, "edcnt = ", edcnt
        STOP
      ENDIF
      

      ! Calculate blocking arrays
      ALLOCATE(edblk(2,npart))
      ALLOCATE(nfblk(2,npart+1))
      ALLOCATE(elblk(2,npart,nel_type))
      ALLOCATE(rnfblk(2,nrblk))
      
!       DO blk = 1,nblk
!         elblk(1,blk) = (blk-1)*(ne/nblk) + 1
!         elblk(2,blk) = blk*(ne/nblk)
!       ENDDO
!       elblk(2,nblk) = ne
      

      
!       tel = 0
!       edblk(1,1) = 1
!       DO blk = 1,npart-1
!         tel = tel + npartel(blk)
!         edblk(2,blk) = tel
!         edblk(1,blk+1) = tel + 1      
!       ENDDO
!       edblk(2,npart) = tel + npartel(npart)
      
      tel = 0 
      DO blk = 1,npart
        elblk(1,blk,1) = tel + 1
        DO et = 1,nel_type-1
          tel = tel + npartet(et,blk)
          elblk(2,blk,et) = tel
          elblk(1,blk,et+1) = tel+1
        ENDDO
        elblk(2,blk,nel_type) = tel + npartet(nel_type,blk)        
      ENDDO


      
      ted = 0 
      nfblk(1,1) = 1
      DO blk = 1,npart
        ted = ted + nparted(blk)
        nfblk(2,blk) = ted
        nfblk(1,blk+1) = ted + 1
      ENDDO
      nfblk(2,npart+1) = ted + nparted(npart+1)
      
      rnfblk(1,1) = ted+1
      DO blk = 1,nrblk-1
        ted = ted + (nparted(npart+1)/nrblk)
        rnfblk(2,blk) = ted  
        rnfblk(1,blk+1) = ted + 1
      ENDDO
      rnfblk(2,nrblk) = nied
      
      
      
      
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", "           Loop Blocking Information         "
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", " "            
      
!       PRINT*, "elblk: "
!       DO blk = 1,nblk
!         PRINT*, elblk(1,blk),elblk(2,blk) 
!       ENDDO
!       PRINT*, " "      
      
      PRINT*, "Number of partitions: ", npart      
      PRINT*, " "
      PRINT*, "edblk: "
      DO blk = 1,npart
        PRINT*, elblk(1,blk,1), elblk(2,blk,1), npartet(1,blk)
        PRINT*, elblk(1,blk,2), elblk(2,blk,2), npartet(2,blk)
        PRINT*, elblk(1,blk,3), elblk(2,blk,3), npartet(3,blk)
        PRINT*, elblk(1,blk,4), elblk(2,blk,4), npartet(4,blk)    
        PRINT*, " "
      ENDDO
      mnpartel = MAXVAL(npartel)
      PRINT*, "Max elements per partition: ", mnpartel      
      
      PRINT*, " "
      PRINT*, "nfblk: "
      DO blk = 1,npart+1
        PRINT*, nfblk(1,blk), nfblk(2,blk), nparted(blk)
      ENDDO
      mnparted = MAXVAL(nparted(1:npart))
      PRINT*, "Max edges per partition: ", mnparted
      PRINT*, " " 
      
      PRINT*, " "
      PRINT*, "nfblk: "
      DO blk = 1,nrblk
        PRINT*, rnfblk(1,blk), rnfblk(2,blk)
      ENDDO      
      
      END SUBROUTINE decomp2
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
      SUBROUTINE point_to_el(ed,ged,el_in,el_ex)
      
      USE globals, ONLY: nqpte,ged2led,gel2ael,gel2part, &
                         Hqpt,Qxqpt,Qyqpt, &
                         xmom,ymom,xymom, &
                         Hflux,Qxflux,Qyflux, &
                         Hi,He,Qxi,Qxe,Qyi,Qye, &
                         xmi,xme,ymi,yme,xymi,xyme, &
                         Hfi,Hfe,Qxfi,Qxfe,Qyfi,Qyfe, &
                         nx_pt,ny_pt,detJe, &
                         inx,iny,detJe_in,detJe_ex
                      
      
      IMPLICIT NONE
      
      INTEGER :: pt
      INTEGER :: ged
      INTEGER :: ed
      INTEGER :: el_in,el_ex
      INTEGER :: led_in,led_ex
      INTEGER :: gp_in,gp_ex
      INTEGER :: ael_in,ael_ex
      
      ael_in = gel2ael(el_in)
      ael_ex = gel2ael(el_ex)
      
!       PRINT "(I8,8x,3(I8),8X,3(I8))", ed,el_in,ael_in,gel2part(el_in), el_ex,ael_ex,gel2part(el_ex)
      
      DO pt = 1,nqpte(1)
      
        led_in = ged2led(1,ged)
        led_ex = ged2led(2,ged)

        gp_in = (led_in-1)*nqpte(1) + pt
        gp_ex = (led_ex-1)*nqpte(1) + nqpte(1) - pt + 1

        Hi(ed,pt)%ptr => Hqpt(ael_in,gp_in)
        He(ed,pt)%ptr => Hqpt(ael_ex,gp_ex)

        Qxi(ed,pt)%ptr => Qxqpt(ael_in,gp_in)
        Qxe(ed,pt)%ptr => Qxqpt(ael_ex,gp_ex)

        Qyi(ed,pt)%ptr => Qyqpt(ael_in,gp_in)
        Qye(ed,pt)%ptr => Qyqpt(ael_ex,gp_ex)

        xmi(ed,pt)%ptr => xmom(ael_in,gp_in)
        xme(ed,pt)%ptr => xmom(ael_ex,gp_ex)

        ymi(ed,pt)%ptr => ymom(ael_in,gp_in)
        yme(ed,pt)%ptr => ymom(ael_ex,gp_ex)

        xymi(ed,pt)%ptr => xymom(ael_in,gp_in)
        xyme(ed,pt)%ptr => xymom(ael_ex,gp_ex)

        Hfi(ed,pt)%ptr => Hflux(ael_in,gp_in)
        Hfe(ed,pt)%ptr => Hflux(ael_ex,gp_ex)

        Qxfi(ed,pt)%ptr => Qxflux(ael_in,gp_in)
        Qxfe(ed,pt)%ptr => Qxflux(ael_ex,gp_ex)

        Qyfi(ed,pt)%ptr => Qyflux(ael_in,gp_in)
        Qyfe(ed,pt)%ptr => Qyflux(ael_ex,gp_ex)
        
!         gp_in = pt
!         gp_ex = nqpte(1) - pt + 1        
        
        inx(ed,pt) = nx_pt(ged,pt)
        iny(ed,pt) = ny_pt(ged,pt)
          
        detJe_in(ed,pt) = detJe(ged,pt)
        detJe_ex(ed,pt) = detJe(ged,pt)        
     
      
      ENDDO
      
!       ! straight sided edges => constant normals
!       inx(ed) = nx_pt(ged,1)
!       iny(ed) = ny_pt(ged,1)
!           
!       ! straight sided edges => constant Jacobian
!       detJe_in(ed) = detJe(ged,1)
!       detJe_ex(ed) = detJe(ged,1)
      
      RETURN
      END SUBROUTINE point_to_el
      
