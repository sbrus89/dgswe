      SUBROUTINE read_partitions()
      
      USE globals, ONLY: npart,tnpel,mnelpp,peln,nprel,preln
      
      IMPLICIT NONE
      
      INTEGER :: part,nd,el,i
      INTEGER :: nelg,nnodg
      INTEGER :: mnnpp
      INTEGER :: ignore1,ignore2
      INTEGER :: pn,nnp,nres
      INTEGER :: neighp_r,neighp_s
      INTEGER :: nrel,nsel
      INTEGER :: el_cnt
      INTEGER :: alloc_status
      
      INTEGER, ALLOCATABLE, DIMENSION(:) :: npresel
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: preseln
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: pnodes
      INTEGER, ALLOCATABLE, DIMENSION(:) :: npsel
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: pseln
      
      OPEN(UNIT=80,FILE="fort.80")
      
      READ(80,"A80")
      READ(80,"A80")
      READ(80,"A80")
      
      READ(80,*) nelg,nnodg ! total number of elements and nodes
      READ(80,*) npart ! number of element partitions
      READ(80,*) mnnpp ! max. number of nodes per partition
      READ(80,*) mnelpp ! max. number of elements per partition
      
      READ(80,*) ignore1
      READ(80,*) ignore1
      READ(80,*) ignore1,ignore2
      
      
      ALLOCATE(pnodes(mnnpp,npart),STAT = alloc_status)
      IF(alloc_status /= 0) THEN
        PRINT*, "Allocation error: pnodes"
      ENDIF
      
      DO part = 1,npart
        READ(80,*) pn,nnp,nres
        READ(80,1130) (pnodes(nd,part), nd = 1,nnp)
      ENDDO
      
      ALLOCATE(tnpel(npart),peln(mnelpp,npart),STAT = alloc_status)
      IF(alloc_status /= 0) THEN
        PRINT*, "Allocation error: tnpel,peln"
      ENDIF
      
      DO part = 1,npart
        READ(80,*) pn,tnpel(part)
        IF (pn /= part-1) THEN
          PRINT*, "Error reading fort.80"
          STOP
        ENDDO
        READ(80,1130) (peln(el,part), el=1,tnpel(part))
      ENDDO
      
      OPEN(UNIT=18,FILE="DG.18")
      
      ALLOCATE(npresel(npart),preseln(mnelpp,npart),STAT = alloc_status)
      IF(alloc_status /= 0) THEN
        PRINT*, "Allocation error: npresel,preseln"
      ENDIF      
      
      ALLOCATE(nprel(npart),preln(mnelpp,npart),STAT = alloc_status)
      IF(alloc_status /= 0) THEN
        PRINT*, "Allocation error: nprel,preln"
      ENDIF  
      nprel = 0
      
      ALLOCATE(npsel(npart),pseln(mnelpp,npart),STAT = alloc_status)
      IF(alloc_status /= 0) THEN
        PRINT*, "Allocation error: npsel,pseln"
      ENDIF  
      npsel = 0      
      
      DO part = 1,npart
      
        READ(18,3010) pn,npresel(part)
        IF (np /= part-1) THEN
          PRINT*, "Error reading DG.18"
          STOP
        ENDIF
        
        READ(80,1130) (preseln(el,part), el=1,npresel(part))
        READ(80,3010) neighp_r,neighp_s
        
        el_cnt = 0
        DO i = 1,neighp_r
          READ(80,3010) pn,nrel
          nprel(npart) = nprel(npart) + nrel
          READ(80,1130) (preln(el+el_cnt,part), el = 1,nrel)
          el_cnt = el_cnt + nrel
        ENDDO
        
        el_cnt = 0
        DO i = 1,neighp_s
          READ(80,3010)pn,nsel
          npsel(npart) = npsel(npart) + nsel
          READ(80,1130) (pseln(el+el_cnt,part), el = 1,nsel)
          el_cnt = el_cnt + nsel
        ENDDO
        
      ENDDO
      
 1130 FORMAT(8X,9I8)  
 3010 FORMAT(8X,2I8)
      
      RETURN
      END SUBROUTINE read_partitions
      
      
      
      
      
      
      
      
      SUBROUTINE align_partitions

      USE globals, ONLY: ndof,npart,tnpel,npel,mnelpp, &
                         lel2gel,gel2ael,ael2gel, &
                         H,Hinit,Qx,Qxinit,Qy,Qyinit, &
                         Hwrite,Qxwrite,Qywrite
      
      IMPLICIT NONE
      
      INTEGER :: part,el
      INTEGER :: el_cnt
      INTEGER :: alloc_status
      
      ALLOCATE(lel2gel(mnelpp,npart),npel(npart),gel2ael(ne),ael2gel(ne),STAT = alloc_status)
      IF(alloc_status /= 0) THEN
        PRINT*, "Allocation error: npel,lel2gel,gel2ael,ael2gel"
      ENDIF       
      
      DO part = 1,npart
      
        el_cnt = 0
        DO el = 1,tnpel(part)
        
          IF (ANY(preln(:,part) == el) THEN
          
          ELSE
            el_cnt = el_cnt + 1
            lel2gel(el_cnt,part) = peln(el,part)
          ENDIF
          
        ENDDO
        npel(part) = el_cnt
      ENDDO
      
      el_cnt = 1
      DO part = 1,npart
        DO el = 1,npel(part)
          H(el_cnt,1:ndof) = Hinit(lel2gel(el,part),1:ndof)
          Qx(el_cnt,1:ndof) = Qxinit(lel2gel(el,part),1:ndof)
          Qy(el_cnt,1:ndof) = Qyinit(lel2gel(el,part),1:ndof)
          
          ael2gel(el_cnt) = (lel2gel(el,part))
          gel2ael(lel2gel(el,part)) = el_cnt
          
          el_cnt = el_cnt + 1
        ENDDO
      ENDDO
      
      ALLOCATE(Hwrite(ne,ndof),Qxwrite(ne,ndof),Qywrite(ne,ndof),STAT = alloc_status)
      IF(alloc_status /= 0) THEN
        PRINT*, "Allocation error: Hwrite,Qxwrite,Qywrite"
      ENDIF      
      
      DO el = 1,ne
        DO dof = 1,ndof
          Hwrite(el,dof)%ptr => H(ael2gel(el),dof)
          Qxwrite(el,dof)%ptr => Qx(ael2gel(el),dof)
          Qywrite(el,dof)%ptr => Qy(ael2gel(el),dof)          
        ENDDO
      ENDDO
      
      RETURN
      END SUBROUTINE align_partitions
      
      
      
      
      
      
      
      
      
      SUBROUTINE edge_partition()
      
      USE globals, ONLY: npart,nied, &
                         ged2el,lel2gel, &
                         inx,iny,len_area_in,len_area_ex, &
                         
      
      IMPLICIT NONE
      
      INTEGER :: part,ed
      INTEGER :: ged,el_in,el_ex
      INTEGER :: edcnt
      
      
      
      ALLOCATE(edflag(nied),nied_part(npart+1))
      
      edflag = 0
      nied_part = 0
      
      edcnt = 0
      
      ! loop through each partition and find interior edges associated with those elements
      DO part = 1,npart 

        DO ed = 1,nied

          ged = iedn(ed)
          el_in = ged2el(1,ged)
          el_ex = ged2el(2,ged)
            
          ! check if both elements are in the partition
          IF(ANY(lel2gel(:,part).eq.el_in) .and. ANY(lel2gel(:,part)).eq.el_ex) THEN
            
              edflag(ed) = 1
              edcnt = edcnt + 1
              
              nied_part(part) = nied_part(part) + 1
!               iedn_part(nied_part(part)) = ged
                
              CALL point_to_el(edcnt,ged,el_in,el_ex)
              
              inx(edcnt) = normal(1,ged)
              iny(edcnt) = normal(2,ged)
              
              len_area_in(edcnt) = edlen_area(1,ged)
              len_area_ex(edcnt) = edlen_area(2,ged)
                
          ELSE
            ! ignore edges that cotain elements from two different partitions (for now)
            edflag(ed) = 0
          ENDIF

        ENDDO
          
      ENDDO
      
      DO ed = 1,nied
        IF(edflag(ed) == 0) THEN
        
          edflag(ed) = 1
          edcnt = edcnt + 1
          
          nied_part(npart+1) = nied_part(npart+1) + 1
          
          ged = iedn(ed)
          el_in = ged2el(1,ged)
          el_ex = ged2el(2,ged)
          
          CALL point_to_el(edcnt,ged,el_in,el_ex)
          
          inx(edcnt) = normal(1,ged)
          iny(edcnt) = normal(2,ged)
          
          len_area_in(edcnt) = edlen_area(1,ged)
          len_area_ex(edcnt) = edlen_area(2,ged)
          
        ENDIF
      ENDDO
      
      RETURN
      END SUBROUTINE edge_partition
      
      
      
      
      
      
      
      
      SUBROUTINE point_to_el(ed,ged,el_in,el_ex)
      
      USE globals, ONLY: nqpte,ged2led,gel2ael, &
                         Hqpt,Qxqpt,Qyqpt, &
                         xmom,ymom,xymom, &
                         Hflux,Qxflux,Qyflux, &
                         Hi,He,Qxi,Qxe,Qyi,Qye, &
                         xmi,xme,ymi,yme,xymi,xyme, &
                         Hfi,Hfe,Qxfi,Qxfe,Qyfi,Qyfe
                      
      
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
      
      DO pt = 1,nqpte
      
        led_in = ged2led(1,ged)
        led_ex = ged2led(2,ged)

        gp_in = (led_in-1)*nqpte + pt
        gp_ex = (led_ex-1)*nqpte + nqpte - pt + 1

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

!         Hni(ed,pt)%ptr => Hn(ael_in,gp_in)
!         Hne(ed,pt)%ptr => Hn(ael_ex,gp_ex)
! 
!         Qxni(ed,pt)%ptr => Qxn(ael_in,gp_in)
!         Qxne(ed,pt)%ptr => Qxn(ael_ex,gp_ex)
! 
!         Qyni(ed,pt)%ptr => Qyn(ael_in,gp_in)
!         Qyne(ed,pt)%ptr => Qyn(ael_ex,gp_ex)
! 
!         EVi(ed,pt)%ptr => egnval(ael_in,gp_in)
!         EVe(ed,pt)%ptr => egnval(ael_ex,gp_ex)      
      
      ENDDO
      
      RETURN
      END SUBROUTINE point_to_el
