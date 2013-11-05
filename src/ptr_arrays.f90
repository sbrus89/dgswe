     SUBROUTINE ptr_arrays()
     
     USE globals, ONLY: ne,ned,nied,nqpte,iedn,ged2led,ged2el,edlen_area,normal, &
                       Hi,He,Qxi,Qxe,Qyi,Qye, &
                       xmi,xme,ymi,yme,xymi,xyme, &
                       Hfi,Hfe,Qxfi,Qxfe,Qyfi,Qyfe, &
                       Hqpt,Qxqpt,Qyqpt, &
                       xmom,ymom,xymom, &
                       Hflux,Qxflux,Qyflux, &
                       inx,iny,len_area_in,len_area_ex,const, &
                       Hhatv,Qxhatv,Qyhatv

     IMPLICIT NONE
     INTEGER :: alloc_status
     INTEGER :: ed,ged,pt
     INTEGER :: gp_in,gp_ex,led_in,led_ex,el_in,el_ex

     ALLOCATE(Hi(nied,nqpte),He(nied,nqpte),Qxi(nied,nqpte),Qxe(nied,nqpte),Qyi(nied,nqpte),Qye(nied,nqpte),STAT=alloc_status)
     IF(alloc_status /= 0) THEN
       PRINT*, "Allocation error: Hi,He,Qxi,Qxe,Qyi,Qye"
     ENDIF

     ALLOCATE(xmi(nied,nqpte),xme(nied,nqpte),ymi(nied,nqpte),yme(nied,nqpte),xymi(nied,nqpte),xyme(nied,nqpte),STAT=alloc_status)
     IF(alloc_status /= 0) THEN
       PRINT*, "Allocation error: xmi,xme,ymi,yme,xymi,xyme"
     ENDIF

     ALLOCATE(Hfi(nied,nqpte),Hfe(nied,nqpte),Qxfi(nied,nqpte),Qxfe(nied,nqpte),Qyfi(nied,nqpte),Qyfe(nied,nqpte),STAT=alloc_status)
     IF(alloc_status /= 0) THEN
       PRINT*, "Allocation error: Hfi,Hfe,Qxfi,Qxfe,Qyfi,Qyfe"
     ENDIF

     ALLOCATE(const(nied),inx(nied),iny(nied),len_area_in(nied),len_area_ex(nied),STAT=alloc_status)
     IF(alloc_status /= 0) THEN
       PRINT*, "Allocation error: const,inx,iny,len_area_in,len_area_ex"
     ENDIF

     ALLOCATE(Hhatv(nied),Qxhatv(nied),Qyhatv(nied),STAT=alloc_status)
     IF(alloc_status /= 0) THEN
       PRINT*, "Allocation error: Hhatv,Qxhatv,Qyhatv"
     ENDIF


!      DO pt = 1,nqpte
!        DO ed = 1,nied
!        
!          ged = iedn(ed)
! 
!          led_in = ged2led(1,ged)
!          led_ex = ged2led(2,ged)
! 
!          gp_in = (led_in-1)*nqpte + pt
!          gp_ex = (led_ex-1)*nqpte + nqpte - pt + 1
!   
!          el_in = ged2el(1,ged)
!          el_ex = ged2el(2,ged)
! 
!          Hi(ed,pt)%ptr => Hqpt(el_in,gp_in)
!          He(ed,pt)%ptr => Hqpt(el_ex,gp_ex)
! 
!          Qxi(ed,pt)%ptr => Qxqpt(el_in,gp_in)
!          Qxe(ed,pt)%ptr => Qxqpt(el_ex,gp_ex)
! 
!          Qyi(ed,pt)%ptr => Qyqpt(el_in,gp_in)
!          Qye(ed,pt)%ptr => Qyqpt(el_ex,gp_ex)
! 
!          xmi(ed,pt)%ptr => xmom(el_in,gp_in)
!          xme(ed,pt)%ptr => xmom(el_ex,gp_ex)
! 
!          ymi(ed,pt)%ptr => ymom(el_in,gp_in)
!          yme(ed,pt)%ptr => ymom(el_ex,gp_ex)
! 
!          xymi(ed,pt)%ptr => xymom(el_in,gp_in)
!          xyme(ed,pt)%ptr => xymom(el_ex,gp_ex)
! 
!          Hfi(ed,pt)%ptr => Hflux(el_in,gp_in)
!          Hfe(ed,pt)%ptr => Hflux(el_ex,gp_ex)
! 
!          Qxfi(ed,pt)%ptr => Qxflux(el_in,gp_in)
!          Qxfe(ed,pt)%ptr => Qxflux(el_ex,gp_ex)
! 
!          Qyfi(ed,pt)%ptr => Qyflux(el_in,gp_in)
!          Qyfe(ed,pt)%ptr => Qyflux(el_ex,gp_ex)
! 
!        ENDDO
!      ENDDO
! 
!      DO ed = 1,nied
! 
!       ged = iedn(ed)
! 
!       inx(ed) = normal(1,ged)
!       iny(ed) = normal(2,ged)
! 
!       len_area_in(ed) = edlen_area(1,ged)
!       len_area_ex(ed) = edlen_area(2,ged)   
! 
!      ENDDO
! 
!      DO ged = 1,ned
! 
!        el_in = ged2el(1,ged)
!        el_ex = ged2el(2,ged)
! 
!        led_in = ged2led(1,ged)
!        led_ex = ged2led(2,ged)
! 
!        nxv(el_in,led_in) = normal(1,ged)
!        nyv(el_in,led_in) = normal(2,ged)
! 
!        IF(el_ex.gt.0 .and. led_ex.gt.0) THEN
!          nxv(el_ex,led_ex) = normal(1,ged)
!          nyv(el_ex,led_ex) = normal(2,ged)
!        ENDIF
!        
!      ENDDO
          

     RETURN
     END SUBROUTINE ptr_arrays
