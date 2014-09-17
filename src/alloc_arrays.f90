      SUBROUTINE alloc_arrays()

      USE globals, ONLY: ne,nied,nfbed,nnfbed,nobed,mndof,mnqpte, &
                         H,Hold,Hinit,rhsH,Qx,Qxold,Qxinit,rhsQx,Qy,Qyold,Qyinit, &
                         rhsQy,Hqpt,Qxqpt,Qyqpt,xmom,ymom,xymom, &
                         Hflux,Qxflux,Qyflux,tau,src_x,src_y,pressa,recipHa, &
                         rHi,rHe,xmomi,xmome,ymomi,ymome,xymomi,xymome, &
                         fbHf,nfbHf,obHf,fbQxf,nfbQxf,obQxf,fbQyf,nfbQyf,obQyf, &
                         Hfluxi,Hfluxe,Qxfluxi,Qxfluxe,Qyfluxi,Qyfluxe, &
                         MirhsH,MirhsQx,MirhsQy
                         

      IMPLICIT NONE
      INTEGER :: alloc_status
    
      
      ! Solution arrays
      
      ALLOCATE(Hinit(ne,mndof),Qxinit(ne,mndof),Qyinit(ne,mndof),STAT = alloc_status)
      IF(alloc_status /= 0) THEN
        PRINT*, 'Allocation error: Hinit,Qxinit,Qyinit'
      ENDIF
      
      ALLOCATE(H(ne,mndof),Qx(ne,mndof),Qy(ne,mndof),STAT = alloc_status)
      IF(alloc_status /= 0) THEN
        PRINT*, 'Allocation error: H,Qx,Qy'
      ENDIF 

      ! Old solution arrays
      ALLOCATE(Hold(ne,mndof),Qxold(ne,mndof),Qyold(ne,mndof),STAT = alloc_status)
      IF(alloc_status /= 0) THEN
        PRINT*, 'Allocation error: Hold,Qxold,Qyold'
      ENDIF 

      ! RHS arrays
      ALLOCATE(rhsH(ne,mndof),rhsQx(ne,mndof),rhsQy(ne,mndof),STAT = alloc_status)
      IF(alloc_status /= 0) THEN
        PRINT*, 'Allocation error: rhsH,rhsQx,rhsQy'
      ENDIF
      
      ALLOCATE(MirhsH(ne,mndof),MirhsQx(ne,mndof), MirhsQy(ne,mndof),STAT = alloc_status)
      IF(alloc_status /= 0) THEN
        PRINT*, 'Allocation error: MirhsH,MirhsQx,MirhsQy'
      ENDIF       

      ALLOCATE(Hqpt(ne,4*mnqpte),Qxqpt(ne,4*mnqpte),Qyqpt(ne,4*mnqpte),STAT = alloc_status)
      IF(alloc_status /= 0) THEN
        PRINT*, 'Allocation error: Hqpt,Qxqpt,Qyqpt'
      ENDIF 

      ALLOCATE(xmom(ne,4*mnqpte),ymom(ne,4*mnqpte),xymom(ne,4*mnqpte),STAT = alloc_status)
      IF(alloc_status /= 0) THEN
        PRINT*, 'Allocation error: xmom,ymom,xymom'
      ENDIF 

      ! Bottom friction array
      ALLOCATE(tau(ne),STAT = alloc_status)
      IF(alloc_status /= 0) THEN
        PRINT*, 'Allocation error: tau'
      ENDIF 

      ALLOCATE(src_x(ne),src_y(ne),STAT = alloc_status)
      IF(alloc_status /= 0) THEN
        PRINT*, 'Allocation error: src_x,src_y'
      ENDIF 

      ALLOCATE(Hflux(ne,4*mnqpte),Qxflux(ne,4*mnqpte),Qyflux(ne,4*mnqpte),STAT = alloc_status)
      IF(alloc_status /= 0) THEN
        PRINT*, 'Allocation error: Hflux,Qxflux,Qyflux'
      ENDIF 

      ALLOCATE(pressa(ne),recipHa(ne))
      IF(alloc_status /= 0) THEN
        PRINT*, 'Allocation error: pressa,recipHa'
      ENDIF
      
!       ALLOCATE(rHi(nied),rHe(nied),STAT = alloc_status)
!       IF(alloc_status /= 0) THEN
!         PRINT*, 'Allocation error: rHi,rHe'
!       ENDIF      
!       
!       ALLOCATE(xmomi(nied),xmome(nied),ymomi(nied),ymome(nied),xymomi(nied),xymome(nied),STAT = alloc_status)
!       IF(alloc_status /= 0) THEN
!         PRINT*, 'Allocation error: xmomi,xmome,ymomi,ymome,xymomi,xymome'
!       ENDIF            

!       ALLOCATE(fbHf(nfbed,nqpte),nfbHf(nnfbed,nqpte),obHf(nobed,nqpte),STAT = alloc_status)
!       IF(alloc_status /= 0) THEN
!         PRINT*, 'Allocation error: fbHf,nfbHf,obHf'
!       ENDIF
!       
!       ALLOCATE(fbQxf(nfbed,nqpte),nfbQxf(nnfbed,nqpte),obQxf(nobed,nqpte),STAT = alloc_status)
!       IF(alloc_status /= 0) THEN
!         PRINT*, 'Allocation error: fbQxf,nfbQxf,obQxf'
!       ENDIF        
!       
!       ALLOCATE(fbQyf(nfbed,nqpte),nfbQyf(nnfbed,nqpte),obQyf(nobed,nqpte),STAT = alloc_status)
!       IF(alloc_status /= 0) THEN
!         PRINT*, 'Allocation error: fbQyf,nfbQyf,obQyf'
!       ENDIF    
      
!       ALLOCATE(Hfluxi(nied,nqpte),Qxfluxi(nied,nqpte),Qyfluxi(nied,nqpte),STAT = alloc_status)
!       IF(alloc_status /= 0) THEN
!         PRINT*, 'Allocation error: Hfluxi,Qxfluxi,Qyfluxi'
!       ENDIF  
!       
!       ALLOCATE(Hfluxe(nied,nqpte),Qxfluxe(nied,nqpte),Qyfluxe(nied,nqpte),STAT = alloc_status)
!       IF(alloc_status /= 0) THEN
!         PRINT*, 'Allocation error: Hfluxe,Qxfluxe,Qyfluxe'
!       ENDIF                
      
      RETURN
      END SUBROUTINE alloc_arrays
