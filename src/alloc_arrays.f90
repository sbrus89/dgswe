      SUBROUTINE alloc_arrays()

      USE globals, ONLY: ne,nied,nfbed,nnfbed,nobed,ndof,nqpte, &
                         H,Hold,Hinit,rhsH,Qx,Qxold,Qxinit,rhsQx,Qy,Qyold,Qyinit, &
                         rhsQy,Hqpt,Qxqpt,Qyqpt,xmom,ymom,xymom, &
                         Hflux,Qxflux,Qyflux,tau,src_x,src_y,pressa,recipHa, &
                         rHi,rHe,xmomi,xmome,ymomi,ymome,xymomi,xymome, &
                         fbHf,nfbHf,obHf,fbQxf,nfbQxf,obQxf,fbQyf,nfbQyf,obQyf, &
                         Hfluxi,Hfluxe,Qxfluxi,Qxfluxe,Qyfluxi,Qyfluxe
                         

      IMPLICIT NONE
      INTEGER :: alloc_status
      
      
      ! Solution arrays
      
      ALLOCATE(Hinit(ne,ndof),Qxinit(ne,ndof),Qyinit(ne,ndof),STAT = alloc_status)
      IF(alloc_status /= 0) THEN
        PRINT*, 'Allocation error: Hinit,Qxinit,Qyinit'
      ENDIF
      
      ALLOCATE(H(ne,ndof),Qx(ne,ndof),Qy(ne,ndof),STAT = alloc_status)
      IF(alloc_status /= 0) THEN
        PRINT*, 'Allocation error: H,Qx,Qy'
      ENDIF 

      ! Old solution arrays
      ALLOCATE(Hold(ne,ndof),Qxold(ne,ndof),Qyold(ne,ndof),STAT = alloc_status)
      IF(alloc_status /= 0) THEN
        PRINT*, 'Allocation error: Hold,Qxold,Qyold'
      ENDIF 

      ! RHS arrays
      ALLOCATE(rhsH(ne,ndof),rhsQx(ne,ndof),rhsQy(ne,ndof),STAT = alloc_status)
      IF(alloc_status /= 0) THEN
        PRINT*, 'Allocation error: rhsH,rhsQx,rhsQy'
      ENDIF 

      ALLOCATE(Hqpt(ne,3*nqpte),Qxqpt(ne,3*nqpte),Qyqpt(ne,3*nqpte),STAT = alloc_status)
      IF(alloc_status /= 0) THEN
        PRINT*, 'Allocation error: Hqpt,Qxqpt,Qyqpt'
      ENDIF 

      ALLOCATE(xmom(ne,3*nqpte),ymom(ne,3*nqpte),xymom(ne,3*nqpte),STAT = alloc_status)
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

      ALLOCATE(Hflux(ne,3*nqpte),Qxflux(ne,3*nqpte),Qyflux(ne,3*nqpte),STAT = alloc_status)
      IF(alloc_status /= 0) THEN
        PRINT*, 'Allocation error: Hflux,Qxflux,Qyflux'
      ENDIF 

      ALLOCATE(pressa(ne),recipHa(ne))
      IF(alloc_status /= 0) THEN
        PRINT*, 'Allocation error: pressa,recipHa'
      ENDIF
      
      ALLOCATE(rHi(nied),rHe(nied),STAT = alloc_status)
      IF(alloc_status /= 0) THEN
        PRINT*, 'Allocation error: rHi,rHe'
      ENDIF      
      
      ALLOCATE(xmomi(nied),xmome(nied),ymomi(nied),ymome(nied),xymomi(nied),xymome(nied),STAT = alloc_status)
      IF(alloc_status /= 0) THEN
        PRINT*, 'Allocation error: xmomi,xmome,ymomi,ymome,xymomi,xymome'
      ENDIF            

      ALLOCATE(fbHf(nfbed,nqpte),nfbHf(nnfbed,nqpte),obHf(nobed,nqpte),STAT = alloc_status)
      IF(alloc_status /= 0) THEN
        PRINT*, 'Allocation error: fbHf,nfbHf,obHf'
      ENDIF
      
      ALLOCATE(fbQxf(nfbed,nqpte),nfbQxf(nnfbed,nqpte),obQxf(nobed,nqpte),STAT = alloc_status)
      IF(alloc_status /= 0) THEN
        PRINT*, 'Allocation error: fbQxf,nfbQxf,obQxf'
      ENDIF        
      
      ALLOCATE(fbQyf(nfbed,nqpte),nfbQyf(nnfbed,nqpte),obQyf(nobed,nqpte),STAT = alloc_status)
      IF(alloc_status /= 0) THEN
        PRINT*, 'Allocation error: fbQyf,nfbQyf,obQyf'
      ENDIF    
      
      ALLOCATE(Hfluxi(nied,nqpte),Qxfluxi(nied,nqpte),Qyfluxi(nied,nqpte),STAT = alloc_status)
      IF(alloc_status /= 0) THEN
        PRINT*, 'Allocation error: Hfluxi,Qxfluxi,Qyfluxi'
      ENDIF  
      
      ALLOCATE(Hfluxe(nied,nqpte),Qxfluxe(nied,nqpte),Qyfluxe(nied,nqpte),STAT = alloc_status)
      IF(alloc_status /= 0) THEN
        PRINT*, 'Allocation error: Hfluxe,Qxfluxe,Qyfluxe'
      ENDIF                
      
      RETURN
      END SUBROUTINE alloc_arrays
