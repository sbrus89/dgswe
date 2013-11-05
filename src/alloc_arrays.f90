      SUBROUTINE alloc_arrays()

      USE globals, ONLY: sp,ne,nied,ndof,nqpte, &
                         H,Hold,Hinit,rhsH,Qx,Qxold,Qxinit,rhsQx,Qy,Qyold,Qyinit, &
                         rhsQy,Hqpt,Qxqpt,Qyqpt,xmom,ymom,xymom, &
                         Hflux,Qxflux,Qyflux,tau,src_x,src_y,pressa,recipHa, &
                         nsp,nsp2,split,split2
                         

      IMPLICIT NONE
      INTEGER :: alloc_status
      
      ALLOCATE(split(2,nsp))
      
      DO sp = 1,nsp
        split(1,sp) = (sp-1)*(ne/nsp) + 1
        split(2,sp) = sp*(ne/nsp)
      ENDDO
      split(2,nsp) = ne

      DO sp = 1,nsp
        PRINT*, split(1,sp),split(2,sp) 
      ENDDO
      PRINT*, ' '
      
      ALLOCATE(split2(2,nsp2))
      
      DO sp = 1,nsp2
        split2(1,sp) = (sp-1)*(nied/nsp2) + 1
        split2(2,sp) = sp*(nied/nsp2)
      ENDDO
      split2(2,nsp2) = nied

      DO sp = 1,nsp2
        PRINT*, split2(1,sp),split2(2,sp) 
      ENDDO
      PRINT*, ' '

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

      
      RETURN
      END SUBROUTINE alloc_arrays
