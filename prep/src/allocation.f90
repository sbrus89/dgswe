      MODULE allocation
      
      USE globals, ONLY: nel_type,ne,nn, &
                         mnepn,mnnds, &
                         ned,nied,nobed,nnfbed,nfbed, &
                         nope,neta,nbou,nvel, &
                         nobfr,nfbfr   
                         
      USE read_dginp, ONLY: p,ctp,npart,lines,hbp                         

      
      IMPLICIT NONE
      
      CONTAINS     
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE sizes()
      
      USE globals, ONLY: ndof,nverts,np,nnds, &
                         mndof,mnp,mnnds,nlines, &
                         hbnds,mhbnds
                         
      
      IMPLICIT NONE
      
      ndof(1) = (p+1)*(p+2)/2
      ndof(2) = (p+1)**2
      ndof(3) = ndof(1)
      ndof(4) = ndof(2)      
      mndof = maxval(ndof)
      
      nverts(1) = 3
      nverts(2) = 4
      nverts(3) = 3
      nverts(4) = 4
      
      np(1) = 1
      np(2) = 1
      np(3) = ctp
      np(4) = ctp   
      mnp = maxval(np)+1

      nnds(1) = 3
      nnds(2) = 4
      nnds(3) = (ctp+1)*(ctp+2)/2
      nnds(4) = (ctp+1)*(ctp+1) 
      mnnds = maxval(nnds)    
      
      hbnds(1) = (hbp+1)*(hbp+2)/2
      hbnds(2) = (hbp+1)*(hbp+1)       
      hbnds(3) = hbnds(1)
      hbnds(4) = hbnds(2)
      mhbnds = maxval(hbnds)
      
      nlines = lines
      
      RETURN
      END SUBROUTINE sizes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE alloc_grid_arrays(stage)
      
      USE globals, ONLY: ect,vct,xy,depth,nelnds,el_type,elxy,elhb,nepn, &
                         obseg,obnds,fbseg,fbnds,mhbnds
      
      IMPLICIT NONE
      INTEGER :: stage
      INTEGER :: n,i
      INTEGER :: alloc_status(3)
      
      alloc_status(:) = 0      
      
      IF (stage == 1) THEN
        n = 3      
        ! Node information
        ALLOCATE(ect(mnnds,ne),vct(4,ne),xy(2,nn),depth(nn),nelnds(ne),el_type(ne),STAT = alloc_status(1))  
        ALLOCATE(elxy(mnnds,ne,2),elhb(mhbnds,ne), STAT = alloc_status(2))
        ALLOCATE(nepn(nn),STAT = alloc_status(3))
      ELSE IF (stage == 2) THEN
        n = 1
        ! Open boundary information
        ALLOCATE(obseg(nope),obnds(neta,nope), STAT = alloc_status(1)) 
      ELSE IF (stage == 3) THEN
        n = 1
        ! Flow boundary information
        ALLOCATE(fbseg(2,nbou),fbnds(nvel,nbou),STAT = alloc_status(1))
      ENDIF
      
      
      
      DO i = 1,n
        IF (alloc_status(i) /= 0) THEN
          PRINT*, "Allocation error: alloc_grid_arrays"
          PRINT*, "Stage = ", stage
          STOP
        ENDIF
      ENDDO        
      
      
      RETURN 
      END SUBROUTINE alloc_grid_arrays
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE alloc_forcing_arrays(stage)
      
      USE globals, ONLY: obtag,obfreq,obnfact,obeq,obtag2,obamp,obph, &
                         fbtag,fbfreq,fbnfact,fbeq,fbtag2,fbamp,fbph
      
      IMPLICIT NONE
      INTEGER :: stage
      INTEGER :: n,i
      INTEGER :: alloc_status(2)
      
      alloc_status(:) = 0      
      
      IF (stage == 1) THEN
        n = 2
        ! Open boundary forcing arrays
        ALLOCATE(obtag(nobfr),obfreq(nobfr),obnfact(nobfr),obeq(nobfr),STAT = alloc_status(1))
        ALLOCATE(obtag2(nobfr),obamp(neta,nobfr),obph(neta,nobfr),STAT = alloc_status(2))
      ELSE IF (stage == 2) THEN
        n = 2
        ! Flow boundary forcing arrays
        ALLOCATE(fbtag(nfbfr),fbfreq(nfbfr),fbnfact(nfbfr),fbeq(nfbfr),STAT = alloc_status(1))
        ALLOCATE(fbtag2(nfbfr),fbamp(nvel,nbou,nfbfr),fbph(nvel,nbou,nfbfr),STAT = alloc_status(2))
      ENDIF
      
      
      
      DO i = 1,n
        IF (alloc_status(i) /= 0) THEN
          PRINT*, "Allocation error: alloc_forcing_arrays"
          PRINT*, "Stage = ", stage
          STOP
        ENDIF
      ENDDO        
      
      
      RETURN 
      END SUBROUTINE alloc_forcing_arrays      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE alloc_connect_arrays(stage)
      
      USE globals, ONLY: epn,ged2nn,ged2el,ged2led,&
                         iedn,obedn,fbedn,nfbedn,bedn
      
      IMPLICIT NONE
      INTEGER :: stage
      INTEGER :: n,i
      INTEGER :: alloc_status(3)
      
      alloc_status(:) = 0      
      
      IF (stage == 1) THEN
        n = 1        
        ! Elements associated with each node
        ALLOCATE(epn(mnepn,nn),STAT = alloc_status(1))
      ELSE IF (stage == 2) THEN
        n = 2
        ! Edge look-up tables
        ALLOCATE(ged2nn(2,ned),ged2el(2,ned),ged2led(2,ned),STAT = alloc_status(1))
        ALLOCATE(bedn(ned),STAT = alloc_status(2))
      ELSE IF (stage == 3) THEN
        n = 3
        ! Edge look-up tables
        ALLOCATE(iedn(nied),STAT = alloc_status(1))
        ALLOCATE(obedn(nobed),STAT = alloc_status(2))
        ALLOCATE(fbedn(nfbed),nfbedn(nnfbed),STAT=alloc_status(3))
      ENDIF
      
      
      
      DO i = 1,n
        IF (alloc_status(i) /= 0) THEN
          PRINT*, "Allocation error: alloc_connect_arrays"
          PRINT*, "Stage = ", stage
          STOP
        ENDIF
      ENDDO          
      
      RETURN
      END SUBROUTINE alloc_connect_arrays

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      END MODULE allocation
