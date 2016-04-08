      MODULE curvilinear_nodes_mod

      
      USE globals, ONLY: rp,pi
      
      
      IMPLICIT NONE

      CONTAINS
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  


      SUBROUTINE shape_functions_linear_at_ctp(nel_type,np,psi)

      USE basis, ONLY: element_nodes
      USE shape_functions_mod, ONLY: shape_functions_area_eval      
      
      IMPLICIT NONE      
      
      INTEGER, INTENT(IN) :: nel_type
      INTEGER, DIMENSION(:), INTENT(IN) :: np
      REAL(rp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(OUT) :: psi
      
      INTEGER :: pt,i,j
      INTEGER :: et,str_typ,crv_typ
      INTEGER :: npts,nnd 
      INTEGER :: mnp,mnnds
      REAL(rp), DIMENSION(:), ALLOCATABLE :: r,s
       

      ! Evaluates linear shape functions at curved element nodal sets
      !
      ! Used to create additional nodes which can then be adjusted to make 
      ! straight elements curved. See edge_coordinates_curved()
      
      mnp = maxval(np)
      mnnds = (mnp+1)**2
      ALLOCATE(r(mnnds),s(mnnds))
      
      ALLOCATE(psi(mnnds,mnnds,2))             
      psi = 0d0 
      
      DO et = 1,2

        IF (mod(et,2) == 1) THEN
          str_typ = 1
          crv_typ = 3
        ELSE IF (mod(et,2) == 0) THEN
          str_typ = 2
          crv_typ = 4
        ENDIF  

        CALL element_nodes(et,1,np(crv_typ),npts,r,s)
        CALL shape_functions_area_eval(et,np(str_typ),nnd,npts,r,s,psi(:,:,et))   
              
!         
!         PRINT*, "psiv"
!         DO i = 1,nnd
!           PRINT "(300(f27.17))", (psi(i,j,et), j = 1,npts)
!         ENDDO
!         PRINT*," "
!         

      ENDDO
      
      END SUBROUTINE shape_functions_linear_at_ctp

      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  


      SUBROUTINE eval_coordinates_curved(ctp,nnds,nverts,el_type,xy,ect,fbseg,fbnds, &
                                         nnfbed,nfbedn,nfbednn,ged2el,ged2led, &
                                         psiv,bndxy,elxy)
     
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: ctp
      INTEGER, DIMENSION(:), INTENT(IN) :: nnds
      INTEGER, DIMENSION(:), INTENT(IN) :: nverts
      INTEGER, DIMENSION(:), INTENT(INOUT) :: el_type
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: xy
      INTEGER, DIMENSION(:,:), INTENT(IN) :: ect
      INTEGER, DIMENSION(:,:), INTENT(IN) :: fbseg
      INTEGER, DIMENSION(:,:), INTENT(IN) :: fbnds
      INTEGER, INTENT(IN) :: nnfbed
      INTEGER, DIMENSION(:), INTENT(IN) :: nfbedn
      INTEGER, DIMENSION(:,:), INTENT(IN) :: nfbednn
      INTEGER, DIMENSION(:,:), INTENT(IN) :: ged2el
      INTEGER, DIMENSION(:,:), INTENT(IN) :: ged2led
      REAL(rp), DIMENSION(:,:,:), INTENT(IN) :: psiv
      REAL(rp), DIMENSION(:,:,:,:), INTENT(IN) :: bndxy
      REAL(rp), DIMENSION(:,:,:), INTENT(INOUT) :: elxy
     
      INTEGER :: i,ed,nd
      INTEGER :: ged,seg,n1,el,led
      INTEGER :: n1ind
      REAL(rp), DIMENSION(2,ctp-1) :: segxy     
      
      IF (ctp == 1) THEN
        RETURN
      ENDIF
      
      
      
      
      DO ed = 1,nnfbed
        ged = nfbedn(ed)
        
        seg = nfbednn(ed,1)
        n1 = nfbednn(ed,2)
        
        el = ged2el(1,ged)
        led = ged2led(1,ged)     
        
        n1ind = 0
 search:DO i = 1,fbseg(1,seg)-1
          IF (fbnds(i,seg) == n1) THEN
            n1ind = i
            EXIT search
          ENDIF
        ENDDO search
        
        IF (n1ind == 0) THEN
          PRINT*, "node number not found"
          STOP
        ENDIF
        
        DO nd = 1,ctp-1
          segxy(1,nd) = bndxy(1,nd,n1ind,seg)
          segxy(2,nd) = bndxy(2,nd,n1ind,seg)
        ENDDO
        
        CALL edge_coordinates_curved(el,ctp,led,nnds,nverts,el_type,xy,ect,segxy,psiv,elxy)
                
        
      ENDDO     
     
     END SUBROUTINE eval_coordinates_curved


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE edge_coordinates_curved(el,ctp,led,nnds,nverts,el_type,xy,ect,segxy,psiv,elxy)
      
      USE transformation, ONLY: element_transformation            
       
      IMPLICIT NONE
       
      INTEGER, INTENT(IN) :: el
      INTEGER, INTENT(IN) :: ctp
      INTEGER, INTENT(IN) :: led      
      INTEGER, DIMENSION(:), INTENT(IN) :: nnds
      INTEGER, DIMENSION(:), INTENT(IN) :: nverts
      INTEGER, DIMENSION(:), INTENT(INOUT) :: el_type
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: xy
      INTEGER, DIMENSION(:,:), INTENT(IN) :: ect
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: segxy
      REAL(rp), DIMENSION(:,:,:), INTENT(IN) :: psiv
      REAL(rp), DIMENSION(:,:,:), INTENT(INOUT) :: elxy
      
      INTEGER :: pt,nd
      INTEGER :: nnd,nv,et,mnnds
      REAL(rp) :: xpt,ypt
      REAL(rp), DIMENSION(:), ALLOCATABLE :: x,y
       
      mnnds = MAXVAL(nnds)
      ALLOCATE(x(mnnds),y(mnnds))
       
       
      et = el_type(el)
      nv = nverts(et)  
        
      IF (et <= 2) THEN            ! only compute extra nodes if element isn't already curved
                                   ! otherwise, just adjust the boundary nodes. 
        IF (mod(et,2) == 1) THEN   ! this is important for elements with two no normal flow boundaries
          el_type(el) = 3 
          nnd = nnds(3)      
        ELSE IF (mod(et,2) == 0) THEN
          el_type(el) = 4
          nnd = nnds(4)                  
        ENDIF         
        
        DO nd = 1,nv
          x(nd) = xy(1,ect(nd,el))
          y(nd) = xy(2,ect(nd,el))
        ENDDO        
        
        DO pt = 1,nnd               

          CALL element_transformation(nv,x,y,psiv(:,pt,et),xpt,ypt)
        
          elxy(pt,el,1) = xpt
          elxy(pt,el,2) = ypt
        ENDDO      
        
      ENDIF
        
      DO nd = 1,ctp-1
        pt = mod(led,nv)*ctp + 1 + nd
          
!         ytest = elxy(pt,el,2)
!         xpt = elxy(pt,el,1) 
!           
!         IF (ytest < 250d0) THEN
!           ypt = 0d0 + 100d0*(1d0/(COSH(4d0*(xpt-2000d0)/500d0)))
!         ELSE IF (ytest > 250d0) THEN
!           ypt = 500d0 - 100d0*(1d0/(COSH(4d0*(xpt-2000d0)/500d0)))
!         ENDIF
!         
!         elxy(pt,el,2) = ypt

        elxy(pt,el,1) = segxy(1,nd)
        elxy(pt,el,2) = segxy(2,nd)
      ENDDO           
       
      RETURN
      END SUBROUTINE edge_coordinates_curved


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
      
!       SUBROUTINE blend_coordinates()
!       
!       
!       IF ( BLENDFLAG == 1 ) THEN
! !C============================================ DW =========================== 
! !     ----- For blending ---  
!         SELECT CASE( LED)
!         CASE (3)
!             Vr = CLNDS(1,:) ; ! r
!         CASE (1)
!             Vr = CLNDS(2,:) ; ! r
!         CASE (2)
!             Vr = CLNDS(2,:) ; ! s
!         END SELECT
! 
!         !IF ( LED == 3 ) Vr = CLNDS(1,:) ; ! r
!         !IF ( LED == 1 ) Vr = CLNDS(2,:) ; ! s
!         !IF ( LED == 2 ) Vr = CLNDS(2,:) ; ! s
!  
!         Fr = Vr(Fmask(:,LED)) ; 
!         CALL Vandermonde1Da( Vface, CLP, Fr ) ; 
!         CALL Vandermonde1Da( Vvol, CLP, Vr  ) ;  
!       
!         !c LU   
!         NI = CLP + 1 ;
!         NRHS = 1 ;
!         CALL DGETRF(NI,NI,Vface,NI,IPIV,IINFO) ;
!  
!         FDXM(:,1) = FDX ;
!         CALL DGETRS( 'N', NI, NRHS, Vface, NI, IPIV, FDXM, NI, IINFO ) ;        
!         Vdx = MATMUL(VVol, FDXM(:,1) ) ;
!  
!         FDYM(:,1) = FDY ;
!         CALL DGETRS( 'N', NI, NRHS, Vface, NI, IPIV, FDYM, NI, IINFO ) ;        
!         Vdy = MATMUL(Vvol, FDYM(:,1) ) ;
! 
!       
!         ids = 0 ; 
!         nids = 0 ;
!         DO PT = 1, NCLNDS
!           IF ( abs(1.D0 - Vr(PT)) > 1.0D-8 ) THEN
!              nids = nids + 1 ;
! 
!              ids(nids) = PT ;
!           END IF
!         END DO  
!     
!         blend = 0.D0 ; 
!         SELECT CASE (LED)
!         CASE (3)
!            blend(ids(1:nids)) = -(CLNDS(1,ids(1:nids)) + CLNDS(2,ids(1:nids)))/ &
!         &   (1.D0 - Vr(ids(1:nids))) ; 
!         CASE (1) 
!            blend(ids(1:nids)) =  (CLNDS(1,ids(1:nids)) + 1.D0)/ &
!         &   (1.D0 - Vr(ids(1:nids))) ; 
!         CASE (2)
!            blend(ids(1:nids)) = -(CLNDS(1,ids(1:nids)) + CLNDS(2,ids(1:nids)))/ &
!         &   (1.D0 - Vr(ids(1:nids))) ; 
!         END SELECT 
! 
!         XCL(:,EDCNT) = XCL(:,EDCNT) + blend*Vdx ;
!         YCL(:,EDCNT) = YCL(:,EDCNT) + blend*Vdy ;
!       ELSE
! 
!        XCL(MOD(LED,3)*CLP+2:MOD(LED,3)*CLP+CLP,EDCNT) = &
!             & XCL(MOD(LED,3)*CLP+2:MOD(LED,3)*CLP+CLP,EDCNT) + FDX(2:CLP) ;
!        YCL(MOD(LED,3)*CLP+2:MOD(LED,3)*CLP+CLP,EDCNT) = &
!             & YCL(MOD(LED,3)*CLP+2:MOD(LED,3)*CLP+CLP,EDCNT) + FDY(2:CLP) ; 
!       END IF
! !C=========================================== DW ============================  
! 
!       END SUBROUTINE blend_coordinates
      
      END MODULE curvilinear_nodes_mod