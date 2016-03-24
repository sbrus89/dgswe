      MODULE bathymetry_interp_mod

      USE globals, ONLY: rp
      
      IMPLICIT NONE
      
      CONTAINS

      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
       
      
      
      SUBROUTINE shape_functions_eltypes_at_hbp(nel_type,np,psi)

      USE basis, ONLY: element_nodes
      USE shape_functions_mod, ONLY: shape_functions_area_eval      
      
      IMPLICIT NONE    
      
      INTEGER, INTENT(IN) :: nel_type
      INTEGER, DIMENSION(:), INTENT(IN) :: np
      REAL(rp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(OUT) :: psi
      
      
      INTEGER :: pt,i,j
      INTEGER :: et,hbp_type
      INTEGER :: npts,nnd
      INTEGER :: mnp,mnnds
      REAL(rp), DIMENSION(:), ALLOCATABLE :: r,s
      
      ! Evaluates linear/curvilinear shape functions at high-order batymetry nodal sets
      !
      ! Used to compute function-specified bathymetry at high-order batymetry nodes
      ! for curved elements. See bathy_coordinates()    
       
      mnp = maxval(np)
      mnnds = (mnp+1)**2
      ALLOCATE(r(mnnds),s(mnnds))
      
      ALLOCATE(psi(mnnds,mnnds,nel_type))             
      psi = 0d0 
      
      
       
      DO et = 1,nel_type     

        IF (mod(et,2) == 1) THEN
          hbp_type = 5    
        ELSE IF (mod(et,2) == 0) THEN
          hbp_type = 6
        ENDIF   
                 

        CALL element_nodes(et,1,np(hbp_type),npts,r,s)
        CALL shape_functions_area_eval(et,np(et),nnd,npts,r,s,psi(:,:,et))             
                
!         PRINT*, "psic"        
!         DO i = 1,nnd
!           PRINT "(20(f27.17))", (psi(i,j,et), j = 1,npts)
!         ENDDO
!         PRINT*," "


      ENDDO     
      
      END SUBROUTINE shape_functions_eltypes_at_hbp      
      
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       



      SUBROUTINE bathy_coordinates(ne,nnds,nverts,el_type,elxy,psic,xyhb,depth,ect,elhb)
      
      USE transformation, ONLY: element_transformation      
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: ne   
      INTEGER, DIMENSION(:), INTENT(IN) :: nnds
      INTEGER, DIMENSION(:), INTENT(IN) :: nverts      
      INTEGER, DIMENSION(:), INTENT(IN) :: el_type
      REAL(rp), DIMENSION(:,:,:), INTENT(IN) :: elxy
      REAL(rp), DIMENSION(:,:,:), INTENT(IN) :: psic
      REAL(rp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(OUT) :: xyhb
      REAL(rp), DIMENSION(:), INTENT(IN), OPTIONAL :: depth
      INTEGER, DIMENSION(:,:), INTENT(IN), OPTIONAL :: ect
      REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT), OPTIONAL :: elhb
      
      INTEGER :: interp
      INTEGER :: el,pt,nd
      INTEGER :: et,npts,mnnds,nv
      INTEGER :: et_linear
      INTEGER :: nnd_coord,nnd_interp
      REAL(rp) :: xpt,ypt,hb      
      REAL(rp) :: hbvert(4)
      REAL(rp), DIMENSION(:), ALLOCATABLE :: x,y
      
      
      interp = 0
      IF ( PRESENT(depth) .AND. PRESENT(ect) .AND. PRESENT(elhb) ) THEN
        interp = 1
      ENDIF
      
      
      mnnds = MAXVAL(nnds)      
      ALLOCATE(x(mnnds),y(mnnds))
      
      ALLOCATE(xyhb(mnnds,ne,2))
      ALLOCATE(elhb(mnnds,ne))

      DO el = 1,ne
      
        et = el_type(el)        
        nv = nverts(et)
        IF (mod(et,2) == 1) THEN
          npts = nnds(5) 
          et_linear = 1
        ELSE IF (mod(et,2) == 0) THEN
          npts = nnds(6) 
          et_linear = 2
        ENDIF                   
      
        nnd_coord = nnds(et)
        nnd_interp = nnds(et_linear)    
      
        DO pt = 1,npts              

          CALL element_transformation(nnd_coord,elxy(:,el,1),elxy(:,el,2),psic(:,pt,et),xpt,ypt)
          
!           elhb(pt,el) = 10d0
!           elhb(pt,el) = 10d0 - 5d0*cos(2d0*pi/500d0*ypt)        
          xyhb(pt,el,1) = xpt
          xyhb(pt,el,2) = ypt 
        
        
          IF (interp) THEN
            DO nd = 1,nv        
              hbvert(nd) = depth(ect(nd,el))         
            ENDDO        
                                              
            CALL bathymetry_interp_eval(nnd_interp,hbvert,psic(:,pt,et_linear),hb)
            
            elhb(pt,el) = hb            
          ENDIF
        
        ENDDO
      
      ENDDO
      
      END SUBROUTINE bathy_coordinates      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
      SUBROUTINE bathymetry_interp_eval(n,hb,l,h,dldr,dlds,drdx,drdy,dsdx,dsdy,Sp,dhdx,dhdy)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: n
      REAL(rp), DIMENSION(:), INTENT(IN) :: hb
      REAL(rp), DIMENSION(:), INTENT(IN) :: l
      REAL(rp), INTENT(OUT) :: h      
      REAL(rp), DIMENSION(:), INTENT(IN), OPTIONAL :: dldr,dlds
      REAL(rp), INTENT(IN), OPTIONAL :: drdx,drdy,dsdx,dsdy
      REAL(rp), INTENT(IN), OPTIONAL :: Sp
      REAL(rp), INTENT(OUT), OPTIONAL :: dhdx,dhdy
      
      INTEGER :: nd
      INTEGER :: calc_deriv
      
      IF( PRESENT(dldr) .AND. PRESENT(dlds) .AND. &
          PRESENT(drdx) .AND. PRESENT(drdy) .AND. &
          PRESENT(dsdx) .AND. PRESENT(dsdy) .AND. & 
          PRESENT(Sp)   .AND.                     &
          PRESENT(dhdx) .AND. PRESENT(dhdy) ) THEN
          
        calc_deriv = 1
          
      ELSE
      
        calc_deriv = 0
          
      ENDIF

      h = 0d0      
      DO nd = 1,n                
        h =  h + l(nd)*hb(nd)                           
      ENDDO         
      
      IF (calc_deriv == 1) THEN
      
        dhdx = 0d0
        dhdy = 0d0
        DO nd = 1,n  
          dhdx = dhdx + (dldr(nd)*drdx + dlds(nd)*dsdx)*hb(nd)*Sp
          dhdy = dhdy + (dldr(nd)*drdy + dlds(nd)*dsdy)*hb(nd)              
        ENDDO      
      
      ENDIF 
      
      
      RETURN
      END SUBROUTINE bathymetry_interp_eval
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      

      END MODULE bathymetry_interp_mod