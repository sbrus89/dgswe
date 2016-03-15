      MODULE curvilinear_nodes_mod

      
      USE globals, ONLY: rp,pi
      USE transformation, ONLY: element_transformation
      
      
      IMPLICIT NONE

      CONTAINS
      
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       



      SUBROUTINE bathy_coordinates_straight(ne,nnds,nverts,el_type,xy,ect,psiv,xyhb)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: ne   
      INTEGER, DIMENSION(:), INTENT(IN) :: nnds
      INTEGER, DIMENSION(:), INTENT(IN) :: nverts
      INTEGER, DIMENSION(:), INTENT(IN) :: el_type
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: xy
      INTEGER, DIMENSION(:,:), INTENT(IN) :: ect
      REAL(rp), DIMENSION(:,:,:), INTENT(IN) :: psiv
      REAL(rp), DIMENSION(:,:,:), INTENT(OUT) :: xyhb
      
      
      INTEGER :: el,pt,nd
      INTEGER :: et,nv,npts,mnnds
      REAL(rp) :: xpt,ypt      
      REAL(rp), DIMENSION(:), ALLOCATABLE :: x,y
      
      mnnds = MAXVAL(nnds)      
      ALLOCATE(x(mnnds),y(mnnds))
      

      DO el = 1,ne
      
        et = el_type(el)        
        nv = nverts(et)
        IF (mod(et,2) == 1) THEN
          et = 5
        ELSE IF (mod(et,2) == 0) THEN
          et = 6
        ENDIF             
      
        DO nd = 1,nv
          x(nd) = xy(1,ect(nd,el))
          y(nd) = xy(2,ect(nd,el))
        ENDDO        
      
        npts = nnds(et)         
      
        DO pt = 1,npts              

          CALL element_transformation(nv,x,y,psiv(:,pt,et),xpt,ypt)
          
!           elhb(pt,el) = 10d0
!           elhb(pt,el) = 10d0 - 5d0*cos(2d0*pi/500d0*ypt)        
          xyhb(pt,el,1) = xpt
          xyhb(pt,el,2) = ypt
        ENDDO   
      
      ENDDO
      
      END SUBROUTINE bathy_coordinates_straight
      
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  


      SUBROUTINE eval_coordinates_curved(ctp,nnds,nverts,el_type,xy,ect,fbseg,fbnds, &
                                         nnfbed,nfbedn,nfbednn,ged2el,ged2led, &
                                         psiv,psic,bndxy,elxy,xyhb)
     
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
      REAL(rp), DIMENSION(:,:,:), INTENT(IN) :: psic
      REAL(rp), DIMENSION(:,:,:,:), INTENT(IN) :: bndxy
      REAL(rp), DIMENSION(:,:,:), INTENT(INOUT) :: elxy
      REAL(rp), DIMENSION(:,:,:), INTENT(INOUT) :: xyhb
     
      INTEGER :: i,ed,nd
      INTEGER :: ged,seg,n1,el,led
      INTEGER :: n1ind
      REAL(rp), DIMENSION(2,ctp-1) :: segxy      
      
      
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
        
        CALL bathy_coordinates_curved(el,nnds,el_type,elxy,psic,xyhb)
                
        
      ENDDO     
     
     END SUBROUTINE eval_coordinates_curved


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE edge_coordinates_curved(el,ctp,led,nnds,nverts,el_type,xy,ect,segxy,psiv,elxy)
       
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
      IF (mod(et,2) == 1) THEN
        el_type(el) = 3 
        nnd = nnds(3)      
      ELSE IF (mod(et,2) == 0) THEN
        el_type(el) = 4
        nnd = nnds(4)                  
      ENDIF
        
      et = el_type(el)  
      nv = nverts(et)
        
      DO nd = 1,nv
        x(nd) = xy(1,ect(nd,el))
        y(nd) = xy(2,ect(nd,el))
      ENDDO        
        
      DO pt = 1,nnd               

        CALL element_transformation(nv,x,y,psiv(:,pt,et),xpt,ypt)
        
        elxy(pt,el,1) = xpt
        elxy(pt,el,2) = ypt
      ENDDO      
        
        
        
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

      SUBROUTINE bathy_coordinates_curved(el,nnds,el_type,elxy,psic,xyhb)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: el
      INTEGER, DIMENSION(:), INTENT(IN) :: nnds
      INTEGER, DIMENSION(:), INTENT(IN) :: el_type
      REAL(rp), DIMENSION(:,:,:), INTENT(IN) :: elxy
      REAL(rp), DIMENSION(:,:,:), INTENT(IN) :: psic
      REAL(rp), DIMENSION(:,:,:), INTENT(INOUT) :: xyhb

      INTEGER :: pt
      INTEGER :: et,nnd,npts
      REAL(rp) :: xpt,ypt
      
      et = el_type(el)
      nnd = nnds(et)        
      IF (mod(et,2) == 1) THEN
        npts = nnds(5)
        et = 5
      ELSE IF (mod(et,2) == 0) THEN
        npts = nnds(6)
        et = 6
      ENDIF
        
      DO pt = 1,npts              
      
        CALL element_transformation(nnd,elxy(:,el,1),elxy(:,el,2),psic(:,pt,et),xpt,ypt)
          
!         elhb(pt,el) = 10d0
!         elhb(pt,el) = 10d0 - 5d0*cos(2d0*pi/500d0*ypt)     
        xyhb(pt,el,1) = xpt
        xyhb(pt,el,2) = ypt          
      ENDDO        
      
      END SUBROUTINE bathy_coordinates_curved
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
      
      END MODULE curvilinear_nodes_mod