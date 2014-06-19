      SUBROUTINE element_data()

      USE globals, ONLY: pres,ect,xy,depth,ne,nn,ned,ndof,nqpta,qpta,nqpte, &
                         el_type,nel_type,p,ctp,nelnds, &
                         area,edlen,edlen_area,normal,ged2nn,ged2el,ged2led, &
                         dhbdx,dhbdy,dhbdx_init,dhbdy_init, &
                         dpdx_init,dpdy_init, &
                         dpdr,dpds,wpta, &
                         drdx,drdy,dsdx,dsdy, &
                         psia,dpsidr,dpsids,detJa,mmi,  &
                         psie,dpsidxi,detJe, &
                         nx_pt,ny_pt
                         
      USE basis, ONLY: 

      IMPLICIT NONE
      INTEGER :: el,ed,led,dof,pt,i,nd
      INTEGER :: ind,et
      INTEGER :: mnqpta,mnnds,mnqpte,mndof
      INTEGER :: np(4),nnds(4)
      INTEGER :: alloc_status
      REAL(pres) :: x1,x2,x3,x4,y1,y2,y3,y4
      REAL(pres) :: dxdr,dxds,dydr,dyds      
      REAL(pres) :: hb1,hb2,hb3,hb,dhbdx1,dhbdy1 
      REAL(pres), ALLOCATABLE, DIMENSION(:) :: r,s
      REAL(pres), ALLOCATABLE, DIMENSION(:,:) :: jac
      

      np(1) = p
      np(2) = p
      np(3) = ctp
      np(4) = ctp          

      nnds(1) = 3
      nnds(2) = 4
      nnds(3) = (ctp+1)*(ctp+2)/2
      nnds(4) = (ctp+1)*(ctp+1)
      
      mnqpta = maxval(nqpta)
      mnnds = maxval(nnds)
      mnqpte = maxval(nqpte)
      mndof = maxval(ndof)
      
      ALLOCATE(dhbdx(ne,mnqpta),dhbdy(ne,mnqpta),STAT = alloc_status)
      IF(alloc_status /= 0) THEN
        PRINT*, 'Allocation error: dhbdx,dhbdy'
      ENDIF
      
      ALLOCATE(dhbdx_init(ne,mnqpta),dhbdy_init(ne,mnqpta),STAT = alloc_status)
      IF(alloc_status /= 0) THEN
        PRINT*, 'Allocation error: dhbdx_init,dhbdy_init'
      ENDIF     
    
      
      ALLOCATE(psia(mnnds,mnqpta+4*mnqpte,nel_type),dpsidr(mnnds,mnqpta+4*mnqpte,nel_type),dpsids(mnnds,mnqpta+4*mnqpte,nel_type))
      ALLOCATE(detJa(ne,mnqpta),mmi(ne,mndof*mndof))
      ALLOCATE(nx_pt(ned,mnqpte),ny_pt(ned,mnqpte))
      nx_pt = 0d0
      ny_pt = 0d0    
      
      dhbdx_init = 0d0
      dhbdy_init = 0d0
      
      DO i = 1,nel_type
        CALL area_transformation(i,np(i),nnds(i),nqpta(i),nqpte(i))
      ENDDO
  
      
      mnnds = maxval(np)+1  
      
      ALLOCATE(psie(mnnds,mnqpte,nel_type),dpsidxi(mnnds,mnqpte,nel_type))
      ALLOCATE(detJe(ned,mnqpte))
      
      DO i = 1,nel_type
        CALL edge_transformation(i,np(i),nqpte(i))
      ENDDO
      

      !!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Calculate element area 
      !!!!!!!!!!!!!!!!!!!!!!!!!!

      ALLOCATE(area(ne),STAT = alloc_status)
      IF(alloc_status /= 0) THEN
        PRINT*, 'Allocation error: area'
      ENDIF
      
      ALLOCATE(r(mnqpta),s(mnqpta),jac(ne,mnqpta))
      DO i = 1,nqpta(2)
        r(i) = qpta(i,1,2)
        s(i) = qpta(i,2,2)
      ENDDO

      PRINT*, " "
      DO el = 1,ne
        x1 = xy(1,ect(1,el))
        y1 = xy(2,ect(1,el))

        x2 = xy(1,ect(2,el))
        y2 = xy(2,ect(2,el))

        x3 = xy(1,ect(3,el))
        y3 = xy(2,ect(3,el))

        area(el) = .5d0*((x2*y3-x3*y2) + (x3*y1-x1*y3) + (x1*y2-x2*y1))
        IF (el_type(el) == 1) THEN
!           PRINT("(I25,15(e25.16))"), el,area(el),(2d0*detJa(el,pt), pt = 1,nqpta(1))
        ELSE
          x4 = xy(1,ect(4,el))
          y4 = xy(2,ect(4,el))
          
          DO pt = 1,nqpta(2)
          
            dxdr = .25d0*((-1d0+s(pt))*x1 + ( 1d0-s(pt))*x2 + (1d0+s(pt))*x3 + (-1d0-s(pt))*x4)
            dxds = .25d0*((-1d0+r(pt))*x1 + (-1d0-r(pt))*x2 + (1d0+r(pt))*x3 + ( 1d0-r(pt))*x4)
            dydr = .25d0*((-1d0+s(pt))*y1 + ( 1d0-s(pt))*y2 + (1d0+s(pt))*y3 + (-1d0-s(pt))*y4)
            dyds = .25d0*((-1d0+r(pt))*y1 + (-1d0-r(pt))*y2 + (1d0+r(pt))*y3 + ( 1d0-r(pt))*y4)
            
            jac(el,pt) = dxdr*dyds - dxds*dydr
          ENDDO
         
!           PRINT("(I5,4(e23.14),10x,4(e23.14))"), el,(jac(el,pt), pt = 1,nqpta(2)),(detJa(el,pt), pt = 1,nqpta(2))         
        ENDIF
      ENDDO
     


      !!!!!!!!!!!!!!!!!!!!!!!!!
      ! Calculate edge length 
      !!!!!!!!!!!!!!!!!!!!!!!!!
      ALLOCATE(edlen(ned),edlen_area(2,ned),STAT = alloc_status)
!       ALLOCATE(edlen(ned),edlen_area(ne,3),STAT = alloc_status)
      IF(alloc_status /= 0) THEN
        PRINT*, 'Allocation error: edlen'
      ENDIF

      DO ed = 1,ned
        x1 = xy(1,ged2nn(1,ed))
        y1 = xy(2,ged2nn(1,ed))

        x2 = xy(1,ged2nn(2,ed))
        y2 = xy(2,ged2nn(2,ed))

        edlen(ed) = sqrt((x2-x1)**2 + (y2-y1)**2)
        
!         PRINT("(I5,e23.14,10x,4(e23.14))"), ed,0.5d0*edlen(ed),(detJe(ed,pt),pt = 1,nqpte(1))

        edlen_area(1,ed) = edlen(ed)/area(ged2el(1,ed))
        IF (ged2el(2,ed) /= 0) THEN
          edlen_area(2,ed) = edlen(ed)/area(ged2el(2,ed))   
        ENDIF    

! Set up for alternative imlementation where all edlen_area
! multiplications are done right before/during the edge integration instead of after numerical flux calculation
! The idea was to vectorize these multiplications.  Didn't show any improvement.
!         el = ged2el(1,ed)     
!         led = ged2led(1,ed)
!         edlen_area(el,led) = edlen(ed)/area(el)
! 
!         el = ged2el(2,ed)
!         led = ged2led(2,ed)
!         IF (el /= 0) THEN
!           edlen_area(el,led) = edlen(ed)/area(el)   
!         ENDIF 
      ENDDO
        

      !!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Calculate edge normals 
      !!!!!!!!!!!!!!!!!!!!!!!!!!

      ALLOCATE(normal(2,ned),STAT = alloc_status)
      IF(alloc_status /= 0) THEN
        PRINT*, 'Allocation error: edlen'
      ENDIF

      DO ed = 1,ned
        x1 = xy(1,ged2nn(1,ed))
        y1 = xy(2,ged2nn(1,ed))

        x2 = xy(1,ged2nn(2,ed))
        y2 = xy(2,ged2nn(2,ed))

        normal(1,ed) =  (y2-y1)/edlen(ed)
        normal(2,ed) = -(x2-x1)/edlen(ed)
        
!         PRINT("(I5,2(e23.14),10x,4(e23.14))"), ed,normal(1,ed),normal(2,ed),(nx_pt(ed,i),i=1,nqpte(1)), (ny_pt(ed,i),i=1,nqpte(1))
      ENDDO
            

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Calculate bathymetry derivatives 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!       DO el = 1,nn
!         depth(el) = 10d0
!       ENDDO
      

 

!       DO el = 1,ne
!         x1 = xy(1,ect(1,el))
!         y1 = xy(2,ect(1,el))
! 
!         x2 = xy(1,ect(2,el))
!         y2 = xy(2,ect(2,el))
! 
!         x3 = xy(1,ect(3,el))
!         y3 = xy(2,ect(3,el))
!         
!         hb1 = depth(ect(1,el))
!         hb2 = depth(ect(2,el))
!         hb3 = depth(ect(3,el))
! 
!         dhbdx_init(el) = ( -(.5d0*(y3-y1)+.5d0*(y1-y2))*hb1 + .5d0*(y3-y1)*hb2 + .5d0*(y1-y2)*hb3 )/area(el)
!         dhbdy_init(el) = ( -(.5d0*(x1-x3)+.5d0*(x2-x1))*hb1 + .5d0*(x1-x3)*hb2 + .5d0*(x2-x1)*hb3 )/area(el)
!         
! !         dhbdx(el) = ( -(.5d0*(y3-y1)+.5d0*(y1-y2))*hb1 + .5d0*(y3-y1)*hb2 + .5d0*(y1-y2)*hb3 )/area(el)
! !         dhbdy(el) = ( -(.5d0*(x1-x3)+.5d0*(x2-x1))*hb1 + .5d0*(x1-x3)*hb2 + .5d0*(x2-x1)*hb3 )/area(el)
!  
! !         WRITE(64,*) dhbdx(el),dhbdy(el)
! 
!       ENDDO

      DO el = 1,ne
        et = el_type(el)
        
        IF (et == 1) THEN
        
          x1 = xy(1,ect(1,el))
          y1 = xy(2,ect(1,el))

          x2 = xy(1,ect(2,el))
          y2 = xy(2,ect(2,el))

          x3 = xy(1,ect(3,el))
          y3 = xy(2,ect(3,el))
        
          hb1 = depth(ect(1,el))
          hb2 = depth(ect(2,el))
          hb3 = depth(ect(3,el))

          dhbdx1= ( -(.5d0*(y3-y1)+.5d0*(y1-y2))*hb1 + .5d0*(y3-y1)*hb2 + .5d0*(y1-y2)*hb3 )/area(el)
          dhbdy1 = ( -(.5d0*(x1-x3)+.5d0*(x2-x1))*hb1 + .5d0*(x1-x3)*hb2 + .5d0*(x2-x1)*hb3 )/area(el)
          
          
!           PRINT("(I5,e23.14,10x,3(e23.14))"), el,dhbdx1,(dhbdx_init(el,i),i=1,nqpta(1))
!           PRINT("(I5,e23.14,10x,3(e23.14))"), el,dhbdy1,(dhbdy_init(el,i),i=1,nqpta(1))          
        
        ENDIF
      ENDDO
      

      RETURN
      END SUBROUTINE element_data