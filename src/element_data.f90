      SUBROUTINE element_data()

      USE globals, ONLY: pres,ect,xy,depth,ne,nn,ned, &
                         area,edlen,edlen_area,normal,ged2nn,ged2el,ged2led, &
                         dhbdx,dhbdy,dhbdx_init,dhbdy_init

      IMPLICIT NONE
      INTEGER :: el,ed,led
      INTEGER :: alloc_status
      REAL(pres) :: x1,x2,x3,y1,y2,y3
      REAL(pres) :: hb1,hb2,hb3

      !!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Calculate element area 
      !!!!!!!!!!!!!!!!!!!!!!!!!!

      ALLOCATE(area(ne),STAT = alloc_status)
      IF(alloc_status /= 0) THEN
        PRINT*, 'Allocation error: area'
      ENDIF

      DO el = 1,ne
        x1 = xy(1,ect(1,el))
        y1 = xy(2,ect(1,el))

        x2 = xy(1,ect(2,el))
        y2 = xy(2,ect(2,el))

        x3 = xy(1,ect(3,el))
        y3 = xy(2,ect(3,el))

        area(el) = .5d0*((x2*y3-x3*y2) + (x3*y1-x1*y3) + (x1*y2-x2*y1))
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
      ENDDO

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Calculate bathymetry derivatives 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!       DO el = 1,nn
!         depth(el) = 10d0
!       ENDDO
      

      ALLOCATE(dhbdx(ne),dhbdy(ne),STAT = alloc_status)
      IF(alloc_status /= 0) THEN
        PRINT*, 'Allocation error: dhbdx,dhbdy'
      ENDIF
      
      ALLOCATE(dhbdx_init(ne),dhbdy_init(ne),STAT = alloc_status)
      IF(alloc_status /= 0) THEN
        PRINT*, 'Allocation error: dhbdx_init,dhbdy_init'
      ENDIF      

      DO el = 1,ne
        x1 = xy(1,ect(1,el))
        y1 = xy(2,ect(1,el))

        x2 = xy(1,ect(2,el))
        y2 = xy(2,ect(2,el))

        x3 = xy(1,ect(3,el))
        y3 = xy(2,ect(3,el))
        
        hb1 = depth(ect(1,el))
        hb2 = depth(ect(2,el))
        hb3 = depth(ect(3,el))

        dhbdx_init(el) = ( -(.5d0*(y3-y1)+.5d0*(y1-y2))*hb1 + .5d0*(y3-y1)*hb2 + .5d0*(y1-y2)*hb3 )/area(el)
        dhbdy_init(el) = ( -(.5d0*(x1-x3)+.5d0*(x2-x1))*hb1 + .5d0*(x1-x3)*hb2 + .5d0*(x2-x1)*hb3 )/area(el)
 
        WRITE(64,*) dhbdx(el),dhbdy(el)

      ENDDO

      RETURN
      END SUBROUTINE element_data