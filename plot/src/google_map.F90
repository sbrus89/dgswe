      MODULE google_map

      USE globals, ONLY:rp,r_earth,pi,deg2rad
      USE ppm

      CONTAINS
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
      SUBROUTINE get_map(xmin,xmax,ymin,ymax,slam0,sphi0,lamc,phic,img,img_height,img_width,img_res)

      IMPLICIT NONE
      
      REAL(rp), INTENT(IN) :: xmin
      REAL(rp), INTENT(IN) :: xmax
      REAL(rp), INTENT(IN) :: ymin
      REAL(rp), INTENT(IN) :: ymax
      REAL(rp), INTENT(IN) :: slam0
      REAL(rp), INTENT(IN) :: sphi0  
      REAL(rp), INTENT(OUT) :: lamc
      REAL(rp), INTENT(OUT) :: phic
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE, INTENT(OUT) :: img
      INTEGER, INTENT(OUT) :: img_width
      INTEGER, INTENT(OUT) :: img_height
      REAL(rp), INTENT(OUT) :: img_res
      
      INTEGER :: i,j
      INTEGER :: avg
      INTEGER, PARAMETER :: tilesize = 256
      INTEGER, PARAMETER :: width = 640
      INTEGER, PARAMETER :: height = 640
      INTEGER, PARAMETER :: scale = 2
      INTEGER :: zoom_level
      REAL(rp) :: xres,yres
      REAL(rp) :: xc,yc
      REAL(rp) :: min_res,init_res
      CHARACTER(8192) :: url      
      CHARACTER(9000) :: command
      CHARACTER(10) :: lamc_char,phic_char
      CHARACTER(2) :: zoom_char
      CHARACTER(3) :: width_char,height_char
      CHARACTER(1) :: scale_char
      
      ! Determine zoom level and image resolution
      xres = (xmax-xmin)/real(width,rp)
      yres = (ymax-ymin)/real(height,rp)
      min_res = max(xres,yres)
      init_res = 2d0*pi*r_earth/real(tilesize,rp)
      zoom_level = floor(log2(init_res/min_res))
      img_res = init_res/2d0**zoom_level/real(scale,rp)      
      
      ! Find center coordinate
      xc = 0.5d0*(xmax+xmin)
      yc = 0.5d0*(ymax+ymin)      
      lamc = xc/(r_earth*cos(sphi0)) + slam0
      phic = yc/r_earth      
      lamc = lamc/deg2rad
      phic = phic/deg2rad
      
      ! Convert API parameters to characters
      WRITE(lamc_char,"(F10.6)") lamc
      WRITE(phic_char,"(F10.6)") phic
      WRITE(zoom_char,"(I2)") zoom_level
      WRITE(width_char,"(I3)") width
      WRITE(height_char,"(I3)") height
      WRITE(scale_char,"(I1)") scale
            
      ! Build API URL
      url = 'https://maps.googleapis.com/maps/api/staticmap?'       
      url = TRIM(url) // 'size=' // width_char // 'x' // height_char 
      url = TRIM(url) // '&scale=' // scale_char
      url = TRIM(url) // '&maptype=satellite'
!       url = TRIM(url) // '&maptype=terrain'      
      url = TRIM(url) // '&center=' // TRIM(ADJUSTL(phic_char)) // "," // TRIM(ADJUSTL(lamc_char))
      url = TRIM(url) // '&zoom=' // TRIM(ADJUSTL(zoom_char))
      
      ! Download the Google Maps image
      command = 'wget "' // TRIM(url) // '" -O map.jpg'
      PRINT("(A)"), TRIM(command)      
      CALL SYSTEM(TRIM(command))
      
      ! Convert to ppm and load into img array
      CALL SYSTEM("convert map.jpg map.ppm")      
      CALL loadppm("map.ppm",img,img_width,img_height)

      CALL lighten_image(1.25d0,img_width,img_height,img)      
!       CALL convert2grayscale(img_width,img_height,img)      
      CALL saveppm("out.ppm",img)

      END SUBROUTINE get_map
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE write_map(file_unit,lamc_deg,phic_deg,xmin,xmax,ymin,ymax,ax,bx,ay,by,slam0,sphi0, &
                           map,map_height,map_width,map_res)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: file_unit
      REAL(rp), INTENT(INOUT) :: lamc_deg      
      REAL(rp), INTENT(INOUT) :: phic_deg
      REAL(rp), INTENT(IN) :: xmin
      REAL(rp), INTENT(IN) :: xmax
      REAL(rp), INTENT(IN) :: ymin
      REAL(rp), INTENT(IN) :: ymax
      REAL(rp), INTENT(IN) :: ax
      REAL(rp), INTENT(IN) :: bx
      REAL(rp), INTENT(IN) :: ay
      REAL(rp), INTENT(IN) :: by
      REAL(rp), INTENT(IN) :: slam0
      REAL(rp), INTENT(IN) :: sphi0
      INTEGER, DIMENSION(:,:,:), INTENT(IN) :: map
      INTEGER, INTENT(IN) :: map_width
      INTEGER, INTENT(IN) :: map_height
      REAL(rp), INTENT(IN) :: map_res    
      
      INTEGER :: i,j,k
      REAL(rp) :: xc,yc   
      INTEGER :: xpixc,ypixc
      INTEGER :: inside
      REAL(rp) :: lamc,phic
      REAL(rp) :: xo,yo
      REAL(rp) :: ii,jj
      REAL(rp) :: x(4),y(4)
      REAL(rp) :: lon(4),lat(4)
      
      ! Determine center x,y coordinate
      lamc = lamc_deg*deg2rad
      phic = phic_deg*deg2rad      
      xc = r_earth*lamc
      yc = r_earth*log(tan(pi/4d0+phic/2d0))  
      
      ! Determine lower left x,y coordinate
      xpixc = map_width/2
      ypixc = map_height/2      
      xo = xc - real(xpixc,rp)*map_res
      yo = yc - real(ypixc,rp)*map_res      
      
      ! Write pixels to PS file
      DO j = 1,map_height
        DO i = 1,map_width    
          ii = real(i,rp)
          jj = real(j,rp)
          
          ! Calculate the corner coordinates of the pixel
          x(1) = (ii-1d0)
          y(1) = (jj-1d0) 
          x(2) = ii 
          y(2) = (jj-1d0) 
          x(3) = ii        
          y(3) = jj
          x(4) = (ii-1d0)
          y(4) = jj          
          DO k = 1,4
            x(k) = x(k)*map_res + xo
            y(k) = y(k)*map_res + yo
          ENDDO
          
          ! Convert to lon,lat using Mercator projection
          DO k = 1,4
            lon(k) = x(k)/r_earth
            lat(k) = 2d0*(atan(exp(y(k)/r_earth))-pi/4d0)
          ENDDO                 
          
          ! Convert to x,y using CPP projection
          DO k = 1,4         
            x(k) = r_earth*(lon(k)-slam0)*cos(sphi0)
            y(k) = lat(k)*r_earth
          ENDDO      
          
          ! Decide if pixel is inside plot box
          inside = 0
          IF (y(4) > ymin .and. y(1) < ymax .and. &
              x(2) > xmin .and. x(1) < xmax ) THEN              
           inside = 1
          ENDIF
          
          ! Write the PS code to fill the pixel
          IF(inside == 1) THEN
            WRITE(file_unit,"(3(F9.5,1x))") (real(map(k,i,j),rp)/255d0, k = 1,3)
            DO k = 1,4
              WRITE(file_unit,"(2(F9.5,1x))") ax*x(k)+bx,ay*y(k)+by
            ENDDO
            WRITE(file_unit,"(A)") "fill-pixel"
          ENDIF
          
        ENDDO
      ENDDO
      
      END SUBROUTINE write_map
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE convert2grayscale(img_width,img_height,img)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: img_width
      INTEGER, INTENT(IN) :: img_height
      INTEGER, DIMENSION(:,:,:), INTENT(INOUT) :: img
      
      INTEGER :: i,j
      REAL(rp) :: avg
      
      DO j = 1,img_height
        DO i = 1,img_width
         avg = INT(REAL(img(1,i,j)+img(2,i,j)+img(3,i,j),rp)/3d0)
        
         img(1,i,j) = avg
         img(2,i,j) = avg
         img(3,i,j) = avg
        ENDDO
      ENDDO      
      
      END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

      SUBROUTINE lighten_image(factor,img_width,img_height,img)
      
      IMPLICIT NONE
      
      REAL(rp), INTENT(IN) :: factor      
      INTEGER, INTENT(IN) :: img_width
      INTEGER, INTENT(IN) :: img_height
      INTEGER, DIMENSION(:,:,:), INTENT(INOUT) :: img     
      
      INTEGER :: i,j,k
      REAL(rp) :: new_color
      
      DO j = 1,img_height
        DO i = 1,img_width
         DO k = 1,3
           
           new_color = INT(real(img(k,i,j),rp)*factor)
           IF (new_color <= 255) THEN
             img(k,i,j) = new_color
           ELSE
             img(k,i,j) = 255
           ENDIF
         ENDDO
        ENDDO
      ENDDO            
      
      RETURN
      END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     REAL(rp) FUNCTION log2(x)
     
     IMPLICIT NONE
     
     REAL(rp), INTENT(IN) :: x
     
     log2 = LOG(x)/LOG(2d0)
     
     END FUNCTION log2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      END MODULE google_map