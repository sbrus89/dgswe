      MODULE plot_mod
      
      USE globals, ONLY: rp
      
      IMPLICIT NONE
      
      INTEGER :: unit_count = 99
      INTEGER :: nlev        
      REAL(rp), DIMENSION(:,:), ALLOCATABLE :: colors        
      
      CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      SUBROUTINE zoom_box(ne,el_type,nplt,xyplt,xbox_min,xbox_max,ybox_min,ybox_max,xmin,xmax,ymin,ymax,el_in)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: ne
      INTEGER, DIMENSION(:), INTENT(IN) :: el_type
      INTEGER, DIMENSION(:), INTENT(IN) :: nplt
      REAL(rp), DIMENSION(:,:,:), INTENT(IN) :: xyplt
      REAL(rp), INTENT(IN) :: xbox_min
      REAL(rp), INTENT(IN) :: xbox_max
      REAL(rp), INTENT(IN) :: ybox_min
      REAL(rp), INTENT(IN) :: ybox_max
      REAL(rp), INTENT(OUT) :: xmin
      REAL(rp), INTENT(OUT) :: xmax
      REAL(rp), INTENT(OUT) :: ymin
      REAL(rp), INTENT(OUT) :: ymax
      INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: el_in
      
      INTEGER :: el,pt
      INTEGER :: et,npts
      INTEGER :: outside
      REAL(rp) :: xpt,ypt
      
      
      xmax = -1d10
      ymax = -1d10
      xmin = 1d10
      ymin = 1d10      
      
      ALLOCATE(el_in(ne))
      el_in = 1
      
      DO el = 1,ne      
        et = el_type(el)                          
        npts = nplt(et)
        outside = 0
        DO pt = 1,npts  
        
          xpt = xyplt(pt,el,1)
          ypt = xyplt(pt,el,2)
        
          IF (xpt > xmax) THEN
            xmax = xpt
          ENDIF
          
          IF (xpt < xmin) THEN
            xmin = xpt
          ENDIF

          IF (ypt > ymax) THEN
            ymax = ypt
          ENDIF
          
          IF (ypt < ymin) THEN
            ymin = ypt
          ENDIF          
          
          IF (xpt < xbox_min .or. xpt > xbox_max) THEN
            outside = 1
            xmin = xbox_min
            xmax = xbox_max
          ENDIF
          
          IF (ypt < ybox_min .or. ypt > ybox_max) THEN
            outside = 1
            ymin = ybox_min
            ymax = ybox_max
          ENDIF

        ENDDO
        IF (outside == 1) THEN
          el_in(el) = 0
        ENDIF
      ENDDO    
            
      
      RETURN
      END SUBROUTINE zoom_box
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      SUBROUTINE scale_coordinates(ne,nn,el_type,nverts,nplt,xmin,xmax,ymin,ymax,xyplt,xy)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: ne
      INTEGER, INTENT(IN) :: nn
      INTEGER, DIMENSION(:), INTENT(IN) :: el_type
      INTEGER, DIMENSION(:), INTENT(IN) :: nverts
      INTEGER, DIMENSION(:), INTENT(IN) :: nplt
      REAL(rp), INTENT(IN) :: xmin
      REAL(rp), INTENT(IN) :: xmax
      REAL(rp), INTENT(IN) :: ymin
      REAL(rp), INTENT(IN) :: ymax
      REAL(rp), DIMENSION(:,:,:), INTENT(INOUT) :: xyplt
      REAL(rp), DIMENSION(:,:), INTENT(INOUT) :: xy
      
      INTEGER :: el,nd
      INTEGER :: et,npts,nv
      REAL(rp) :: rmin,rmax
      REAL(rp) :: smin,smax
      REAL(rp) :: ax,bx  
      REAL(rp) :: ay,by
      
      
      
!!!!!  Letter Size  !!!!      
!       rmin = 0d0
!       rmax = 612d0
!       smin = 0d0
!       smax = 792d0
      
      rmin = 20d0
      rmax = 592d0
      smin = 10d0
      smax = 772d0
      
      ax = (rmin/(xmin-xmax)+rmax/(xmax-xmin))
      bx = -(rmin*xmax/(xmin-xmax)+rmax*xmin/(xmax-xmin))
      
      ay = ax     ! axis equal
      by = smin-ax*ymin            
      
      DO el = 1,ne
        et = el_type(el)
        npts = nplt(et)
        nv = nverts(et)
        DO nd = 1,npts
          xyplt(nd,el,1) = ax*xyplt(nd,el,1) + bx
          xyplt(nd,el,2) = ay*xyplt(nd,el,2) + by                       
        ENDDO
      ENDDO
      
      DO nd = 1,nn
        xy(1,nd) = ax*xy(1,nd) + bx
        xy(2,nd) = ay*xy(2,nd) + by
      ENDDO      
      
      RETURN
      END SUBROUTINE scale_coordinates

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
      
      
      SUBROUTINE write_psheader(file_name,unit_number)
      
      IMPLICIT NONE
      
      CHARACTER(*), INTENT(IN) :: file_name
      INTEGER, INTENT(OUT) :: unit_number
      
      unit_count = unit_count + 1
      unit_number = unit_count
      
      OPEN(UNIT=unit_number,FILE=file_name)      
      WRITE(unit_number,"(A)") "%!PS-Adobe-3.0"
      
      WRITE(unit_number,"(A)") "/draw-element {"
      WRITE(unit_number,"(A)") "moveto"
      WRITE(unit_number,"(A)") "lineto"
      WRITE(unit_number,"(A)") "lineto"
      WRITE(unit_number,"(A)") "closepath"
      WRITE(unit_number,"(A)") ".5 setlinewidth 2 setlinejoin"
      WRITE(unit_number,"(A)") "stroke"      
      WRITE(unit_number,"(A)") "} def" 
      
      WRITE(unit_number,"(A)") "/trifill {"
      WRITE(unit_number,"(A)") "/A exch def"
      WRITE(unit_number,"(A)") "/B exch def"
      WRITE(unit_number,"(A)") "/C exch def"     
      WRITE(unit_number,"(A)") "/ds ["
      WRITE(unit_number,"(A)") "  0 A aload pop" 
      WRITE(unit_number,"(A)") "  0 B aload pop"
      WRITE(unit_number,"(A)") "  0 C aload pop"
      WRITE(unit_number,"(A)") "] def"
      WRITE(unit_number,"(A)") "newpath"
      WRITE(unit_number,"(A)") "<<"
      WRITE(unit_number,"(A)") "  /ShadingType 4"
      WRITE(unit_number,"(A)") "  /ColorSpace [ /DeviceRGB ]"
      WRITE(unit_number,"(A)") "  /DataSource ds" 
!       WRITE(unit_number,"(A)")   " /AntiAlias true"
      WRITE(unit_number,"(A)") ">>"
      WRITE(unit_number,"(A)") "shfill"
      WRITE(unit_number,"(A)") "} def"        
      
      RETURN
      END SUBROUTINE write_psheader
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!              
      
      SUBROUTINE close_ps(unit_number)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: unit_number

      WRITE(unit_number,"(A)") "showpage"
      CLOSE(unit_number)            
      
      RETURN
      END SUBROUTINE close_ps           
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      SUBROUTINE evaluate_depth_solution(ne,el_type,el_in,nplt,ndof,phi,snap,H,H_val)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: ne
      INTEGER, DIMENSION(:), INTENT(IN) :: el_type
      INTEGER, DIMENSION(:), INTENT(IN) :: el_in
      INTEGER, DIMENSION(:), INTENT(IN) :: nplt
      INTEGER, DIMENSION(:), INTENT(IN) :: ndof
      REAL(rp), DIMENSION(:,:,:), INTENT(IN) :: phi
      INTEGER, INTENT(IN) :: snap
      REAL(rp), DIMENSION(:,:,:), INTENT(IN) :: H
      REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: H_val
      
      INTEGER :: el,nd,dof
      INTEGER :: et,npts,ndf
      INTEGER :: mnpp
      
      IF ( .NOT. ALLOCATED(H_val)) THEN
        mnpp = MAXVAL(nplt)
        ALLOCATE(H_val(mnpp,ne)) 
      ENDIF

 elem:DO el = 1,ne
 
        IF(el_in(el) == 0) THEN
          CYCLE elem
        ENDIF
        
        et = el_type(el)
        npts = nplt(et)
        ndf = ndof(et)
        DO nd = 1,npts
          H_val(nd,el) = 0d0
          DO dof = 1,ndf
            H_val(nd,el) = H_val(nd,el) + H(dof,el,snap)*phi(dof,nd,et)           
          ENDDO            
        ENDDO
        
      ENDDO elem      
      
      
      RETURN
      END SUBROUTINE evaluate_depth_solution      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE evaluate_velocity_solution(ne,el_type,el_in,nplt,ndof,phi,snap,Qx,Qy,Z_val,hb_val,vel_val)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: ne
      INTEGER, DIMENSION(:), INTENT(IN) :: el_type
      INTEGER, DIMENSION(:), INTENT(IN) :: el_in
      INTEGER, DIMENSION(:), INTENT(IN) :: nplt
      INTEGER, DIMENSION(:), INTENT(IN) :: ndof
      REAL(rp), DIMENSION(:,:,:), INTENT(IN) :: phi
      INTEGER, INTENT(IN) :: snap
      REAL(rp), DIMENSION(:,:,:), INTENT(IN) :: Qx
      REAL(rp), DIMENSION(:,:,:), INTENT(IN) :: Qy      
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: Z_val
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: hb_val  
      REAL(rp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: vel_val      
      
      REAL(rp) :: Qx_val,Qy_val,H_val
      INTEGER :: el,nd,dof
      INTEGER :: et,npts,ndf
      INTEGER :: mnpp
      
      IF ( .NOT. ALLOCATED(vel_val)) THEN
        mnpp = MAXVAL(nplt)
        ALLOCATE(vel_val(mnpp,ne)) 
      ENDIF

 elem:DO el = 1,ne
 
        IF(el_in(el) == 0) THEN
          CYCLE elem
        ENDIF
        
        et = el_type(el)
        npts = nplt(et)
        ndf = ndof(et)
        DO nd = 1,npts
          Qx_val = 0d0
          Qy_val = 0d0
          DO dof = 1,ndf
            Qx_val = Qx_val + Qx(dof,el,snap)*phi(dof,nd,et)
            Qy_val = Qy_val + Qy(dof,el,snap)*phi(dof,nd,et)            
          ENDDO            
          H_val = Z_val(nd,el) + hb_val(nd,el)
          vel_val(nd,el) = sqrt((Qx_val/H_val)**2 + (Qy_val/H_val)**2)           
        ENDDO
        
      ENDDO elem      
      
      
      RETURN
      END SUBROUTINE evaluate_velocity_solution      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE read_colormap(cmap_file)
      
      IMPLICIT NONE
      
      CHARACTER(*), INTENT(IN) :: cmap_file
      
      INTEGER :: j,lev
      REAL(rp) :: cmax      
      LOGICAL :: file_exists
      
      INQUIRE(file=TRIM(ADJUSTL(cmap_file)),exist=file_exists)
      IF (file_exists == .FALSE.) THEN
        PRINT*, "colormap file does not exist"
        STOP
      ENDIF
      
      OPEN(UNIT=101,FILE=TRIM(ADJUSTL(cmap_file)))
      READ(101,*) nlev
      ALLOCATE(colors(nlev,3))
      DO lev = 1,nlev
        READ(101,*) (colors(lev,j), j=1,3)
      ENDDO
      CLOSE(101)
      
      cmax = MAXVAL(colors)
      IF (cmax > 1d0+1d-12) THEN
        colors = colors/255d0
      ENDIF      
      
      RETURN
      END SUBROUTINE read_colormap

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE plot_contours(file_unit,nplt,ntri,rect,ne,el_type,el_in,xy,sol_val)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: file_unit
      INTEGER, DIMENSION(:), INTENT(IN) :: nplt
      INTEGER, DIMENSION(:), INTENT(IN) :: ntri
      INTEGER, DIMENSION(:,:,:), INTENT(IN) :: rect
      INTEGER, INTENT(IN) :: ne
      INTEGER, DIMENSION(:), INTENT(IN) :: el_type
      INTEGER, DIMENSION(:), INTENT(IN) :: el_in
      REAL(rp), DIMENSION(:,:,:), INTENT(IN) :: xy  
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: sol_val
      
      INTEGER :: i,j,v
      INTEGER :: el,nd,dof,lev,tri
      INTEGER :: et,nv    
      REAL(rp) :: sol_min,sol_max,sol_lev
      REAL(rp) :: dc
   
      REAL(rp) :: color_val(3)




      sol_min = 999d0
      sol_max = -999d0
      
      DO el = 1,ne
        et = el_type(el)
        nv = nplt(et)
        DO nd = 1,nv
          IF (sol_val(nd,el) < sol_min) THEN
            sol_min = sol_val(nd,el)
          ENDIF
        
          IF (sol_val(nd,el) > sol_max) THEN
            sol_max = sol_val(nd,el)
          ENDIF          
        ENDDO
        
      ENDDO
      
      PRINT*, sol_min,sol_max
      
  


      
      dc = (sol_max-sol_min)/real(nlev-1,rp)
            
      
 elem:DO el = 1,ne
        et = el_type(el)
        IF (el_in(el) == 0) THEN
          CYCLE elem
        ENDIF
        
        DO tri = 1,ntri(et)

        DO v = 1,3  
   
          nd = rect(v,tri,et)
          
          sol_lev = sol_min
          lev = nlev
  levels: DO i = 1,nlev-1
            IF ((sol_val(nd,el) >= sol_lev) .and. (sol_val(nd,el) < sol_lev+dc)) THEN
              lev = i
              CALL interp_colors(lev,sol_lev,dc,colors,sol_val(nd,el),color_val)
              EXIT levels
            ENDIF
            sol_lev = sol_lev + dc
          ENDDO levels
          
          IF (sol_val(nd,el) <= sol_min) THEN
           lev = 1
         ELSE IF (sol_val(nd,el) > sol_max) THEN
           lev = nlev
         ENDIF          

          WRITE(file_unit,"(A,5(F9.5,1x),A)") "[",xy(nd,el,1),xy(nd,el,2),color_val(1),color_val(2),color_val(3),"]"  
        ENDDO        

        WRITE(file_unit,"(A)") "trifill"        
        
        ENDDO
      ENDDO elem


      END SUBROUTINE plot_contours
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   

      SUBROUTINE plot_mesh(file_unit,ne,nverts,el_type,el_in,xy,ect)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: file_unit
      INTEGER, INTENT(IN) :: ne
      INTEGER, DIMENSION(:), INTENT(IN) :: nverts
      INTEGER, DIMENSION(:), INTENT(IN) :: el_type
      INTEGER, DIMENSION(:), INTENT(IN) :: el_in      
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: xy
      INTEGER, DIMENSION(:,:), INTENT(IN) :: ect

      INTEGER :: el,nd
      INTEGER :: et
      
 elem:DO el = 1,ne
        et = el_type(el)
        IF (el_in(el) == 0) THEN
          CYCLE elem
        ENDIF
        DO nd = 1,nverts(et)
          WRITE(file_unit,"(2(F9.5,1x))") xy(1,ect(nd,el)),xy(2,ect(nd,el))    
        ENDDO        

        WRITE(file_unit,"(A)") "draw-element"
      ENDDO elem

      END SUBROUTINE plot_mesh
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   

      SUBROUTINE interp_colors(lev,sol_lev,dc,colors,sol_val,color_val)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: lev
      REAL(rp), INTENT(IN) :: sol_lev
      REAL(rp), INTENT(IN) :: dc
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: colors
      REAL(rp), INTENT(IN) :: sol_val
      REAL(rp), DIMENSION(:), INTENT(OUT) :: color_val
      
      INTEGER :: i
      REAL(rp) :: l1,l2
      
      l1 = (sol_val-(sol_lev+dc))/(-dc)
      l2 = (sol_val-sol_lev)/dc
      
      DO i = 1,3
        color_val(i) = l1*colors(lev,i) + l2*colors(lev+1,i)
      ENDDO      
      
      
!       DO i = 1,3
!         color_val(i) = colors(lev,i)
!       ENDDO
!       
      RETURN
      END SUBROUTINE interp_colors

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      END MODULE plot_mod