      MODULE plot_mod
      
      USE globals, ONLY: rp
      
      IMPLICIT NONE
      
      INTEGER :: unit_count = 99
      INTEGER :: nlev        
      REAL(rp), DIMENSION(:,:), ALLOCATABLE :: colors        
      
      REAL(rp) :: lr_margin 
      REAL(rp) :: cscale_width
      REAL(rp) :: rmin,rmax
      REAL(rp) :: smin,smax      
      REAL(rp) :: ax,bx  
      REAL(rp) :: ay,by      
      
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

      SUBROUTINE scale_coordinates(ne,nn,el_type,nverts,nplt,figure_width,xmin,xmax,ymin,ymax,xyplt,xy)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: ne
      INTEGER, INTENT(IN) :: nn
      INTEGER, DIMENSION(:), INTENT(IN) :: el_type
      INTEGER, DIMENSION(:), INTENT(IN) :: nverts
      INTEGER, DIMENSION(:), INTENT(IN) :: nplt
      REAL(rp), INTENT(IN) :: figure_width
      REAL(rp), INTENT(IN) :: xmin
      REAL(rp), INTENT(IN) :: xmax
      REAL(rp), INTENT(IN) :: ymin
      REAL(rp), INTENT(IN) :: ymax
      REAL(rp), DIMENSION(:,:,:), INTENT(INOUT) :: xyplt
      REAL(rp), DIMENSION(:,:), INTENT(INOUT) :: xy
      
      INTEGER :: el,nd
      INTEGER :: et,npts,nv


      
      
      
!!!!!  Letter Size  !!!!      
!       rmin = 0d0
!       rmax = 612d0
!       smin = 0d0
!       smax = 792d0
      
      cscale_width = figure_width*.05d0
      lr_margin =(612d0 - figure_width + 2d0*cscale_width)/2d0
      
      rmin = 0d0 + lr_margin
      rmax = 612d0 - lr_margin
      smin = 0d0 + lr_margin
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
      
      smax = ay*ymax + by
      
      RETURN
      END SUBROUTINE scale_coordinates

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
      
      
      SUBROUTINE write_psheader(file_name,unit_number)
      
      IMPLICIT NONE
      
      CHARACTER(*), INTENT(IN) :: file_name
      INTEGER, INTENT(OUT) :: unit_number
      
      CHARACTER(20) :: font
      INTEGER :: fontsize
      
      
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
      
      WRITE(unit_number,"(A)") "/draw-box {"
      WRITE(unit_number,"(A)") "moveto"
      WRITE(unit_number,"(A)") "lineto"
      WRITE(unit_number,"(A)") "lineto"
      WRITE(unit_number,"(A)") "lineto"      
      WRITE(unit_number,"(A)") "closepath"
      WRITE(unit_number,"(A)") ".5 setlinewidth 2 setlinejoin"
      WRITE(unit_number,"(A)") "stroke"      
      WRITE(unit_number,"(A)") "} def"   
      
      WRITE(unit_number,"(A)") "/draw-line {"
      WRITE(unit_number,"(A)") "moveto"
      WRITE(unit_number,"(A)") "lineto"
      WRITE(unit_number,"(A)") ".5 setlinewidth"
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
      
      WRITE(unit_number,"(A)") "/recfill {"
      WRITE(unit_number,"(A)") "/A exch def"
      WRITE(unit_number,"(A)") "/B exch def"
      WRITE(unit_number,"(A)") "/C exch def"
      WRITE(unit_number,"(A)") "/D exch def"      
      WRITE(unit_number,"(A)") "<<"
      WRITE(unit_number,"(A)") " /ShadingType 2" 
      WRITE(unit_number,"(A)") " /ColorSpace /DeviceRGB"
      WRITE(unit_number,"(A)") " /BBox A "      
      WRITE(unit_number,"(A)") " /Coords B "
      WRITE(unit_number,"(A)") " /Function"
      WRITE(unit_number,"(A)") " <<"
      WRITE(unit_number,"(A)") "  /FunctionType 2"
      WRITE(unit_number,"(A)") "  /Domain [0 1]"
      WRITE(unit_number,"(A)") "  /C0 C "
      WRITE(unit_number,"(A)") "  /C1 D "
      WRITE(unit_number,"(A)") "  /N 1"
      WRITE(unit_number,"(A)") " >>"
      WRITE(unit_number,"(A)") ">>" 
      WRITE(unit_number,"(A)") "shfill"
      WRITE(unit_number,"(A)") "} def"  
      
      WRITE(unit_number,"(A)") "/vertcenter-leftalign {"
      WRITE(unit_number,"(A)") "moveto gsave dup false charpath pathbbox grestore exch pop sub 2 div 0 exch rmoveto pop show"
      WRITE(unit_number,"(A)") "} def"
      
      WRITE(unit_number,"(A)") "/lowalign-horzcenter {"
      WRITE(unit_number,"(A)") "moveto gsave dup false charpath pathbbox grestore pop exch pop sub 2 div 0 rmoveto show"
      WRITE(unit_number,"(A)") "} def"      
      
      WRITE(unit_number,"(A)") "/upalign-horzcenter {"
      WRITE(unit_number,"(A)") "moveto gsave dup false charpath pathbbox grestore"
      WRITE(unit_number,"(A)") "3 1 roll"      
      WRITE(unit_number,"(A)") "4 1 roll"    
      WRITE(unit_number,"(A)") "sub neg" 
      WRITE(unit_number,"(A)") "3 1 roll" 
      WRITE(unit_number,"(A)") "sub 2 div neg" 
      WRITE(unit_number,"(A)") "exch" 
      WRITE(unit_number,"(A)") "rmoveto" 
      WRITE(unit_number,"(A)") "show" 
      WRITE(unit_number,"(A)") "} def"         
      
      WRITE(unit_number,"(A)") "/vertcenter-rightalign {"
      WRITE(unit_number,"(A)") "moveto gsave dup false charpath pathbbox grestore"
      WRITE(unit_number,"(A)") "3 1 roll"      
      WRITE(unit_number,"(A)") "4 1 roll"    
      WRITE(unit_number,"(A)") "sub 2 div neg" 
      WRITE(unit_number,"(A)") "3 1 roll" 
      WRITE(unit_number,"(A)") "sub neg" 
      WRITE(unit_number,"(A)") "exch" 
      WRITE(unit_number,"(A)") "rmoveto" 
      WRITE(unit_number,"(A)") "show" 
      WRITE(unit_number,"(A)") "} def"      
      
      font = "Times-Roman"
      fontsize = 12
      WRITE(unit_number,"(A)") "/"//TRIM(ADJUSTL(font))//" findfont"      
      WRITE(unit_number,"(I3,A)") fontsize," scalefont"  
      WRITE(unit_number,"(A)") "setfont"
      
      
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

      SUBROUTINE colorscale(nplt,ne,el_type,sol_val,sol_min,sol_max,dc)
      
      IMPLICIT NONE
      
      INTEGER, DIMENSION(:), INTENT(IN) :: nplt
      INTEGER, INTENT(IN) :: ne
      INTEGER, DIMENSION(:), INTENT(IN) :: el_type
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: sol_val
      REAL(rp), INTENT(OUT) :: sol_min
      REAL(rp), INTENT(OUT) :: sol_max
      REAL(rp), INTENT(OUT) :: dc
      
      INTEGER :: el,nd
      INTEGER :: et,nv
      
      
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
      
      
      RETURN
      END SUBROUTINE colorscale

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE write_colorscale(file_unit,sol_min,sol_max)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: file_unit
      REAL(rp), INTENT(IN) :: sol_min
      REAL(rp), INTENT(IN) :: sol_max

      
      INTEGER :: lev
      INTEGER :: tick      
      INTEGER :: nctick
      REAL(rp) :: r0,r1,s0,s1
      REAL(rp) :: rbox0,rbox1
      REAL(rp) :: ds
      REAL(rp) :: rdash
      REAL(rp) :: cval
      REAL(rp) :: dc      
      CHARACTER(20) :: cchar
      
      nctick = 10
      
      r0 = rmax
      s0 = smin
      r1 = rmax            

      rbox0 = r0 + cscale_width
      rbox1 = rbox0 + cscale_width
      
      ds = (smax-smin)/(nlev-1)
      DO lev = 1,nlev-1
      
        s1 = s0 + ds
      
        WRITE(file_unit,"(A,3(F9.5,1x),A)") "[",colors(lev+1,1),colors(lev+1,2),colors(lev+1,3),"]"       
        WRITE(file_unit,"(A,3(F9.5,1x),A)") "[",colors(lev,1),colors(lev,2),colors(lev,3),"]"  
        WRITE(file_unit,"(A,4(F9.5,1x),A)") "[",r0,s0,r1,s1,"]" 
        WRITE(file_unit,"(A,4(F9.5,1x),A)") "[",rbox0,s0,rbox1,s1,"]"         
        WRITE(file_unit,"(A)") "recfill" 
        
        s0 = s1
        
      ENDDO
                     
      WRITE(file_unit,"(2(F9.5,1x))") rbox0,s1
      WRITE(file_unit,"(2(F9.5,1x))") rbox1,s1
      WRITE(file_unit,"(2(F9.5,1x))") rbox1,smin
      WRITE(file_unit,"(2(F9.5,1x))") rbox0,smin      
      WRITE(file_unit,"(A)") "draw-box"       
      
      rdash = 0.2d0*cscale_width
      
      s0 = smin 
      cval = sol_min
      dc = (sol_max-sol_min)/(nctick-1)
      ds = (smax-smin)/(nctick-1)      
      DO tick = 1,nctick
              
        WRITE(file_unit,"(2(F9.5,1x))") rbox0,s0 
        WRITE(file_unit,"(2(F9.5,1x))") rbox0+rdash,s0     
        WRITE(file_unit,"(A)") "draw-line" 
        
        WRITE(file_unit,"(2(F9.5,1x))") rbox1,s0
        WRITE(file_unit,"(2(F9.5,1x))") rbox1-rdash,s0     
        WRITE(file_unit,"(A)") "draw-line"  
        
        WRITE(cchar,"(F20.2)") cval       
        WRITE(file_unit,"(A,2(1x,F9.5))")  "( "//TRIM(ADJUSTL(cchar))//")",rbox1,s0       
        WRITE(file_unit,"(A)") "vertcenter-leftalign"
        
        cval = cval + dc
        s0 = s0 + ds
      ENDDO       
      
      RETURN
      END SUBROUTINE write_colorscale

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE write_axis(file_unit)
      
      IMPLICIT NONE
            
      INTEGER, INTENT(IN) :: file_unit
      
      INTEGER :: nxtick
      INTEGER :: i
      INTEGER :: expnt
      REAL(rp) :: dr,ds
      REAL(rp):: dash
      REAL(rp) :: r0,r1,s0,s1
      REAL(rp) :: ramax,ramin
      REAL(rp) :: samax,samin
      REAL(rp) :: xval,yval
      REAL(rp) :: xticklabel_pad,yticklabel_pad
      REAL(rp) :: xlabel_pad,ylabel_pad      
      CHARACTER(20) :: xchar,ychar
      
      ! x-axis line
      WRITE(file_unit,"(2(F9.5,1x))") rmin,smin
      WRITE(file_unit,"(2(F9.5,1x))") rmax,smin     
      WRITE(file_unit,"(A)") "draw-line"  
      
      ! y-axis line
      WRITE(file_unit,"(2(F9.5,1x))") rmin,smin
      WRITE(file_unit,"(2(F9.5,1x))") rmin,smax     
      WRITE(file_unit,"(A)") "draw-line"        
      
      nxtick = 10

      dash = 0.2d0*cscale_width
      xticklabel_pad = 2d0*dash
      yticklabel_pad = 2d0*dash
      xlabel_pad = 10d0*dash
      ylabel_pad = 10d0*dash      
      
      ramax = rmax
      ramin = rmin
      
      dr = (ramax-ramin)/(real(nxtick,rp)-1d0)
      ds = dr
      
      s0 = smin
      s1 = smin + dash
      r0 = ramin
      
      xval = (r0-bx)/ax
      expnt = INT(LOG10(xval))
      IF (expnt <= 3) THEN
        expnt = 0
      ENDIF
      
      DO i = 1,nxtick        
        WRITE(file_unit,"(2(F9.5,1x))") r0,s1
        WRITE(file_unit,"(2(F9.5,1x))") r0,s0
        WRITE(file_unit,"(A)") "draw-line" 
        
        xval = (r0-bx)/ax
        WRITE(xchar,"(F20.2)") xval/(10d0**expnt)       
        WRITE(file_unit,"(A,2(1x,F9.5))")  "("//TRIM(ADJUSTL(xchar))//")",r0,s0-xticklabel_pad    
        WRITE(file_unit,"(A)") "upalign-horzcenter"
        r0 = r0 + dr
      ENDDO
      
      
      WRITE(file_unit,"(A,2(1x,F9.5))")  "(x)",(rmax+rmin)/2d0,s0-xlabel_pad    
      WRITE(file_unit,"(A)") "upalign-horzcenter"      
      
      
      r0 = rmin
      r1 = rmin + dash
      
      yval = (s0-by)/ay
      expnt = INT(LOG10(yval))
      IF (expnt <= 3) THEN
        expnt = 0
      ENDIF      
      
      DO 
      
        IF (s0 > smax) THEN
          EXIT        
        ENDIF
        
        WRITE(file_unit,"(2(F9.5,1x))") r0,s0
        WRITE(file_unit,"(2(F9.5,1x))") r1,s0
        WRITE(file_unit,"(A)") "draw-line" 
        
        yval = (s0-by)/ay
        WRITE(ychar,"(F20.2)") yval/(10d0**expnt)       
        WRITE(file_unit,"(A,2(1x,F9.5))")  "("//TRIM(ADJUSTL(ychar))//")",r0-yticklabel_pad,s0        
        WRITE(file_unit,"(A)") "vertcenter-rightalign"        
        s0 = s0 + dr        
      
      ENDDO 
      
!       WRITE(file_unit,"(A,2(1x,F9.5))")  "(x)",(rmax-rmin)/2d0,s0-xlabel_pad    
!       WRITE(file_unit,"(A)") "upalign-horzcenter"         
      
      RETURN
      END SUBROUTINE write_axis

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



      CALL colorscale(nplt,ne,el_type,sol_val,sol_min,sol_max,dc)
      CALL write_colorscale(file_unit,sol_min,sol_max)
      CALL write_axis(file_unit)
            
      
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