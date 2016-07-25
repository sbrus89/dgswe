      MODULE plot_mod
      
      USE globals, ONLY: rp
      
      IMPLICIT NONE
      
      INTEGER :: unit_count = 99
      INTEGER :: tex_output_unit = 10
      INTEGER :: nlev        
      REAL(rp), DIMENSION(:,:), ALLOCATABLE :: colors     
      INTEGER :: nline_texfile      
      
      REAL(rp) :: lr_margin 
      REAL(rp) :: cscale_width
      REAL(rp) :: axes_width
      REAL(rp) :: rmin_axes,rmax_axes
      REAL(rp) :: smin_axes,smax_axes    
      REAL(rp) :: rmin_cbar,rmax_cbar
      REAL(rp) :: smin_cbar,smax_cbar   
      REAL(rp) :: rmin_tbar,rmax_tbar
      REAL(rp) :: smin_tbar,smax_tbar        
      REAL(rp) :: ax,bx  
      REAL(rp) :: ay,by      
      
      REAL(rp) :: dash      
      REAL(rp) :: xticklabel_pad,yticklabel_pad
      REAL(rp) :: xlabel_pad,ylabel_pad  
      REAL(rp) :: cticklabel_pad
      REAL(rp) :: clabel_pad     
      REAL(rp) :: dr_xlabel,ds_ylabel,ds_clabel
      
      REAL(rp) :: rmin_page = 0d0
      REAL(rp) :: rmax_page = 612d0
      REAL(rp) :: smin_page = 0d0
      REAL(rp) :: smax_page = 792d0      
      
      INTEGER :: nxtick
      INTEGER :: nytick
      INTEGER :: nctick
      
      INTEGER :: nxdec
      INTEGER :: nydec
      INTEGER :: ncdec
      INTEGER :: ntdec
      
      CHARACTER(100) :: main_font = "/Times-Roman"
      CHARACTER(100) :: math_font = "/Times-Italic"    
!       CHARACTER(100) :: main_font = "(/usr/share/fonts/type1/gsfonts/cmr10.pfb)"
      INTEGER :: fontsize     
      
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
            outside = outside + 1
            xmin = xbox_min
            xmax = xbox_max
          ENDIF
          
          IF (ypt < ybox_min .or. ypt > ybox_max) THEN
            outside = outside + 1
            ymin = ybox_min
            ymax = ybox_max
          ENDIF

        ENDDO
        IF (outside >= npts) THEN
          el_in(el) = 0
        ENDIF
      ENDDO    
      
!       CALL in_element(xy,el_type,elxy,el_found,rs)
            
      
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
               
      lr_margin =(612d0 - figure_width)/2d0 
      axes_width = figure_width/1.37d0
      
      cscale_width = axes_width*.05d0
      
      dash = 0.2d0*cscale_width
      xticklabel_pad = 2d0*dash
      yticklabel_pad = 2d0*dash
      xlabel_pad = 10d0*dash
      ylabel_pad = 15d0*dash        
      cticklabel_pad = 2d0*dash
      clabel_pad = 15d0*dash 
      

         
      
      
      rmin_axes = rmin_page + lr_margin + ylabel_pad
      rmax_axes = rmax_page - lr_margin - 2d0*cscale_width - clabel_pad
      smin_axes = smin_page + lr_margin
      smax_axes = smax_page      
      
    
            
      ax = (rmin_axes/(xmin-xmax)+rmax_axes/(xmax-xmin))
      bx = -(rmin_axes*xmax/(xmin-xmax)+rmax_axes*xmin/(xmax-xmin))
      
      ay = ax     ! axis equal
      by = smin_axes-ax*ymin            
      
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
      
      smax_axes = ay*ymax + by
      
      IF (smax_axes > smax_page) THEN
        PRINT*, "Error: Equal axis scaling puts y axis off page"
        STOP
      ENDIF
      
      rmin_tbar = rmax_axes + cscale_width
      rmax_tbar = rmin_tbar + clabel_pad
      smin_tbar = smax_axes - 2d0*cscale_width
      smax_tbar = smin_tbar + cscale_width   
      
      RETURN
      END SUBROUTINE scale_coordinates

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
      
      
      SUBROUTINE write_psheader(file_name,file_unit)
      
      IMPLICIT NONE
      
      CHARACTER(*), INTENT(IN) :: file_name
      INTEGER, INTENT(OUT) :: file_unit
      
      CHARACTER(1000) :: line
      INTEGER :: read_stat       

      
      
      unit_count = unit_count + 1
      file_unit = unit_count
      
      OPEN(UNIT=file_unit,FILE=file_name)      
      WRITE(file_unit,"(A)") "%!PS-Adobe-3.0"
      
      WRITE(file_unit,"(A)") "/draw-element {"
      WRITE(file_unit,"(A)") "newpath"        
      WRITE(file_unit,"(A)") "moveto"
      WRITE(file_unit,"(A)") "lineto"
      WRITE(file_unit,"(A)") "lineto"
      WRITE(file_unit,"(A)") "closepath"
      WRITE(file_unit,"(A)") ".5 setlinewidth 2 setlinejoin"
      WRITE(file_unit,"(A)") "stroke"      
      WRITE(file_unit,"(A)") "} def" 
      
      WRITE(file_unit,"(A)") "/draw-box {"
      WRITE(file_unit,"(A)") "newpath"        
      WRITE(file_unit,"(A)") "moveto"
      WRITE(file_unit,"(A)") "lineto"
      WRITE(file_unit,"(A)") "lineto"
      WRITE(file_unit,"(A)") "lineto"      
      WRITE(file_unit,"(A)") "closepath"
      WRITE(file_unit,"(A)") ".5 setlinewidth 2 setlinejoin"
      WRITE(file_unit,"(A)") "stroke"      
      WRITE(file_unit,"(A)") "} def"   
      
      WRITE(file_unit,"(A)") "/draw-line {"
      WRITE(file_unit,"(A)") "newpath"      
      WRITE(file_unit,"(A)") "moveto"
      WRITE(file_unit,"(A)") "lineto"
      WRITE(file_unit,"(A)") ".5 setlinewidth"
      WRITE(file_unit,"(A)") "stroke"      
      WRITE(file_unit,"(A)") "} def"       
      
      WRITE(file_unit,"(A)") "/trifill {"
      WRITE(file_unit,"(A)") "/A exch def"
      WRITE(file_unit,"(A)") "/B exch def"
      WRITE(file_unit,"(A)") "/C exch def"     
      WRITE(file_unit,"(A)") "/ds ["
      WRITE(file_unit,"(A)") "  0 A aload pop" 
      WRITE(file_unit,"(A)") "  0 B aload pop"
      WRITE(file_unit,"(A)") "  0 C aload pop"
      WRITE(file_unit,"(A)") "] def"
      WRITE(file_unit,"(A)") "newpath"
      WRITE(file_unit,"(A)") "<<"
      WRITE(file_unit,"(A)") "  /ShadingType 4"
      WRITE(file_unit,"(A)") "  /ColorSpace [ /DeviceRGB ]"
      WRITE(file_unit,"(A)") "  /DataSource ds" 
!       WRITE(file_unit,"(A)")   " /AntiAlias true"
      WRITE(file_unit,"(A)") ">>"
      WRITE(file_unit,"(A)") "closepath"      
      WRITE(file_unit,"(A)") "shfill"
      WRITE(file_unit,"(A)") "} def"  
      
      WRITE(file_unit,"(A)") "/recfill {"
      WRITE(file_unit,"(A)") "/A exch def"
      WRITE(file_unit,"(A)") "/B exch def"
      WRITE(file_unit,"(A)") "/C exch def"
      WRITE(file_unit,"(A)") "/D exch def"      
      WRITE(file_unit,"(A)") "<<"
      WRITE(file_unit,"(A)") " /ShadingType 2" 
      WRITE(file_unit,"(A)") " /ColorSpace /DeviceRGB"
      WRITE(file_unit,"(A)") " /BBox A "      
      WRITE(file_unit,"(A)") " /Coords B "
      WRITE(file_unit,"(A)") " /Function"
      WRITE(file_unit,"(A)") " <<"
      WRITE(file_unit,"(A)") "  /FunctionType 2"
      WRITE(file_unit,"(A)") "  /Domain [0 1]"
      WRITE(file_unit,"(A)") "  /C0 C "
      WRITE(file_unit,"(A)") "  /C1 D "
      WRITE(file_unit,"(A)") "  /N 1"
      WRITE(file_unit,"(A)") " >>"
      WRITE(file_unit,"(A)") ">>" 
      WRITE(file_unit,"(A)") "shfill"
      WRITE(file_unit,"(A)") "} def"  
      
      ! center vertically, left align horizontally
      WRITE(file_unit,"(A)") "/caxis-tick-labels {"
      WRITE(file_unit,"(A)") "moveto gsave dup false charpath pathbbox grestore exch pop sub 2 div 0 exch rmoveto pop show"
      WRITE(file_unit,"(A)") "} def"
      

      ! top align vertically, center horizontally
      WRITE(file_unit,"(A)") "/xaxis-labels {"
      WRITE(file_unit,"(A)") "moveto gsave dup false charpath pathbbox grestore"
      WRITE(file_unit,"(A)") "3 1 roll"      
      WRITE(file_unit,"(A)") "4 1 roll"    
      WRITE(file_unit,"(A)") "sub neg" 
      WRITE(file_unit,"(A)") "3 1 roll" 
      WRITE(file_unit,"(A)") "sub 2 div neg" 
      WRITE(file_unit,"(A)") "exch" 
      WRITE(file_unit,"(A)") "rmoveto" 
      WRITE(file_unit,"(A)") "show" 
      WRITE(file_unit,"(A)") "} def"         
      
      ! center vertically, right align horizontally
      WRITE(file_unit,"(A)") "/yaxis-tick-labels {"
      WRITE(file_unit,"(A)") "moveto gsave dup false charpath pathbbox grestore"
      WRITE(file_unit,"(A)") "3 1 roll"      
      WRITE(file_unit,"(A)") "4 1 roll"    
      WRITE(file_unit,"(A)") "sub 2 div neg" 
      WRITE(file_unit,"(A)") "3 1 roll" 
      WRITE(file_unit,"(A)") "sub neg" 
      WRITE(file_unit,"(A)") "exch" 
      WRITE(file_unit,"(A)") "rmoveto" 
      WRITE(file_unit,"(A)") "show" 
      WRITE(file_unit,"(A)") "} def"  
      
      ! rotate text, center vertically, right align horizontally
      WRITE(file_unit,"(A)") "/yaxis-labels {"   
      WRITE(file_unit,"(A)") "moveto gsave dup false charpath pathbbox grestore"
      WRITE(file_unit,"(A)") "3 1 roll"      
      WRITE(file_unit,"(A)") "4 1 roll"    
      WRITE(file_unit,"(A)") "sub neg" 
      WRITE(file_unit,"(A)") "3 1 roll" 
      WRITE(file_unit,"(A)") "sub 2 div neg" 
      WRITE(file_unit,"(A)") "rmoveto" 
      WRITE(file_unit,"(A)") "90 rotate"        
      WRITE(file_unit,"(A)") "show" 
      WRITE(file_unit,"(A)") "-90 rotate"      
      WRITE(file_unit,"(A)") "} def"  
      
      ! rotate text, center vertically, left align horizontally      
      WRITE(file_unit,"(A)") "/caxis-labels {"   
      WRITE(file_unit,"(A)") "moveto gsave dup false charpath pathbbox grestore"
      WRITE(file_unit,"(A)") "3 1 roll"      
      WRITE(file_unit,"(A)") "4 1 roll"    
      WRITE(file_unit,"(A)") "sub " 
      WRITE(file_unit,"(A)") "3 1 roll" 
      WRITE(file_unit,"(A)") "sub 2 div neg" 
      WRITE(file_unit,"(A)") "rmoveto"  
      WRITE(file_unit,"(A)") "90 rotate"        
      WRITE(file_unit,"(A)") "show" 
      WRITE(file_unit,"(A)") "-90 rotate"      
      WRITE(file_unit,"(A)") "} def"      
      
      WRITE(file_unit,"(A)") "/choosefont {"
      WRITE(file_unit,"(A)") "findfont"      
      WRITE(file_unit,"(I3,A)") fontsize," scalefont"  
      WRITE(file_unit,"(A)") "setfont"
      WRITE(file_unit,"(A)") "} def"   
      
      WRITE(file_unit,"(A)") TRIM(ADJUSTL(main_font))//" choosefont"      
!       WRITE(file_unit,"(A)") "(/usr/share/fonts/type1/gsfonts/n021003l.pfb)     choosefont"    
!       WRITE(file_unit,"(A)") "(/usr/share/fonts/type1/gsfonts/cmr10.pfb)     choosefont"        
      
      
      
      OPEN(UNIT=tex_output_unit, FILE="axes.ps")
      WRITE(file_unit,"(A)") "%!PS-Adobe-3.0"      
      nline_texfile = 0
head: DO
        READ(tex_output_unit,"(A)",IOSTAT=read_stat) line
        nline_texfile = nline_texfile + 1        
        IF (read_stat < 0 ) THEN
          PRINT*, "Error in LaTeX file"
          STOP
        ELSEIF (TRIM(ADJUSTL(line)) == "%%EndSetup") THEN
          WRITE(file_unit,"(A)") line
          EXIT head
        ELSEIF (nline_texfile > 1) THEN
          WRITE(file_unit,"(A)") line                        
        ENDIF
                
      ENDDO head       

     
      WRITE(file_unit,"(A)") "gsave newpath"        
      WRITE(file_unit,"(2(F9.5,1x),A)") rmin_axes,smin_axes,"moveto"
      WRITE(file_unit,"(2(F9.5,1x),A)") rmax_axes,smin_axes,"lineto"
      WRITE(file_unit,"(2(F9.5,1x),A)") rmax_axes,smax_axes,"lineto"      
      WRITE(file_unit,"(2(F9.5,1x),A)") rmin_axes,smax_axes,"lineto"      
      WRITE(file_unit,"(A)") "closepath"     
      WRITE(file_unit,"(A)") ".75 .75 .75 setrgbcolor fill grestore"      
  
      RETURN
      END SUBROUTINE write_psheader
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE write_texheader(file_name,unit_number)
      
      IMPLICIT NONE
      
      CHARACTER(*), INTENT(IN) :: file_name      
      INTEGER, INTENT(INOUT) :: unit_number 
      
      unit_count = unit_count + 1
      unit_number = unit_count      
      
      OPEN(UNIT=unit_number,FILE=file_name)           
      WRITE(unit_number,"(A,I2,A)") "\documentclass[",fontsize,"pt]{article}"
      WRITE(unit_number,"(A)") "\pagestyle{empty}"
!       WRITE(unit_number,"(A)") "\usepackage[absolute,showboxes]{textpos}"
      WRITE(unit_number,"(A)") "\usepackage[absolute]{textpos}"       
      WRITE(unit_number,"(A)") "\usepackage{graphicx}"  
      WRITE(unit_number,"(A)") "\usepackage{xcolor}"      
      WRITE(unit_number,"(A)") "\setlength{\TPHorizModule}{1pt}"  
      WRITE(unit_number,"(A)") "\setlength{\TPVertModule}{1pt}"
      WRITE(unit_number,"(A)") "\textblockorigin{0pt}{0pt}"
      WRITE(unit_number,"(A)") "\setlength{\parindent}{0pt}"      
      WRITE(unit_number,"(A)") "\begin{document}"
       
      
      RETURN
      END SUBROUTINE write_texheader
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!              
      
      SUBROUTINE close_ps(filename,unit_number)
      
      IMPLICIT NONE
      
      CHARACTER(*), INTENT(IN) :: filename
      INTEGER, INTENT(IN) :: unit_number
      
      CHARACTER(1000) :: line
      INTEGER :: read_stat   
      
tail: DO
        READ(10,"(A)",IOSTAT=read_stat) line
        nline_texfile = nline_texfile + 1        
        IF (read_stat < 0 ) THEN
          EXIT tail                
        ENDIF
                
        WRITE(unit_number,"(A)") line
                
      ENDDO tail           
   
 
      CLOSE(tex_output_unit)      

!       WRITE(unit_number,"(A)") "showpage"
      CLOSE(unit_number)            
!       CALL SYSTEM("ps2eps -B -C < "//filename//".ps > "//filename//".eps")
      
      RETURN
      END SUBROUTINE close_ps           
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE close_tex(unit_number)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(INOUT) :: unit_number      
      
      WRITE(unit_number,"(A)") "\end{document}"
      CLOSE(unit_number)
      
      RETURN
      END SUBROUTINE close_tex      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      SUBROUTINE evaluate_depth_solution(ne,el_type,el_in,nplt,ndof,phi,snap,H,H_val,H_min,H_max)
      
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
      REAL(rp), INTENT(OUT) :: H_min
      REAL(rp), INTENT(OUT) :: H_max
      
      INTEGER :: el,nd,dof
      INTEGER :: et,npts,ndf
      INTEGER :: mnpp
      
      IF ( .NOT. ALLOCATED(H_val)) THEN
        mnpp = MAXVAL(nplt)
        ALLOCATE(H_val(mnpp,ne)) 
      ENDIF
      
      H_min = 1d10
      H_max = -1d10

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
          
          IF (H_val(nd,el) < H_min) THEN
            H_min = H_val(nd,el)
          ENDIF
          
          IF (H_val(nd,el) > H_max) THEN
            H_max = H_val(nd,el)
          ENDIF
        ENDDO
        
      ENDDO elem      
      
      
      RETURN
      END SUBROUTINE evaluate_depth_solution      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE evaluate_velocity_solution(ne,el_type,el_in,nplt,ndof,phi,snap,Qx,Qy,Z_val,hb_val,vel_val,vel_min,vel_max)
      
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
      REAL(rp), INTENT(OUT) :: vel_min
      REAL(rp), INTENT(OUT) :: vel_max
      
      REAL(rp) :: Qx_val,Qy_val,H_val
      INTEGER :: el,nd,dof
      INTEGER :: et,npts,ndf
      INTEGER :: mnpp
      
      IF ( .NOT. ALLOCATED(vel_val)) THEN
        mnpp = MAXVAL(nplt)
        ALLOCATE(vel_val(mnpp,ne)) 
      ENDIF
      
      vel_min = 1d10
      vel_max = -1d10

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
          
          IF(vel_val(nd,el) < vel_min) THEN
            vel_min = vel_val(nd,el)
          ENDIF
          
          IF (vel_val(nd,el) > vel_max) THEN
            vel_max = vel_val(nd,el)
          ENDIF
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

      SUBROUTINE write_colorscale(file_unit)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: file_unit

      
      INTEGER :: lev
      INTEGER :: tick         
      REAL(rp) :: r0,r1,s0,s1
      REAL(rp) :: ds
      REAL(rp) :: cval
      REAL(rp) :: dc   
      CHARACTER(20) :: cchar
      

      
      r0 = rmax_axes
      s0 = smin_axes
      r1 = rmax_axes            

      
      ds = (smax_cbar-smin_cbar)/(nlev-1)
      DO lev = 1,nlev-1
      
        s1 = s0 + ds
      
        WRITE(file_unit,"(A,3(F9.5,1x),A)") "[",colors(lev+1,1),colors(lev+1,2),colors(lev+1,3),"]"       
        WRITE(file_unit,"(A,3(F9.5,1x),A)") "[",colors(lev,1),colors(lev,2),colors(lev,3),"]"  
        WRITE(file_unit,"(A,4(F9.5,1x),A)") "[",r0,s0,r1,s1,"]" 
        WRITE(file_unit,"(A,4(F9.5,1x),A)") "[",rmin_cbar,s0,rmax_cbar,s1,"]"         
        WRITE(file_unit,"(A)") "recfill" 
        
        s0 = s1
        
      ENDDO
                     
      WRITE(file_unit,"(2(F9.5,1x))") rmin_cbar,s1
      WRITE(file_unit,"(2(F9.5,1x))") rmax_cbar,s1
      WRITE(file_unit,"(2(F9.5,1x))") rmax_cbar,smin_axes
      WRITE(file_unit,"(2(F9.5,1x))") rmin_cbar,smin_axes      
      WRITE(file_unit,"(A)") "draw-box"       
      

      
      s0 = smin_axes 
      ds = (smax_cbar-smin_cbar)/(nctick-1)      
      DO tick = 1,nctick
              
        WRITE(file_unit,"(2(F9.5,1x))") rmin_cbar,s0 
        WRITE(file_unit,"(2(F9.5,1x))") rmin_cbar+dash,s0     
        WRITE(file_unit,"(A)") "draw-line" 
        
        WRITE(file_unit,"(2(F9.5,1x))") rmax_cbar,s0
        WRITE(file_unit,"(2(F9.5,1x))") rmax_cbar-dash,s0     
        WRITE(file_unit,"(A)") "draw-line"  

        s0 = s0 + ds
      ENDDO       
           
      
      RETURN
      END SUBROUTINE write_colorscale
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE write_caxis_labels(file_unit,time_bar,sol_min,sol_max,sol_label)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: file_unit
      INTEGER, INTENT(IN) :: time_bar
      REAL(rp), INTENT(IN) :: sol_min
      REAL(rp), INTENT(IN) :: sol_max
      CHARACTER(*), INTENT(IN) :: sol_label

      
      INTEGER :: lev
      INTEGER :: tick         
      REAL(rp) :: r0,r1,s0,s1
      REAL(rp) :: ds
      REAL(rp) :: dash
      REAL(rp) :: cval
      REAL(rp) :: dc   
      CHARACTER(20) :: cchar
      
      rmin_cbar = rmax_axes + cscale_width
      rmax_cbar = rmin_cbar + cscale_width
      smin_cbar = smin_axes
      IF (time_bar == 1) THEN      
!         smax_cbar = .75d0*smax_axes     
        smax_cbar = smin_tbar - 2d0*cscale_width            
      ELSE 
        smax_cbar = smax_axes 
      ENDIF     
            

      r0 = rmax_axes
      s0 = smin_axes
      r1 = rmax_axes            
      
      cval = sol_min
      dc = (sol_max-sol_min)/(nctick-1)
      ds = (smax_cbar-smin_cbar)/(nctick-1)      
      DO tick = 1,nctick               
        
        CALL format_number(ncdec,cval,0,cchar)         
!         WRITE(file_unit,"(A,2(1x,F9.5))")  "("//TRIM(ADJUSTL(cchar))//")",rmax_cbar+cticklabel_pad,s0       
!         WRITE(file_unit,"(A)") "caxis-tick-labels"
        
        WRITE(file_unit,"(A,F9.5,A,F9.5,A)") "\begin{textblock}{400}[0,0.5](",rmax_cbar+cticklabel_pad,",",smax_page-s0,")"
        WRITE(file_unit,"(A)") TRIM(ADJUSTL(cchar))        
        WRITE(file_unit,"(A)") "\end{textblock}"                 
        
        cval = cval + dc
        s0 = s0 + ds
      ENDDO       
      
!         WRITE(file_unit,"(A,2(1x,F9.5))")  "("//TRIM(ADJUSTL(sol_label))//")",rmax_cbar+clabel_pad,(smin_axes+smax_axes)/2d0       
!         WRITE(file_unit,"(A)") "caxis-labels"     
        
      WRITE(file_unit,"(A,F9.5,A,F9.5,A)") "\begin{textblock}{50}[0,0.5](",-50+rmax_cbar+clabel_pad,",",smax_page-(smax_cbar+smin_cbar)/2d0,")"
      WRITE(file_unit,"(A)") "\hfill \rotatebox[origin=c]{90}{"//TRIM(ADJUSTL(sol_label))//"}\vskip-\TPboxrulesize"       
      WRITE(file_unit,"(A)") "\end{textblock}"         
      
      RETURN
      END SUBROUTINE write_caxis_labels     

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE write_xyaxis_labels(file_name,file_unit)
      
      IMPLICIT NONE
            
      CHARACTER(*), INTENT(IN) :: file_name      
      INTEGER, INTENT(INOUT) :: file_unit

      INTEGER :: i
      INTEGER :: expnt
      REAL(rp) :: r0,r1,s0,s1
      REAL(rp) :: xval,yval
    
      CHARACTER(20) :: xchar,ychar     
      
      CALL write_texheader(file_name,file_unit)      
      

      dr_xlabel = (rmax_axes-rmin_axes)/(real(nxtick,rp)-1d0)
      

      r0 = rmin_axes
      
      xval = (r0-bx)/ax
      expnt = INT(LOG10(xval))
      IF (expnt <= 3) THEN
        expnt = 0
      ENDIF
      
      
      DO i = 1,nxtick             
        
        xval = (r0-bx)/ax
        
        CALL format_number(nxdec,xval,expnt,xchar)        
      
                  
        WRITE(file_unit,"(A,F9.5,A,F9.5,A)") "\begin{textblock}{400}[0.5,0](",r0,",",smax_page-smin_axes+xticklabel_pad,")"
        WRITE(file_unit,"(A)") "\centerline{"//TRIM(ADJUSTL(xchar))//"}"
        WRITE(file_unit,"(A)") "\end{textblock}"        
        
!         WRITE(file_unit,"(A,2(1x,F9.5))")  "("//TRIM(ADJUSTL(xchar))//")",r0,smin_axes-xticklabel_pad    
!         WRITE(file_unit,"(A)") "xaxis-labels"           
        
        r0 = r0 + dr_xlabel
      ENDDO
      
      
      WRITE(file_unit,"(A,F9.5,A,F9.5,A)") "\begin{textblock}{400}[0.5,0](",(rmax_axes+rmin_axes)/2d0,",",smax_page-smin_axes+xlabel_pad,")"
      WRITE(file_unit,"(A)") "\centerline{$x$}"        
      WRITE(file_unit,"(A)") "\end{textblock}"  
      
!       WRITE(file_unit,"(A)") TRIM(ADJUSTL(math_font))//" choosefont"         
!       WRITE(file_unit,"(A,2(1x,F9.5))")  "(x)",(rmax_axes+rmin_axes)/2d0,smin_axes-xlabel_pad    
!       WRITE(file_unit,"(A)") "xaxis-labels"      
!       WRITE(file_unit,"(A)") TRIM(ADJUSTL(main_font))//" choosefont"        
 
      
      IF (expnt /= 0) THEN      
        WRITE(file_unit,"(A,F9.5,A,F9.5,A)") "\begin{textblock}{400}[0,0](",rmin_axes+(nxtick-1)*dr_xlabel,",",smax_page-smin_axes+xticklabel_pad+xlabel_pad/2d0,")"
        WRITE(file_unit,"(A,I2,A)") "$\times10^{",expnt,"}$"        
        WRITE(file_unit,"(A)") "\end{textblock}"         
      ENDIF
      
      
      
      
      
      
      IF (nytick >= 10000) THEN
        ds_ylabel = dr_xlabel      
      ELSE
        ds_ylabel = (smax_axes-smin_axes)/(real(nytick,rp)-1d0)      
      ENDIF
      
      
      s0 = smin_axes      
      
      yval = (s0-by)/ay
      expnt = INT(LOG10(yval))
      IF (expnt <= 3) THEN
        expnt = 0
      ENDIF      
          
      
      DO i = 1,nytick
      
        IF (s0 > smax_axes) THEN
          EXIT        
        ENDIF
        
        yval = (s0-by)/ay
        
        CALL format_number(nydec,yval,expnt,ychar) 
        
        WRITE(file_unit,"(A,F9.5,A,F9.5,A)") "\begin{textblock}{400}[1,0.5](",rmin_axes-yticklabel_pad,",",smax_page-s0,")"
        WRITE(file_unit,"(A)") "\hfill "//TRIM(ADJUSTL(ychar))        
        WRITE(file_unit,"(A)") "\end{textblock}"  
        
!         WRITE(file_unit,"(A,2(1x,F9.5))")  "("//TRIM(ADJUSTL(ychar))//")",rmin_axes-yticklabel_pad,s0        
!         WRITE(file_unit,"(A)") "yaxis-tick-labels"       

        s0 = s0 + ds_ylabel        
      
      ENDDO         
      
      WRITE(file_unit,"(A,F9.5,A,F9.5,A)") "\begin{textblock}{50}[0,0.5](",rmin_axes-ylabel_pad,",",smax_page-(smax_axes+smin_axes)/2d0,")"
      WRITE(file_unit,"(A)") "\rotatebox[origin=c]{90}{$y$}\vskip-\TPboxrulesize"        
      WRITE(file_unit,"(A)") "\end{textblock}"  
      
      
!       WRITE(file_unit,"(A)") TRIM(ADJUSTL(math_font))//" choosefont"        
!       WRITE(file_unit,"(A,2(1x,F9.5))")  "(y)",rmin_axes-ylabel_pad,(smax_axes+smin_axes)/2d0    
!       WRITE(file_unit,"(A)") "yaxis-labels"  
!       WRITE(file_unit,"(A)") TRIM(ADJUSTL(main_font))//" choosefont"        
      
      IF (expnt /= 0) THEN      
        WRITE(file_unit,"(A,F9.5,A,F9.5,A)") "\begin{textblock}{400}[0,1](",rmin_axes,",",smax_page-(smax_axes + yticklabel_pad),")"
        WRITE(file_unit,"(A,I2,A)") "$\times10^{",expnt,"}$"        
        WRITE(file_unit,"(A)") "\end{textblock}"         
      ENDIF
      
      RETURN
      END SUBROUTINE write_xyaxis_labels

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE write_xyaxis(file_unit)
      
      IMPLICIT NONE
            
      INTEGER, INTENT(IN) :: file_unit
      
      INTEGER :: i
      INTEGER :: expnt
      REAL(rp) :: dr,ds
      REAL(rp) :: r0,r1,s0,s1
      REAL(rp) :: xval,yval
    
      CHARACTER(20) :: xchar,ychar
      
      WRITE(file_unit,"(A)")  "newpath"     
      WRITE(file_unit,"(2(F9.5,1x),A)") rmin_page,smin_page," moveto"
      WRITE(file_unit,"(2(F9.5,1x),A)") rmax_page,smin_page," lineto"
      WRITE(file_unit,"(2(F9.5,1x),A)") rmax_page,smax_page," lineto"  
      WRITE(file_unit,"(2(F9.5,1x),A)") rmin_page,smax_page," lineto"    
      WRITE(file_unit,"(2(F9.5,1x),A)") rmin_page,smin_page," lineto"        
      WRITE(file_unit,"(2(F9.5,1x),A)") rmin_axes,smin_axes," moveto"   
      WRITE(file_unit,"(2(F9.5,1x),A)") rmin_axes,smax_axes," lineto"  
      WRITE(file_unit,"(2(F9.5,1x),A)") rmax_axes,smax_axes," lineto"
      WRITE(file_unit,"(2(F9.5,1x),A)") rmax_axes,smin_axes," lineto"
      WRITE(file_unit,"(2(F9.5,1x),A)") rmin_axes,smin_axes," lineto"     
      WRITE(file_unit,"(A)") "gsave 1 1 1 setrgbcolor fill grestore"      
      WRITE(file_unit,"(A)") "clear"      
      
      
      WRITE(file_unit,"(2(F9.5,1x))") rmin_axes,smin_axes
      WRITE(file_unit,"(2(F9.5,1x))") rmax_axes,smin_axes
      WRITE(file_unit,"(2(F9.5,1x))") rmax_axes,smax_axes
      WRITE(file_unit,"(2(F9.5,1x))") rmin_axes,smax_axes      
      WRITE(file_unit,"(A)") "draw-box"  
      
!       ! x-axis line
!       WRITE(file_unit,"(2(F9.5,1x))") rmin_axes,smin_axes
!       WRITE(file_unit,"(2(F9.5,1x))") rmax_axes,smin_axes     
!       WRITE(file_unit,"(A)") "draw-line"  
!       
!       ! y-axis line
!       WRITE(file_unit,"(2(F9.5,1x))") rmin_axes,smin_axes
!       WRITE(file_unit,"(2(F9.5,1x))") rmin_axes,smax_axes     
!       WRITE(file_unit,"(A)") "draw-line"        

      
      dr = (rmax_axes-rmin_axes)/(real(nxtick,rp)-1d0)
      

      r0 = rmin_axes
      
      xval = (r0-bx)/ax
      expnt = INT(LOG10(xval))
      IF (expnt <= 3) THEN
        expnt = 0
      ENDIF
      
      DO i = 1,nxtick        
        WRITE(file_unit,"(2(F9.5,1x))") r0,smin_axes+dash
        WRITE(file_unit,"(2(F9.5,1x))") r0,smin_axes
        WRITE(file_unit,"(A)") "draw-line" 
        
        WRITE(file_unit,"(2(F9.5,1x))") r0,smax_axes
        WRITE(file_unit,"(2(F9.5,1x))") r0,smax_axes-dash
        WRITE(file_unit,"(A)") "draw-line"         

        r0 = r0 + dr
      ENDDO
      
      
      
      
      
      IF (nytick >= 10000) THEN
        ds = dr   
      ELSE
        ds = (smax_axes-smin_axes)/(real(nytick,rp)-1d0)      
      ENDIF      

      s0 = smin_axes      
      
      yval = (s0-by)/ay
      expnt = INT(LOG10(yval))
      IF (expnt <= 3) THEN
        expnt = 0
      ENDIF      
      
      DO i = 1,nytick
      
        IF (s0 > smax_axes) THEN
          EXIT        
        ENDIF
        
        WRITE(file_unit,"(2(F9.5,1x))") rmin_axes,s0
        WRITE(file_unit,"(2(F9.5,1x))") rmin_axes+dash,s0
        WRITE(file_unit,"(A)") "draw-line" 
        
        WRITE(file_unit,"(2(F9.5,1x))") rmax_axes,s0
        WRITE(file_unit,"(2(F9.5,1x))") rmax_axes-dash,s0
        WRITE(file_unit,"(A)") "draw-line"         
 
        s0 = s0 + ds        
      
      ENDDO 
      
      
      RETURN
      END SUBROUTINE write_xyaxis

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE write_tbar_labels(file_unit,t_snap)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: file_unit
      REAL(rp), INTENT(IN) :: t_snap
      
      REAL(rp) :: tday
      CHARACTER(20) :: tchar
      
      tday = t_snap/86400d0
      
      CALL format_number(ntdec,tday,0,tchar)       
      
      WRITE(file_unit,"(A,F9.5,A,F9.5,A)") "\begin{textblock}{400}[0.5,0](",(rmin_tbar+rmax_tbar)/2d0,",",smax_page-(smax_axes),")"
      WRITE(file_unit,"(A,F9.2,A)") "\centerline{$t="//TRIM(ADJUSTL(tchar))//"$ days}"        
      WRITE(file_unit,"(A)") "\end{textblock}"      
      
      RETURN
      END SUBROUTINE write_tbar_labels
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE write_tbar(file_unit,t_snap,t_start,t_end)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: file_unit      
      REAL(rp), INTENT(IN) :: t_snap
      REAL(rp), INTENT(IN) :: t_start
      REAL(rp), INTENT(IN) :: t_end
      
      REAL(rp) ::  rmax_snap
      
      WRITE(file_unit,"(2(F9.5,1x))") rmin_tbar,smin_tbar
      WRITE(file_unit,"(2(F9.5,1x))") rmax_tbar,smin_tbar
      WRITE(file_unit,"(2(F9.5,1x))") rmax_tbar,smax_tbar
      WRITE(file_unit,"(2(F9.5,1x))") rmin_tbar,smax_tbar      
      WRITE(file_unit,"(A)") "draw-box"  
      
      rmax_snap = (t_snap-t_end)/(t_start-t_end)*rmin_tbar + (t_snap-t_start)/(t_end-t_start)*rmax_tbar
      
      WRITE(file_unit,"(A)")  "newpath"     
      WRITE(file_unit,"(2(F9.5,1x),A)") rmin_tbar,smin_tbar," moveto"
      WRITE(file_unit,"(2(F9.5,1x),A)") rmax_snap,smin_tbar," lineto"
      WRITE(file_unit,"(2(F9.5,1x),A)") rmax_snap,smax_tbar," lineto"  
      WRITE(file_unit,"(2(F9.5,1x),A)") rmin_tbar,smax_tbar," lineto"    
      WRITE(file_unit,"(A)") "closepath"            
      WRITE(file_unit,"(A)") "gsave 0 0 0 setrgbcolor fill grestore"        
      
      RETURN
      END SUBROUTINE write_tbar      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      SUBROUTINE latex_all_labels(sol_min,sol_max,sol_label,t_snap,t_start,t_end)
      
      IMPLICIT NONE
      
      REAL(rp), INTENT(IN), OPTIONAL :: sol_min
      REAL(rp), INTENT(IN), OPTIONAL :: sol_max
      CHARACTER(*), INTENT(IN), OPTIONAL :: sol_label
      REAL(rp), INTENT(IN), OPTIONAL :: t_snap
      REAL(rp), INTENT(IN), OPTIONAL :: t_start
      REAL(rp), INTENT(IN), OPTIONAL :: t_end 
      
      INTEGER :: time_bar
      INTEGER :: color_bar
      INTEGER :: tex_unit  
      
      color_bar = 0
      IF (PRESENT(sol_min) .AND. PRESENT(sol_max) .AND. PRESENT(sol_label)) THEN
        color_bar = 1
      ENDIF          
      
      time_bar = 0
      IF (PRESENT(t_snap) .AND. PRESENT(t_start) .AND. PRESENT(t_end)) THEN
        time_bar = 1
      ENDIF            
      
      CALL write_xyaxis_labels("axes.tex",tex_unit)  
      IF (color_bar == 1) THEN      
        CALL write_caxis_labels(tex_unit,time_bar,sol_min,sol_max,sol_label)
      ENDIF
      IF (time_bar == 1) THEN
        CALL write_tbar_labels(tex_unit,t_snap)
      ENDIF      
      CALL close_tex(tex_unit)
      
!       CALL SYSTEM("latex axes.tex")
!       CALL SYSTEM("dvips -o axes.ps axes.dvi")
      
      CALL SYSTEM("latex axes.tex >/dev/null")
      CALL SYSTEM("dvips -q -o axes.ps axes.dvi")      
      
          
      
      RETURN
      END SUBROUTINE latex_all_labels
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      SUBROUTINE write_all_axes(file_unit,t_snap,t_start,t_end)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: file_unit      
      REAL(rp), OPTIONAL :: t_snap
      REAL(rp), OPTIONAL :: t_start
      REAL(rp), OPTIONAL :: t_end 
      
      INTEGER :: time_bar
      
      time_bar = 0
      IF (PRESENT(t_snap) .AND. PRESENT(t_start) .AND. PRESENT(t_end)) THEN
        time_bar = 1
      ENDIF            
      
      CALL write_xyaxis(file_unit)        
      CALL write_colorscale(file_unit)
      IF (time_bar == 1) THEN
        CALL write_tbar(file_unit,t_snap,t_start,t_end)
      ENDIF            
      
      RETURN
      END SUBROUTINE write_all_axes
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE plot_contours(file_unit,nplt,ntri,rect,ne,el_type,el_in,xy,sol_val,sol_min,sol_max)
      
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
      REAL(rp), INTENT(IN) :: sol_min
      REAL(rp), INTENT(IN) :: sol_max

      
      INTEGER :: i,j,v
      INTEGER :: el,nd,dof,lev,tri
      INTEGER :: et,nv    
      REAL(rp) :: sol_lev
      REAL(rp) :: dc

   
      REAL(rp) :: color_val(3)

      
      
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
           color_val(1) = colors(lev,1)
           color_val(2) = colors(lev,2)
           color_val(3) = colors(lev,3)
         ELSE IF (sol_val(nd,el) > sol_max) THEN
           lev = nlev
           color_val(1) = colors(lev,1)
           color_val(2) = colors(lev,2)
           color_val(3) = colors(lev,3)           
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
      
      SUBROUTINE format_number(ndec,val,expnt,val_char)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: ndec
      REAL(rp), INTENT(IN) :: val
      INTEGER, INTENT(IN) :: expnt
      CHARACTER(*), INTENT(OUT) :: val_char
      
      CHARACTER(:), ALLOCATABLE :: frmt   
      CHARACTER(1) :: ndec_char
      
      IF (ndec > 0) THEN             
        WRITE(ndec_char,"(I1)") ndec      
        frmt = "(F20."//ndec_char//")"      
      ELSE
        frmt = "(I20)"
      ENDIF      

      IF (ndec > 0) THEN
        WRITE(val_char,frmt) val/(10d0**expnt)      
      ELSE
        WRITE(val_char,frmt) NINT(val/(10d0**expnt))
      ENDIF      
      
      RETURN
      END SUBROUTINE format_number
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   

      SUBROUTINE convert_ps(filename,frmt,density,rm_ps)
      
      IMPLICIT NONE
      
      CHARACTER(*), INTENT(IN) :: filename
      CHARACTER(*), INTENT(IN) :: frmt
      CHARACTER(*), INTENT(IN) :: density
      INTEGER, INTENT(IN) :: rm_ps
      
      CHARACTER(:), ALLOCATABLE :: command

      
      IF (TRIM(ADJUSTL(frmt)) == "png" .or. TRIM(ADJUSTL(frmt)) == "jpg" .or. TRIM(ADJUSTL(frmt)) == "tif") THEN
        command = "convert -trim -density "//density//" "//filename//".ps "//filename//"."//frmt
        PRINT("(A)"), "  Converting to "//frmt//" format..."
        CALL SYSTEM(command)
        IF (rm_ps == 1) THEN
          PRINT("(A)"), "  Removing PostScript file..."
          CALL SYSTEM("rm "//filename//".ps")
        ENDIF        
      ENDIF
      
      
      RETURN
      END SUBROUTINE convert_ps
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      END MODULE plot_mod