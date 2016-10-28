      MODULE plot_mod
      
      USE globals, ONLY: rp
      USE plot_globals, ONLY: plot_type,cscale_width,lr_margin,dash,fontsize, &
                              axes_width, &
                              rmin_page,rmax_page,smin_page,smax_page, &
                              rmin_axes,rmax_axes,smin_axes,smax_axes, &
                              rmin_cbar,rmax_cbar,smin_cbar,smax_cbar, &
                              rmin_tbar,rmax_tbar,smin_tbar,smax_tbar, &
                              xticklabel_pad,yticklabel_pad,cticklabel_pad, &
                              xlabel_pad,ylabel_pad,clabel_pad, &
                              nxtick,nytick,nctick, &
                              nxdec,nydec,ncdec,ntdec, &
                              dr_xlabel,ds_ylabel,ds_clabel,main_font, &
                              ncolors,colors, &
                              mnpp
      
      IMPLICIT NONE
      
      INTEGER :: unit_count = 99   
      INTEGER :: nline_texfile      
      INTEGER :: tex_output_unit = 11             
      
      CONTAINS
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE make_plot(snap,t_snap,fig,H1,H2,Qx,Qy)
            
      USE globals, ONLY: ne,nn,nverts,el_type,xy,ect, &
                         nbou,fbseg,fbnds, &
                         nnfbed,nfbedn, &
                         nfbed,fbedn, &
                         nobed,obedn, &
                         ged2el,ged2led      
      USE plot_globals, ONLY: el_in,t_start,t_end,xyplt,npplt,nptri,rect,r,s, &
                              frmt,density,pc
      USE evaluate_mod, ONLY: evaluate_depth_solution,evaluate_velocity_solution  
      USE labels_mod, ONLY: latex_axes_labels,run_latex, & 
                            latex_element_labels,latex_node_labels
      USE axes_mod, ONLY: write_all_axes
            
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: snap
      REAL(rp), INTENT(IN) :: t_snap
      TYPE(plot_type), INTENT(INOUT) :: fig      
      REAL(rp), DIMENSION(:,:), INTENT(INOUT), OPTIONAL :: H1           
      REAL(rp), DIMENSION(:,:), INTENT(INOUT), OPTIONAL :: H2      
      REAL(rp), DIMENSION(:,:), INTENT(INOUT), OPTIONAL :: Qx     
      REAL(rp), DIMENSION(:,:), INTENT(INOUT), OPTIONAL :: Qy
      
      CHARACTER(:), ALLOCATABLE :: filename
      CHARACTER(4) :: snap_char
      
  
      
      
      IF ((fig%type_flag == 2 .and. fig%plot_sol_option == 1) .or. fig%type_flag == 3) THEN
        PRINT("(A)"), "  Evaluating depth solution at additional plotting points..."
        CALL evaluate_depth_solution(ne,el_type,el_in,npplt,H1,fig)       
      ENDIF            
      
      IF (fig%plot_sol_option /= 1) THEN
        RETURN
      ENDIF     
      
      IF (fig%type_flag == 4) THEN
        PRINT("(A)"), "  Evaluating velocity solution at additional plotting points..."
        CALL evaluate_velocity_solution(ne,el_type,el_in,npplt,Qx,Qy,H1,H2,fig)
      ENDIF             
      
      
      IF (fig%cscale_option == "auto-snap") THEN      
        fig%sol_min = fig%snap_min
        fig%sol_max = fig%snap_max                 
      ELSE IF (fig%cscale_option == "auto-all") THEN      
        ! use existing values        
      ELSE IF (fig%cscale_option == "file") THEN      
!         DO i = 1,num_cscale_zeta_vals
!           IF (ABS(cscale_zeta_vals(i,1)-t_snap) < 1d-8) THEN
!             Z_min = cscale_zeta_vals(i,2)
!             Z_max = cscale_zeta_vals(i,3)              
!             EXIT
!          ENDIF
!        ENDDO            

        fig%sol_min = fig%cscale_vals(snap-1,2)
        fig%sol_max = fig%cscale_vals(snap-1,3)                                      
      ELSE     
        fig%sol_min = fig%cscale_min
        fig%sol_max = fig%cscale_max        
      ENDIF
      
      
          
      PRINT("(2(A,F10.5))"), "  min value = ", fig%sol_min, "  max value = ", fig%sol_max
      WRITE(fig%cscale_unit,"(3(e24.17,1x))") t_snap,fig%sol_min,fig%sol_max
       
      IF (snap > 0) THEN 
        WRITE(snap_char,"(I4.4)") snap-1           
        filename = TRIM(ADJUSTL(fig%name))//"_"//snap_char
      ELSE
        filename = TRIM(ADJUSTL(fig%name))
      ENDIF
      
      
          
      CALL latex_axes_labels(fig,t_snap,t_start,t_end)  
      IF (fig%el_label_option /= "off") THEN
        CALL latex_element_labels(ne,fig%el_label_option,el_type,el_in,nverts,xy,ect,nnfbed,nfbedn,ged2el)
      ENDIF
      IF (fig%nd_label_option /= "off") THEN
        CALL latex_node_labels(nn,fig%nd_label_option,xy,nbou,fbseg,fbnds)
      ENDIF
      CALL run_latex()    
      
      CALL write_psheader(filename//".ps",fig%ps_unit)          
      IF (fig%cbar_flag == 1) THEN
        CALL plot_filled_contours(fig%ps_unit,ne,el_type,el_in,nptri,rect,xyplt,fig)             
      ENDIF
      IF (fig%plot_lines_option == 1) THEN
        CALL plot_line_contours(fig%ps_unit,ne,el_type,el_in,nptri,rect,xyplt,r,s,snap,fig)          
      ENDIF
      IF (fig%plot_mesh_option == 1) THEN
        CALL fill_elements(fig%ps_unit,ne,nverts,pc,el_type,el_in,xy,ect,xyplt)      
!         CALL plot_mesh(fig%ps_unit,ne,nverts,pc,el_type,el_in,xy,ect,xyplt)
        CALL plot_boundaries(fig%ps_unit,nverts,pc,nobed,obedn,ged2el,ged2led,el_type,el_in,ect,xy,xyplt)       
        CALL plot_boundaries(fig%ps_unit,nverts,pc,nnfbed,nfbedn,ged2el,ged2led,el_type,el_in,ect,xy,xyplt)        
        CALL plot_boundaries(fig%ps_unit,nverts,pc,nfbed,fbedn,ged2el,ged2led,el_type,el_in,ect,xy,xyplt)          
      ENDIF 
      
      CALL write_all_axes(fig%ps_unit,fig%cbar_flag,fig%tbar_flag,t_snap,t_start,t_end)               
      CALL close_ps(filename,fig%ps_unit)
      CALL convert_ps(filename,frmt,density,fig%rm_ps)      
      
      RETURN
      END SUBROUTINE make_plot


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      SUBROUTINE zoom_box(ne,el_type,npplt,xyplt,xbox_min,xbox_max,ybox_min,ybox_max,xmin,xmax,ymin,ymax,el_in)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: ne
      INTEGER, DIMENSION(:), INTENT(IN) :: el_type
      INTEGER, DIMENSION(:), INTENT(IN) :: npplt
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
      INTEGER :: mnpts
      INTEGER, DIMENSION(:), ALLOCATABLE :: outside
      INTEGER :: in_flag
      REAL(rp) :: xpt,ypt     
      
      
      xmax = -1d10
      ymax = -1d10
      xmin = 1d10
      ymin = 1d10      
      
      mnpts = maxval(npplt)
      ALLOCATE(outside(mnpts))
      
      ALLOCATE(el_in(ne))
      el_in = 1
      
      DO el = 1,ne      
        et = el_type(el)                          
        npts = npplt(et)
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
            outside(pt) = 1
            xmin = xbox_min
            xmax = xbox_max
          ENDIF
          
          IF (ypt < ybox_min .or. ypt > ybox_max) THEN
            outside(pt) = 1
            ymin = ybox_min
            ymax = ybox_max
          ENDIF

        ENDDO
        
        in_flag = 0
        DO pt = 1,npts
          IF (outside(pt) == 0) THEN
            in_flag = 1
          ENDIF
        ENDDO
        
        IF (in_flag == 0) THEN
          el_in(el) = 0
        ENDIF
      ENDDO    
      
!       CALL in_element(xy,el_type,elxy,el_found,rs)
            
      
      RETURN
      END SUBROUTINE zoom_box
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      SUBROUTINE scale_factors(ne,nn,el_type,nverts,nnds,npplt,figure_width,xmin,xmax,ymin,ymax,ax,bx,ay,by)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: ne
      INTEGER, INTENT(IN) :: nn
      INTEGER, DIMENSION(:), INTENT(IN) :: el_type
      INTEGER, DIMENSION(:), INTENT(IN) :: nverts
      INTEGER, DIMENSION(:), INTENT(IN) :: nnds
      INTEGER, DIMENSION(:), INTENT(IN) :: npplt
      REAL(rp), INTENT(IN) :: figure_width
      REAL(rp), INTENT(IN) :: xmin
      REAL(rp), INTENT(IN) :: xmax
      REAL(rp), INTENT(IN) :: ymin
      REAL(rp), INTENT(IN) :: ymax
      REAL(rp), INTENT(OUT) :: ax
      REAL(rp), INTENT(OUT) :: bx
      REAL(rp), INTENT(OUT) :: ay
      REAL(rp), INTENT(OUT) :: by
!       REAL(rp), DIMENSION(:,:,:), INTENT(INOUT) :: xyplt
!       REAL(rp), DIMENSION(:,:), INTENT(INOUT) :: xy
!       REAL(rp), DIMENSION(:,:,:), INTENT(INOUT) :: elxy      
      
      INTEGER :: el,nd
      INTEGER :: et,npts,nv,nnd
               
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
      
!       DO el = 1,ne
!         et = el_type(el)
!         npts = npplt(et)
!         nv = nverts(et)
!         nnd = nnds(et)
!         
!         DO nd = 1,npts
!           xyplt(nd,el,1) = ax*xyplt(nd,el,1) + bx
!           xyplt(nd,el,2) = ay*xyplt(nd,el,2) + by                       
!         ENDDO
!         
!         DO nd = 1,nnd
!           elxy(nd,el,1) = ax*elxy(nd,el,1) + bx
!           elxy(nd,el,2) = ay*elxy(nd,el,2) + by
!         ENDDO
!       ENDDO
!       
!       DO nd = 1,nn
!         xy(1,nd) = ax*xy(1,nd) + bx
!         xy(2,nd) = ay*xy(2,nd) + by
!       ENDDO      
      
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
!       WRITE(file_unit,"(A)") "%!PS-Adobe-3.0 EPSF-3.0"      
!       WRITE(file_unit,"(A,4(F9.5,1x))") "%%BoundingBox: ",rmin_axes,smin_axes,rmax_axes,smax_axes                
      
      WRITE(file_unit,"(A)") "/draw-line {"
      WRITE(file_unit,"(A)") "newpath"        
      WRITE(file_unit,"(A)") "moveto"
      WRITE(file_unit,"(A)") "lineto"
      WRITE(file_unit,"(A)") ".5 setlinewidth 2 setlinejoin"
      WRITE(file_unit,"(A)") "stroke"      
      WRITE(file_unit,"(A)") "} def"       
      
      WRITE(file_unit,"(A)") "/draw-tri-element {"
      WRITE(file_unit,"(A)") "newpath"        
      WRITE(file_unit,"(A)") "moveto"
      WRITE(file_unit,"(A)") "lineto"
      WRITE(file_unit,"(A)") "lineto"
      WRITE(file_unit,"(A)") "closepath"
      WRITE(file_unit,"(A)") ".5 setlinewidth 2 setlinejoin"
      WRITE(file_unit,"(A)") "stroke"      
      WRITE(file_unit,"(A)") "} def" 
      
      WRITE(file_unit,"(A)") "/draw-quad-element {"
      WRITE(file_unit,"(A)") "newpath"        
      WRITE(file_unit,"(A)") "moveto"
      WRITE(file_unit,"(A)") "lineto"
      WRITE(file_unit,"(A)") "lineto"
      WRITE(file_unit,"(A)") "lineto"      
      WRITE(file_unit,"(A)") "closepath"
      WRITE(file_unit,"(A)") ".5 setlinewidth 2 setlinejoin"
      WRITE(file_unit,"(A)") "stroke"      
      WRITE(file_unit,"(A)") "} def"    
      
      WRITE(file_unit,"(A)") "/fill-tri-element {"
      WRITE(file_unit,"(A)") "newpath"        
      WRITE(file_unit,"(A)") "moveto"
      WRITE(file_unit,"(A)") "lineto"
      WRITE(file_unit,"(A)") "lineto"
      WRITE(file_unit,"(A)") "closepath"
      WRITE(file_unit,"(A)") "gsave 0 0 1 setrgbcolor fill grestore"     
      WRITE(file_unit,"(A)") "} def" 
      
      WRITE(file_unit,"(A)") "/fill-quad-element {"
      WRITE(file_unit,"(A)") "newpath"        
      WRITE(file_unit,"(A)") "moveto"
      WRITE(file_unit,"(A)") "lineto"
      WRITE(file_unit,"(A)") "lineto"
      WRITE(file_unit,"(A)") "lineto"      
      WRITE(file_unit,"(A)") "closepath"
      WRITE(file_unit,"(A)") "gsave 0 0 1 setrgbcolor fill grestore"      
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
      WRITE(file_unit,"(A)") ".5 setlinewidth 2 setlinejoin"
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
      
      
      
      OPEN(UNIT=tex_output_unit, FILE="labels.ps")
      WRITE(file_unit,"(A)") "%!PS-Adobe-3.0"      
!       WRITE(file_unit,"(A)") "%!PS-Adobe-3.0 EPSF-3.0"
!       WRITE(file_unit,"(A,4(F9.5,1x))") "%%BoundingBox: ",rmin_axes,smin_axes,rmax_axes,smax_axes           
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
      
      SUBROUTINE close_ps(filename,unit_number)
      
      IMPLICIT NONE
      
      CHARACTER(*), INTENT(IN) :: filename
      INTEGER, INTENT(IN) :: unit_number
      
      CHARACTER(1000) :: line
      INTEGER :: read_stat   
      
tail: DO
        READ(tex_output_unit,"(A)",IOSTAT=read_stat) line
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

      SUBROUTINE plot_filled_contours(file_unit,ne,el_type,el_in,nptri,rect,xy,fig)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: file_unit
      INTEGER, INTENT(IN) :: ne
      INTEGER, DIMENSION(:), INTENT(IN) :: el_type
      INTEGER, DIMENSION(:), INTENT(IN) :: el_in
      INTEGER, DIMENSION(:), INTENT(IN) :: nptri
      INTEGER, DIMENSION(:,:,:), INTENT(IN) :: rect 
      REAL(rp), DIMENSION(:,:,:), INTENT(IN) :: xy  
      TYPE(plot_type), INTENT(IN) :: fig

      
      INTEGER :: i,j,v
      INTEGER :: el,nd,dof,lev,tri
      INTEGER :: et,nv    
      REAL(rp) :: sol_lev
      REAL(rp) :: dc

   
      REAL(rp) :: color_val(3)

      
      
      dc = (fig%sol_max-fig%sol_min)/real(ncolors-1,rp)      
            
      
 elem:DO el = 1,ne
        et = el_type(el)
        IF (el_in(el) == 0) THEN
          CYCLE elem
        ENDIF
        
        DO tri = 1,nptri(et)

        DO v = 1,3  
   
          nd = rect(v,tri,et)
          
          sol_lev = fig%sol_min
          lev = ncolors
  levels: DO i = 1,ncolors-1
            IF ((fig%sol_val(nd,el) >= sol_lev) .and. (fig%sol_val(nd,el) < sol_lev+dc)) THEN
              lev = i
              CALL interp_colors(lev,sol_lev,dc,colors,fig%sol_val(nd,el),color_val)
              EXIT levels
            ENDIF
            sol_lev = sol_lev + dc
          ENDDO levels
          
         IF (fig%sol_val(nd,el) <= fig%sol_min) THEN
           lev = 1
           color_val(1) = colors(lev,1)
           color_val(2) = colors(lev,2)
           color_val(3) = colors(lev,3)
         ELSE IF (fig%sol_val(nd,el) > fig%sol_max) THEN
           lev = ncolors
           color_val(1) = colors(lev,1)
           color_val(2) = colors(lev,2)
           color_val(3) = colors(lev,3)           
         ENDIF          

          WRITE(file_unit,"(A,5(F9.5,1x),A)") "[",xy(nd,el,1),xy(nd,el,2),color_val(1),color_val(2),color_val(3),"]"  
        ENDDO        

        WRITE(file_unit,"(A)") "trifill"        
        
        ENDDO
      ENDDO elem


      END SUBROUTINE plot_filled_contours
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   

      SUBROUTINE plot_line_contours(file_unit,ne,el_type,el_in,nptri,rect,xyplt,rre,sre,snap,fig)
      
      USE globals, ONLY: mndof,elxy,xy,ect,np,mnnds
      USE plot_globals, ONLY: Z,hb,Qx,Qy
      USE read_dginp, ONLY: p,hbp
      USE basis, ONLY: element_basis,linear_basis
      USE transformation, ONLY: element_transformation
      USE shape_functions_mod, ONLY: shape_functions_area_eval, shape_functions_edge_eval
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: file_unit
      INTEGER, INTENT(IN) :: ne
      INTEGER, DIMENSION(:), INTENT(IN) :: el_type
      INTEGER, DIMENSION(:), INTENT(IN) :: el_in
      INTEGER, DIMENSION(:), INTENT(IN) :: nptri
      INTEGER, DIMENSION(:,:,:), INTENT(IN) :: rect
      REAL(rp), DIMENSION(:,:,:), INTENT(IN) :: xyplt  
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: rre
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: sre
      INTEGER, INTENT(IN) :: snap
      TYPE(plot_type), INTENT(IN) :: fig

      
      INTEGER :: i,j,vrt,it
      INTEGER :: el,nd,dof,lev,tri,ig
      INTEGER :: et,nv,ndf,ntry,nnd    
      REAL(rp) :: sol_lev
      REAL(rp) :: dc
      INTEGER :: above(3)
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: el_below
      INTEGER :: n1,n2
      INTEGER :: nd1,nd2
      REAL(rp) :: x(2),y(2)
      REAL(rp) :: r,s,re(2),se(2)
      REAL(rp) :: r1,r2,s1,s2
      REAL(rp) :: t,f,dfdr,dfds,dfdt,dt,t0
      INTEGER :: max_it
      INTEGER :: fail_flag
      REAL(rp) :: phi(mndof,2),dpdr(mndof,1),dpds(mndof,1)
      REAL(rp) :: l(mnnds,1)
      REAL(rp) :: xv(3),yv(3)
      REAL(rp) :: rv(3),sv(3)
      REAL(rp) :: zeta,bathy,xmom,ymom,H
      REAL(rp) :: dzdr,dzds,dbdr,dbds,dHdr,dHds
      REAL(rp) :: dxmdr,dxmds,dymdr,dymds
      REAL(rp) :: uvel,vvel
      REAL(rp) :: tol
      REAL(rp) :: hbe(2),Ze(2),ve(2),fe(2)
      REAL(rp) :: e,a,b,c,w,u,v
      
      tol = 1d-8
      ntry = 20
      max_it = 1000
      


      ALLOCATE(el_below(nptri(4),ne))
      el_below = 0
      
       dc = (fig%sol_max-fig%sol_min)/real(nctick-1,rp)                  
       sol_lev = fig%sol_min      
       
levels:DO lev = 1,nctick
      
    elem:DO el = 1,ne
           et = el_type(el)
           IF (el_in(el) == 0) THEN
             CYCLE elem
           ENDIF
     
        
    subtri:DO tri = 1,nptri(et)
    
             IF (el_below(tri,el) == 1) THEN
               CYCLE subtri
             ENDIF 

!              PRINT("(3(A,I9))"), "lev = ", lev, " el = ", el, " tri = ",tri
             
             DO vrt = 1,3  
   
               nd = rect(vrt,tri,et)
          
               IF (fig%sol_val(nd,el) >= sol_lev) THEN
                 above(vrt) = 1
               ELSE
                 above(vrt) = 0
               ENDIF
                
             ENDDO        
 
 
             IF (above(1) == 1 .and. above(2) == 1 .and. above(3) == 1) THEN
               CYCLE subtri
             ENDIF
             
             IF (above(1) == 0 .and. above(2) == 0 .and. above(3) == 0) THEN
               el_below(tri,el) = 1             
               CYCLE subtri
             ENDIF       
             
             
             DO i = 1,3
               nd = rect(i,tri,et)
               rv(i) = rre(nd,et)
               sv(i) = sre(nd,et)
             ENDDO

             
             i = 0             
             
       edge: DO vrt = 1,3
               n1 = mod(vrt+0,3)+1
               n2 = mod(vrt+1,3)+1
                              
               IF ((above(n1) == 1 .and. above(n2) == 0) .or. &
                   (above(n1) == 0 .and. above(n2) == 1)) THEN                                      
                   
                 r1 = rv(n1)
                 r2 = rv(n2)
                 s1 = sv(n1)
                 s2 = sv(n2)   
                 
                 dt = 1d0/real(ntry,rp)
                 t0 = 0d0
                 fail_flag = 1
        guesses: DO ig = 1,2*ntry+1
                   
                 t = t0                 
                 it = 0
         newton: DO 
         
                   r = .5d0*((1d0-t)*r1 + (1d0+t)*r2)
                   s = .5d0*((1d0-t)*s1 + (1d0+t)*s2)
                   re(1) = r
                   se(1) = s            
                   
                   IF (fig%type_flag == 2 .or. fig%type_flag == 4) THEN
#ifndef adcirc                   
                     CALL element_basis(et,hbp,ndf,1,re,se,phi,dpdr,dpds)                          
#else
                     CALL linear_basis(1,re,se,phi,dpdr,dpds)
                     ndf = 3
#endif                     
                     bathy = 0d0
                     dbdr = 0d0
                     dbds = 0d0
                     DO dof = 1,ndf
                       bathy = bathy + phi(dof,1)*hb(dof,el,1)
                       dbdr = dbdr + dpdr(dof,1)*hb(dof,el,1)
                       dbds = dbds + dpds(dof,1)*hb(dof,el,1)
                     ENDDO
                     IF (fig%type_flag == 2) THEN
                       f = bathy - sol_lev
                       dfdr = dbdr
                       dfds = dbds
                     ENDIF
                   ENDIF
                   
                   IF (fig%type_flag == 3 .or. fig%type_flag == 4) THEN
#ifndef adcirc
                     CALL element_basis(et,p,ndf,1,re,se,phi,dpdr,dpds)                          
#else                     
                     CALL linear_basis(1,re,se,phi,dpdr,dpds)
                     ndf = 3
#endif                     
                     zeta = 0d0
                     dzdr = 0d0
                     dzds = 0d0
                     DO dof = 1,ndf
                       zeta = zeta + phi(dof,1)*Z(dof,el,snap)
                       dzdr = dzdr + dpdr(dof,1)*Z(dof,el,snap)
                       dzds = dzds + dpds(dof,1)*Z(dof,el,snap)
                     ENDDO
                     IF (fig%type_flag == 3) THEN
                       f = zeta - sol_lev
                       dfdr = dzdr
                       dfds = dzds
                     ENDIF
                   ENDIF    
                   
                   IF (fig%type_flag == 4) THEN
                     xmom = 0d0
                     ymom = 0d0
                     dxmdr = 0d0
                     dxmds = 0d0
                     dymdr = 0d0
                     dymds = 0d0
                     DO dof = 1,ndf
                       xmom = xmom + phi(dof,1)*Qx(dof,el,snap)
                       ymom = ymom + phi(dof,1)*Qy(dof,el,snap)
                       dxmdr = dxmdr + dpdr(dof,1)*Qx(dof,el,snap)
                       dxmds = dxmds + dpds(dof,1)*Qx(dof,el,snap)
                       dymdr = dymdr + dpdr(dof,1)*Qy(dof,el,snap)
                       dymds = dymds + dpds(dof,1)*Qy(dof,el,snap)
                     ENDDO
                     H = zeta + bathy          
                     dHdr = dzdr + dbdr
                     dHds = dzds + dbds
                     uvel = xmom/H
                     vvel = ymom/H
                     f = uvel**2 + vvel**2 - sol_lev**2
                     dfdr = 2d0*(uvel*(dxmdr*H-dHdr*xmom)/H**2 + vvel*(dymdr*H-dHdr*ymom)/H**2)
                     dfds = 2d0*(uvel*(dxmds*H-dHds*xmom)/H**2 + vvel*(dymds*H-dHds*ymom)/H**2)
                   ENDIF
                   
                   
                   dfdt = .5d0*(dfdr*(r2-r1) + dfds*(s2-s1))
                   
                   t = t - f/dfdt
                   it = it + 1
                   
                   
                   IF (abs(f) < tol) THEN


                     IF (t <= 1d0+tol .and. t >= -(1d0+tol)) THEN
                       r = .5d0*((1d0-t)*r1 + (1d0+t)*r2)
                       s = .5d0*((1d0-t)*s1 + (1d0+t)*s2)                   
                       fail_flag = 0           
!                        PRINT*, "Iteration successful"
                       EXIT guesses
                     ELSE                 
!                        PRINT("(4(A,I9))"), "  Point outside edge     lev = ", lev, " el = ", el, " tri = ",tri, " vrt = ",vrt
!                        PRINT("(A,F12.6,A,I5,A,E24.17,A,E24.17)"), "      t0 = ", t0, "  it = ", it,"  t = ", t, "  f = ", f
!                        PRINT*, ""     
                       EXIT newton
                     ENDIF
                   ENDIF
                   
                   IF (it > max_it) THEN     
!                    PRINT*, "  Max iterations exceeded"
                     EXIT newton
                   ENDIF
                   
                 ENDDO newton
                 
                 t0 = t0*-1d0                 
                 IF (mod(ig,2) == 1) THEN
                   t0 = t0 + dt
                 ENDIF
                 
                 ENDDO guesses
                 
                 i = i + 1
                 IF (fail_flag == 0) THEN
                   re(1) = r
                   se(1) = s
                 ELSE                
                   PRINT*, "  iteration failed ", "r = ",r, " s = ",s 
                   re(1) = .5d0*(r1+r2)
                   se(1) = .5d0*(s1+s2)
                 ENDIF
 
                 CALL shape_functions_area_eval(et,np(et),nnd,1,re,se,l)
                 CALL element_transformation(nnd,elxy(:,el,1),elxy(:,el,2),l(:,1),x(i),y(i))   
                 
 
               !!!! Edge midpoint !!!!
!                  i = i + 1
!                  nd1 = rect(n1,tri,et)
!                  nd2 = rect(n2,tri,et)                 
!                  x(i) = 0.5d0*(xyplt(nd1,el,1) + xyplt(nd2,el,1))
!                  y(i) = 0.5d0*(xyplt(nd1,el,2) + xyplt(nd2,el,2))                   
                 
               ENDIF
             ENDDO edge

             IF (i > 2) THEN
               PRINT*, "Error: "
               STOP
             ENDIF 
             

             WRITE(file_unit,"(2(F9.5,1x))") x(1),y(1)
             WRITE(file_unit,"(2(F9.5,1x))") x(2),y(2)     
             WRITE(file_unit,"(A)") "draw-line"  

             
           ENDDO subtri
         ENDDO elem
         sol_lev = sol_lev + dc
      ENDDO levels

      END SUBROUTINE plot_line_contours

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE plot_mesh(file_unit,ne,nverts,ctp,el_type,el_in,xy,ect,xyplt)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: file_unit
      INTEGER, INTENT(IN) :: ne
      INTEGER, DIMENSION(:), INTENT(IN) :: nverts
      INTEGER, INTENT(IN) :: ctp
      INTEGER, DIMENSION(:), INTENT(IN) :: el_type
      INTEGER, DIMENSION(:), INTENT(IN) :: el_in      
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: xy
      INTEGER, DIMENSION(:,:), INTENT(IN) :: ect
      REAL(rp), DIMENSION(:,:,:), INTENT(IN) :: xyplt

      INTEGER :: el,nd
      INTEGER :: et

      
 elem:DO el = 1,ne
        
        IF (el_in(el) == 0) THEN
          CYCLE elem
        ENDIF
        
        et = el_type(el)
        
        IF (et == 1) THEN
          DO nd = 1,nverts(et)
            WRITE(file_unit,"(2(F9.5,1x))") xy(1,ect(nd,el)),xy(2,ect(nd,el))    
          ENDDO                
          WRITE(file_unit,"(A)") "draw-tri-element"         
        ELSE IF (et == 2) THEN
          DO nd = 1,nverts(et)
            WRITE(file_unit,"(2(F9.5,1x))") xy(1,ect(nd,el)),xy(2,ect(nd,el))    
          ENDDO                
          WRITE(file_unit,"(A)") "draw-quad-element"  
        ELSE
          WRITE(file_unit,"(A)") "newpath"
          WRITE(file_unit,"(2(F9.5,1x),A)") xyplt(1,el,1),xyplt(1,el,2),"moveto" 
          DO nd = 2,nverts(et)*ctp
            WRITE(file_unit,"(2(F9.5,1x),A)") xyplt(nd,el,1),xyplt(nd,el,2), "lineto"
          ENDDO
          WRITE(file_unit,"(A)") "closepath"
          WRITE(file_unit,"(A)") ".5 setlinewidth 2 setlinejoin"          
          WRITE(file_unit,"(A)") "stroke"          
        ENDIF
      ENDDO elem

      END SUBROUTINE plot_mesh

      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   

     SUBROUTINE plot_boundaries(file_unit,nverts,ctp,nbed,bedn,ged2el,ged2led,el_type,el_in,ect,xy,xyplt)
     
     IMPLICIT NONE
     
     INTEGER, INTENT(IN) :: file_unit
     INTEGER, DIMENSION(:), INTENT(IN) :: nverts
     INTEGER, INTENT(IN) :: ctp
     INTEGER, INTENT(IN) :: nbed
     INTEGER, DIMENSION(:), INTENT(IN) :: bedn
     INTEGER, DIMENSION(:,:), INTENT(IN) :: ged2el
     INTEGER, DIMENSION(:,:), INTENT(IN) :: ged2led
     INTEGER, DIMENSION(:), INTENT(IN) :: el_type
     INTEGER, DIMENSION(:), INTENT(IN) :: el_in
     INTEGER, DIMENSION(:,:), INTENT(IN) :: ect
     REAL(rp), DIMENSION(:,:), INTENT(IN) :: xy
     REAL(rp), DIMENSION(:,:,:), INTENT(IN) :: xyplt
     
     INTEGER :: ed,nd
     INTEGER :: el,et,nv
     INTEGER :: ged,led
     INTEGER :: n1,n2
     
edge:DO ed = 1,nbed
       ged = bedn(ed)
       
       el = ged2el(1,ged)
       led = ged2led(1,ged)
       
       et = el_type(el)
       nv = nverts(et)
       
       IF (el_in(el) == 0) THEN
         CYCLE
       ENDIF
       
       IF (et <= 2) THEN
       
          n1 = mod(led+0,nv)+1
          n2 = mod(led+1,nv)+1
          
          WRITE(file_unit,"(2(F9.5,1x))") xy(1,ect(n1,el)),xy(2,ect(n1,el))
          WRITE(file_unit,"(2(F9.5,1x))") xy(1,ect(n2,el)),xy(2,ect(n2,el)) 
          WRITE(file_unit,"(A)") "draw-line"          
        ELSE
          
          n1 = mod(led,nv)*ctp + 1
          n2 = n1 + ctp
          WRITE(file_unit,"(A)") "newpath"
          WRITE(file_unit,"(2(F9.5,1x),A)") xyplt(n1,el,1),xyplt(n1,el,2),"moveto" 
          DO nd = n1+1,n2
            WRITE(file_unit,"(2(F9.5,1x),A)") xyplt(nd,el,1),xyplt(nd,el,2), "lineto"
          ENDDO
          WRITE(file_unit,"(A)") ".5 setlinewidth 2 setlinejoin"          
          WRITE(file_unit,"(A)") "stroke"              
        ENDIF
       
     ENDDO edge
     
     END SUBROUTINE plot_boundaries
     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE fill_elements(file_unit,ne,nverts,ctp,el_type,el_in,xy,ect,xyplt)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: file_unit
      INTEGER, INTENT(IN) :: ne
      INTEGER, DIMENSION(:), INTENT(IN) :: nverts
      INTEGER, INTENT(IN) :: ctp
      INTEGER, DIMENSION(:), INTENT(IN) :: el_type
      INTEGER, DIMENSION(:), INTENT(IN) :: el_in      
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: xy
      INTEGER, DIMENSION(:,:), INTENT(IN) :: ect
      REAL(rp), DIMENSION(:,:,:), INTENT(IN) :: xyplt

      INTEGER :: el,nd
      INTEGER :: et
      INTEGER :: i,n
      INTEGER :: fill
      INTEGER, ALLOCATABLE, DIMENSION(:) :: fill_list
      LOGICAL :: file_exists
      
      ALLOCATE(fill_list(ne))
      fill_list = 0
      INQUIRE(file='element.fill',exist=file_exists)
      IF (file_exists) THEN
        OPEN(unit=101,file='element.fill')
        READ(101,*) n
        DO i = 1,n
          READ(101,*) el
          fill_list(el) = 1
        ENDDO
        CLOSE(101)
      ENDIF
      
 elem:DO el = 1,ne
        
        IF (el_in(el) == 0) THEN
          CYCLE elem
        ENDIF
        
        et = el_type(el)
        fill = fill_list(el)        
        
        IF (et == 1) THEN
          IF (fill == 1) THEN
            DO nd = 1,nverts(et)
              WRITE(file_unit,"(2(F9.5,1x))") xy(1,ect(nd,el)),xy(2,ect(nd,el))    
            ENDDO                
          WRITE(file_unit,"(A)") "fill-tri-element"            
          ENDIF
        ELSE IF (et == 2) THEN
          IF (fill == 1) THEN
            DO nd = 1,nverts(et)
              WRITE(file_unit,"(2(F9.5,1x))") xy(1,ect(nd,el)),xy(2,ect(nd,el))    
            ENDDO                
          WRITE(file_unit,"(A)") "fill-quad-element"            
          ENDIF
       ELSE
          IF (fill == 1) THEN       
            WRITE(file_unit,"(A)") "newpath"
            WRITE(file_unit,"(2(F9.5,1x),A)") xyplt(1,el,1),xyplt(1,el,2),"moveto" 
            DO nd = 2,nverts(et)*ctp
              WRITE(file_unit,"(2(F9.5,1x),A)") xyplt(nd,el,1),xyplt(nd,el,2), "lineto"
            ENDDO
            WRITE(file_unit,"(A)") "gsave 0 0 1 setrgbcolor fill grestore"
          ENDIF        
        ENDIF
      ENDDO elem

      END SUBROUTINE fill_elements     

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
      READ(101,*) ncolors
      ALLOCATE(colors(ncolors,3))
      DO lev = 1,ncolors
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

      SUBROUTINE setup_cbounds(fig,snap_start,snap_end)
      
      IMPLICIT NONE
      
      TYPE(plot_type), INTENT(INOUT) :: fig
      INTEGER, INTENT(IN)  :: snap_start,snap_end         
      
      INTEGER :: i  
      INTEGER :: start_snap,end_snap
      CHARACTER(:), ALLOCATABLE :: filename
      
      filename = TRIM(ADJUSTL(fig%name))//".cscale"      
      
      IF (fig%plot_sol_option == 1) THEN   
        OPEN(unit=fig%cscale_unit,file=filename//".out")  
        WRITE(fig%cscale_unit,"(2I5)") snap_start-1,snap_end-1
      ENDIF      
      
      IF (fig%cscale_option == "file") THEN
        OPEN(unit=fig%cscale_unit,file=filename)
        READ(fig%cscale_unit,*) start_snap,end_snap
        fig%num_cscale_vals = end_snap - start_snap + 1
        ALLOCATE(fig%cscale_vals(fig%num_cscale_vals,3))
        DO i = 1,fig%num_cscale_vals
          READ(fig%cscale_unit,*) fig%cscale_vals(i,1),fig%cscale_vals(i,2),fig%cscale_vals(i,3)
        ENDDO
        CLOSE(fig%cscale_unit)
      ENDIF      
      
      RETURN
      END SUBROUTINE setup_cbounds
      
            
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
!         PRINT*, command
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

       SUBROUTINE make_movie(fig,frmt)
       
       IMPLICIT NONE
       
       TYPE(plot_type), INTENT(IN) :: fig
       CHARACTER(*), INTENT(IN) :: frmt
       
        IF (fig%plot_sol_option == 1 .and. fig%movie_flag == 1) THEN
          CALL SYSTEM("ffmpeg -i "//TRIM(ADJUSTL(fig%name))//"_%04d."//frmt//" -y "//TRIM(ADJUSTL(fig%name))//".mp4")
        ENDIF          
       
       
       RETURN
       END SUBROUTINE make_movie

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      END MODULE plot_mod