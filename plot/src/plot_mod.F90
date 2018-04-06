      MODULE plot_mod
      
      USE globals, ONLY: rp,nel_type
      USE plot_globals, ONLY: plot_type,solution_type, &
                              cscale_width,lr_margin,top_margin,dash,fontsize, &
                              axes_width,axes_height, &
                              rmin_page,rmax_page,smin_page,smax_page, &
                              rmin_axes,rmax_axes,smin_axes,smax_axes, &
                              rmin_cbar,rmax_cbar,smin_cbar,smax_cbar, &
                              rmin_tbar,rmax_tbar,smin_tbar,smax_tbar, &
                              rmin_scale,rmax_scale,smin_scale,scale_loc, &                             
                              xticklabel_pad,yticklabel_pad,cticklabel_pad, &
                              xlabel_pad,ylabel_pad,clabel_pad,scale_pad, &
                              nxtick,nytick,nctick, &
                              nxdec,nydec,ncdec,ntdec, &
                              dr_xlabel,ds_ylabel,ds_clabel,main_font, &
                              ncolors,colors, &
                              mnpp,ps,pc,nord,adapt_option, &
                              phi_sol,ndof_sol, &
                              sol_diff_option,ho_diff_option, &
                              ax,ay,bx,by, &
                              xbox_min,xbox_max,ybox_min,ybox_max, &
                              snap_start,snap_end, &
                              max_init,min_init, &
                              p_low,p_high,p_skip                             
      
      IMPLICIT NONE
      
      INTEGER :: unit_count = 99   
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: elpt_diff
      REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: rspt_diff
      
      CONTAINS
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE make_plot(snap,t_snap,fig,sol1,sol2)
            
      USE globals, ONLY: nsta,xysta
      USE plot_globals, ONLY: t_start,t_end,xyplt,pplt,npplt,nptri,rect,r,s, &
                              frmt,density,pc,el_area, &
                              map,map_width,map_height,map_res, &
                              lamc,phic,plot_google_map,spherical_flag, &
                              region_box_option,region_box
      USE labels_mod, ONLY: latex_axes_labels,run_latex,read_latex, & 
                            latex_element_labels,latex_node_labels, &
                            write_latex_ps_body,remove_latex_files, &
                            write_char_array
      USE axes_mod, ONLY: write_all_axes
      USE google_map
            
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: snap
      REAL(rp), INTENT(IN) :: t_snap
      TYPE(plot_type), INTENT(INOUT) :: fig   
      TYPE(solution_type), INTENT(INOUT) :: sol1
      TYPE(solution_type), INTENT(INOUT) :: sol2

      INTEGER :: el,et,i
      INTEGER :: pe
      INTEGER, DIMENSION(:), ALLOCATABLE :: el_in_all
      CHARACTER(:), ALLOCATABLE :: filename
      CHARACTER(4) :: snap_char
      
      ! Return if solution is not being plotted
      IF (fig%plot_sol_option < 1) THEN
        RETURN
      ENDIF     
                
      PRINT("(A)"), "Ploting "//fig%name         
                
      ! Initialize element plotting order          
      IF (.not. ALLOCATED(fig%el_plt)) THEN
        ALLOCATE(fig%el_plt(sol1%ne))
      ENDIF                                 
      DO el = 1,sol1%ne
        et = sol1%el_type(el)
        fig%el_plt(el) = (et-1)*nord + nord
      ENDDO          
      
      ! Initialize things for max plot
      IF (fig%plot_max_option > 0 .and. snap > snap_end) THEN
        IF (fig%cscale_option /= "spec") THEN
          fig%cscale_vals(snap,2) = fig%min_maxval
          fig%cscale_vals(snap,3) = fig%max_maxval
        ENDIF
        fig%sol_label = "max "//TRIM(ADJUSTL(fig%sol_label))      
      ENDIF

      ! Set minimum and maximum value range for colorbar
      IF (fig%name /= "mesh")  THEN
        fig%sol_min = fig%cscale_vals(snap,2)
        fig%sol_max = fig%cscale_vals(snap,3)                                                            
        PRINT("(2(A,F10.5))"), "  min value = ", fig%sol_min, "  max value = ", fig%sol_max
        IF (fig%type_flag == 3 .or. fig%type_flag == 4) THEN
          WRITE(fig%cscale_unit,"(3(e24.17,1x))") t_snap,fig%sol_min,fig%sol_max              
        ENDIF
      ENDIF
       
      ! Create name for PS output file
      IF (fig%type_flag == 3 .or. fig%type_flag == 4) THEN
        IF (snap <= snap_end) THEN
          WRITE(snap_char,"(I4.4)") snap           
          filename = TRIM(ADJUSTL(fig%name))//"_"//snap_char  
        ELSE
          filename = TRIM(ADJUSTL(fig%name))//"_max"          
        ENDIF
      ELSE
        filename = TRIM(ADJUSTL(fig%name))
      ENDIF
      
      
      ! Initialize axes labels    
      CALL latex_axes_labels(fig,t_snap,t_start,t_end)  
      IF (fig%el_label_option /= "off") THEN
        CALL latex_element_labels(fig%el_label_option,sol1%ne,sol1%el_type,sol1%el_in,sol1%nverts,sol1%xy,sol1%ect,sol1%nnfbed,sol1%nfbedn,sol1%ged2el)
      ENDIF
      IF (fig%nd_label_option /= "off") THEN
        CALL latex_node_labels(fig%nd_label_option,sol1%nn,sol1%xy,sol1%nbou,sol1%fbseg,sol1%fbnds)
      ENDIF
      CALL run_latex(fig%tex_file_exists)      
      CALL read_latex(fig%tex_file_exists,fig%latex_header,fig%nline_header,fig%latex_body,fig%nline_body)
      CALL remove_latex_files()      
      
      ! Initialize PS file with PS functions, latex header, and backgound
      CALL write_psheader(filename//".ps",fig%ps_unit)  
      CALL write_char_array(fig%ps_unit,fig%nline_header,fig%latex_header)  
      CALL plot_background(fig%ps_unit,.75d0,.75d0,.75d0)
      
      ! Plot satellite image 
      IF (plot_google_map == 1 .and. spherical_flag == 1) THEN
        CALL write_map(fig%ps_unit,lamc,phic,xbox_min,xbox_max,ybox_min,ybox_max,ax,bx,ay,by,sol1%slam0,sol1%sphi0, &
                       map,map_height,map_width,map_res)
      ENDIF
      
      ! Plot solution
      IF (fig%type_flag > 1) THEN       
        CALL plot_filled_contours(fig%ps_unit,snap,fig,sol1,sol2)
!         CALL plot_filled_contours_adapt(fig%ps_unit,sol1,xyplt,snap,fig)             
      ENDIF
      IF (fig%plot_lines_option == 1) THEN
        CALL plot_line_contours(fig%ps_unit,sol1,sol1%nptri,sol1%rect,xyplt,sol1%r,sol1%s,snap,fig)          
      ENDIF      
      

      ! Plot adaptive visualization mesh
      IF (adapt_option == 1 .and. fig%type_flag > 1) THEN
        CALL SYSTEM("cp "//filename//".ps "//filename//"_pltmesh.ps")
        OPEN(UNIT=999,FILE=filename//"_pltmesh.ps",POSITION="APPEND")
        CALL plot_vis_mesh(999,sol1%ne,sol1%el_type,sol1%el_in,xyplt,sol1%nptri,sol1%rect,fig)         
        CALL write_all_axes(999,fig%axis_label_flag,fig%cbar_flag,fig%tbar_flag,t_snap,t_start,t_end)
!         CALL plot_elxy_nodes(999,sol1%ne,sol1%el_type,sol1%el_in,nnds,sol1%elxy)              
        CALL write_char_array(999,fig%nline_body,fig%latex_body)        
        CALL convert_ps(filename//"_pltmesh",frmt,density,fig%rm_ps)
        CLOSE(999)
      ENDIF            
      
      
      ! Plot mesh 
      IF (fig%name == "mesh") THEN
        CALL fill_elements(fig%ps_unit,sol1%ne,sol1%nverts,sol1%nnds,sol1%pplt,sol1%el_type,sol1%el_in,sol1%psic,sol1%elxy)        
      ENDIF        
      
      IF (fig%plot_mesh_option == 1) THEN              
       
        CALL plot_mesh(fig%ps_unit,sol1%ne,sol1%nverts,sol1%nnds,sol1%pplt,sol1%el_type,sol1%el_in,sol1%psic,sol1%elxy,sol1%mesh_line_color)
        IF (fig%sol_diff_option == 1) THEN
          CALL plot_mesh(fig%ps_unit,sol2%ne,sol2%nverts,sol2%nnds,sol2%pplt,sol2%el_type,sol2%el_in,sol2%psic,sol2%elxy,sol2%mesh_line_color)         
        ENDIF
!         CALL plot_nodestring(fig%ps_unit,sol1%nbou,sol1%nvel,sol1%fbseg,sol1%fbnds,sol1%xy)
!         CALL plot_qpts(fig%ps_unit,sol1%ne,sol1%p,sol1%ctp,sol1%np,sol1%el_type,sol1%el_in,sol1%elxy,sol1%nverts,sol1%ned,sol1%ed_type,sol1%ged2el,sol1%ged2led)

      ELSE IF (fig%plot_mesh_option == 2) THEN
            
        CALL plot_boundaries(fig%ps_unit,sol1%nverts,sol1%nnds,sol1%pplt,sol1%nobed,sol1%obedn,sol1%ged2el,sol1%ged2led,sol1%el_type,sol1%el_in,sol1%psic,sol1%elxy)       
        CALL plot_boundaries(fig%ps_unit,sol1%nverts,sol1%nnds,sol1%pplt,sol1%nnfbed,sol1%nfbedn,sol1%ged2el,sol1%ged2led,sol1%el_type,sol1%el_in,sol1%psic,sol1%elxy)        
        CALL plot_boundaries(fig%ps_unit,sol1%nverts,sol1%nnds,sol1%pplt,sol1%nfbed,sol1%fbedn,sol1%ged2el,sol1%ged2led,sol1%el_type,sol1%el_in,sol1%psic,sol1%elxy)          
      
        
      ELSE IF (fig%plot_mesh_option == 3) THEN
             
        CALL plot_decomp_mesh(fig%ps_unit,sol1%nverts)
        
      ENDIF 
      
      ! Plot region box
      IF (fig%name == "mesh" .and. region_box_option == 1) THEN
        CALL plot_region_box(fig%ps_unit,region_box)
      ENDIF
      
      ! Plot stations
      IF (fig%name == "mesh" .and. fig%plot_sta_option == 1) THEN
        CALL plot_stations(fig%ps_unit,fig%sta_start,fig%sta_end,xysta)        
      ENDIF        
      
!         ! Plot nodes      
!         CALL plot_cb_nodes(fig%ps_unit,sol1%ctp,sol1%nbou,sol1%fbseg,sol1%fbnds,sol1%xy,sol1%bndxy,sol1%nepn,sol1%epn,sol1%el_in)
!         CALL plot_elxy_nodes(fig%ps_unit,sol1%ne,sol1%el_type,sol1%el_in,sol1%nnds,sol1%elxy) 
!         CALL plot_qpts(fig%ps_unit,sol1%ne,sol1%p,sol1%ctp,sol1%np,sol1%el_type,sol1%el_in,sol1%elxy,sol1%nverts,sol1%ned,sol1%ed_type,sol1%ged2el,sol1%ged2led)
              
      
      ! Write axes and latexed labels
      CALL write_all_axes(fig%ps_unit,fig%axis_label_flag,fig%cbar_flag,fig%tbar_flag,t_snap,t_start,t_end)             
      CALL write_char_array(fig%ps_unit,fig%nline_body,fig%latex_body)       
      
      ! Finish up 
      CALL close_ps(filename,fig%ps_unit)
      CALL convert_ps(filename,frmt,density,fig%rm_ps)    
      
      ! Create separate color scales 
      IF (fig%name /= "mesh" .and. fig%cbar_flag == 0) THEN
        CALL make_cscale_horz(fig,filename)
!         CALL make_cscale_vert(fig,filename,t_snap,t_start,t_end)        
      ENDIF
      
      RETURN
      END SUBROUTINE make_plot


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      SUBROUTINE make_cscale_horz(fig,file)
      
      USE plot_globals, ONLY: frmt,density 
      USE labels_mod, ONLY: write_caxis_labels_horz, &
                            write_texheader,run_latex,read_latex, &                             
                            close_tex,remove_latex_files, &
                            write_char_array
      USE axes_mod, ONLY: write_colorscale_horz      
 
      IMPLICIT NONE
      
      TYPE(plot_type), INTENT(INOUT) :: fig   
      CHARACTER(*), INTENT(IN) :: file
      CHARACTER(:), ALLOCATABLE :: filename   
        
      filename = TRIM(ADJUSTL(file))//"_horz_cscale"         
      CALL write_texheader()
      CALL write_caxis_labels_horz(fig%sol_min,fig%sol_max,fig%sol_label)
      CALL close_tex()
      CALL run_latex(fig%tex_file_exists)
      CALL read_latex(fig%tex_file_exists,fig%latex_header,fig%nline_header,fig%latex_body,fig%nline_body)
      CALL remove_latex_files()   
      CALL write_psheader(filename//".ps",fig%ps_unit)
      CALL write_char_array(fig%ps_unit,fig%nline_header,fig%latex_header)   
      CALL write_colorscale_horz(fig%ps_unit)          
      CALL write_char_array(fig%ps_unit,fig%nline_body,fig%latex_body)       
      CALL close_ps(filename,fig%ps_unit)
      CALL convert_ps(filename,frmt,density,fig%rm_ps)
      
      RETURN
      END SUBROUTINE make_cscale_horz
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      SUBROUTINE make_cscale_vert(fig,file,t_snap,t_start,t_end)
      
      USE plot_globals, ONLY: frmt,density 
      USE labels_mod, ONLY: write_caxis_labels_horz,write_caxis_labels, &
                            write_texheader,run_latex,read_latex, &                             
                            close_tex,remove_latex_files, &
                            write_char_array,write_tbar_labels
      USE axes_mod, ONLY: write_colorscale_horz,write_colorscale,write_tbar      
 
      IMPLICIT NONE
      
      TYPE(plot_type), INTENT(INOUT) :: fig   
      CHARACTER(*), INTENT(IN) :: file
      CHARACTER(:), ALLOCATABLE :: filename   
      REAL(rp), INTENT(IN) :: t_snap
      REAL(rp), INTENT(IN) :: t_start
      REAL(rp), INTENT(IN) :: t_end      
        
      filename = TRIM(ADJUSTL(file))//"_vert_cscale"         
      CALL write_texheader()
      CALL write_caxis_labels(fig%tbar_flag,fig%sol_min,fig%sol_max,fig%sol_label)
      IF (fig%tbar_flag == 1) THEN
        CALL write_tbar_labels(t_snap)      
      ENDIF 
      CALL close_tex()
      CALL run_latex(fig%tex_file_exists)
      CALL read_latex(fig%tex_file_exists,fig%latex_header,fig%nline_header,fig%latex_body,fig%nline_body)
      CALL remove_latex_files()   
      CALL write_psheader(filename//".ps",fig%ps_unit)
      CALL write_char_array(fig%ps_unit,fig%nline_header,fig%latex_header)    
      CALL write_colorscale(fig%ps_unit)    
      IF (fig%tbar_flag == 1) THEN
        CALL write_tbar(fig%ps_unit,t_snap,t_start,t_end)      
      ENDIF
      CALL write_char_array(fig%ps_unit,fig%nline_body,fig%latex_body)       
      CALL close_ps(filename,fig%ps_unit)
      CALL convert_ps(filename,frmt,density,fig%rm_ps)
      
      RETURN
      END SUBROUTINE make_cscale_vert      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE prep_station_plot(fig,axis_label_option,snap,xmin,xmax,ymin,ymax,filename)
      
      USE plot_globals, ONLY: figure_width,figure_height,frmt,density, &
                              ax,bx,ay,by
      USE labels_mod, ONLY: write_texheader,run_latex,read_latex, & 
                            write_xyaxis_labels,remove_latex_files, &
                            write_char_array
      USE axes_mod, ONLY: write_xyaxis                             
      
      IMPLICIT NONE
            
      TYPE(plot_type), INTENT(INOUT) :: fig   
      CHARACTER(*), INTENT(IN) :: axis_label_option
      INTEGER, INTENT(IN) :: snap
      REAL(rp), INTENT(IN) :: xmax,xmin,ymin,ymax
      
      CHARACTER(:), ALLOCATABLE, INTENT(IN) :: filename     

           
      CALL scale_factors_line_plot(figure_width,figure_height,xmin,xmax,ymin,ymax,ax,bx,ay,by)      
      
      CALL write_texheader()  
      CALL write_xyaxis_labels(axis_label_option)      
      
      CALL run_latex(fig%tex_file_exists)
      CALL read_latex(fig%tex_file_exists,fig%latex_header,fig%nline_header,fig%latex_body,fig%nline_body)
      CALL remove_latex_files()     
        
      CALL write_psheader(filename//".ps",fig%ps_unit)  
      CALL write_char_array(fig%ps_unit,fig%nline_header,fig%latex_header)  
      CALL plot_background(fig%ps_unit,1d0,1d0,1d0)      
               
      
      RETURN
      END SUBROUTINE prep_station_plot
      


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  


      SUBROUTINE finish_station_plot(fig,filename)
      
      USE plot_globals, ONLY: frmt,density
      USE labels_mod, ONLY: write_char_array
      USE axes_mod, ONLY: write_xyaxis           
      
      IMPLICIT NONE
      
      TYPE(plot_type), INTENT(INOUT) :: fig      
      CHARACTER(*), INTENT(IN) :: filename
      
      CALL write_xyaxis(fig%ps_unit,1)       
      CALL write_char_array(fig%ps_unit,fig%nline_body,fig%latex_body)  
      
      CALL close_ps(filename,fig%ps_unit)
      CALL convert_ps(filename,frmt,density,fig%rm_ps)            
                  
      RETURN
      END SUBROUTINE finish_station_plot
      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      SUBROUTINE zoom_box(ne,npplt,nnds,el_type,psic,elxy,xbox_min,xbox_max,ybox_min,ybox_max,xmin,xmax,ymin,ymax,el_in)
      
      USE transformation, ONLY: element_transformation
      USE basis, ONLY: element_nodes             
      USE shape_functions_mod, ONLY: shape_functions_area_eval
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: ne
      INTEGER, DIMENSION(:), INTENT(IN) :: npplt    
      INTEGER, DIMENSION(:), INTENT(IN) :: nnds
      INTEGER, DIMENSION(:), INTENT(IN) :: el_type
      REAL(rp), DIMENSION(:,:,:), INTENT(IN) :: psic
      REAL(rp), DIMENSION(:,:,:), INTENT(IN) :: elxy
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
      INTEGER :: i,et,npts,nnd
      INTEGER, DIMENSION(:), ALLOCATABLE :: outside
      INTEGER :: in_flag
      REAL(rp) :: xpt,ypt     
         
      
      xmax = max_init
      ymax = max_init
      xmin = min_init
      ymin = min_init      
      
      ALLOCATE(outside(maxval(npplt)))
      
      ALLOCATE(el_in(ne))
      el_in = 1
      
      DO el = 1,ne      
        et = el_type(el)                
        npts = npplt(et)
        nnd = nnds(et)
        outside = 0
        DO pt = 1,npts  
        
          CALL element_transformation(nnd,elxy(:,el,1),elxy(:,el,2),psic(:,pt,et),xpt,ypt)             
                
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

      SUBROUTINE scale_factors(figure_width,figure_height,xmin,xmax,ymin,ymax,ax,bx,ay,by)
      
      IMPLICIT NONE
      
      REAL(rp), INTENT(IN) :: figure_width
      REAL(rp), INTENT(IN) :: figure_height
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
               
      lr_margin = (rmax_page - figure_width)/2d0 
      axes_width = figure_width/1.37d0
      
      cscale_width = axes_width*.05d0
      
      dash = 0.2d0*cscale_width
      xticklabel_pad = 2d0*dash
      yticklabel_pad = 2d0*dash
      xlabel_pad = 10d0*dash
      ylabel_pad = 15d0*dash        
      cticklabel_pad = 2d0*dash
      clabel_pad = 15d0*dash 
      scale_pad = 5d0*dash
      

         
      
      
      rmin_axes = rmin_page + lr_margin + ylabel_pad
      rmax_axes = rmax_page - lr_margin - 2d0*cscale_width - clabel_pad
      smin_axes = smin_page + lr_margin
      smax_axes = smax_page      
      
                
      ax = ((rmax_axes-rmin_axes)/(xmax-xmin))
      bx = (rmin_axes*xmax-rmax_axes*xmin)/(xmax-xmin)
      
      IF (figure_height < 0d0) THEN
        ay = ax     ! axis equal            
      ELSE
!         axes_height = (rmax_axes-rmin_axes)/1.25d0

        axes_height = figure_height
        smax_axes = smin_axes+axes_height
        ay = ((smax_axes-smin_axes)/(ymax-ymin))      
      ENDIF
      by = smin_axes-ay*ymin              
      
      smax_axes = ay*ymax + by
      
      IF (smax_axes > smax_page) THEN
        PRINT*, "Error: Equal axis scaling puts y axis off page"
        STOP
      ENDIF
      
      rmin_tbar = rmax_axes + cscale_width
      rmax_tbar = rmin_tbar + clabel_pad
      smin_tbar = smax_axes - 2d0*cscale_width
      smax_tbar = smin_tbar + cscale_width   
      
      IF (scale_loc == "SW") THEN
        rmin_scale = rmin_axes + scale_pad !3*scale_pad
        smin_scale = smin_axes + scale_pad
      ELSE IF (scale_loc == "W") THEN
        rmin_scale = rmin_axes + scale_pad
        smin_scale = .5d0*(smax_axes+smin_axes)
      ELSE IF (scale_loc == "NW") THEN
        rmin_scale = rmin_axes + scale_pad
        smin_scale = smax_axes - scale_pad      
      ELSE IF (scale_loc == "NE") THEN
        rmax_scale = rmax_axes - scale_pad
        smin_scale = smax_axes - scale_pad      
      ELSE IF (scale_loc == "E") THEN
        rmax_scale = rmax_axes - scale_pad
        smin_scale = .5d0*(smax_axes+smin_axes)      
      ELSE IF (scale_loc == "SE") THEN
        rmax_scale = rmax_axes - scale_pad
        smin_scale = smin_axes + scale_pad      
      ENDIF
      
      RETURN
      END SUBROUTINE scale_factors
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      SUBROUTINE scale_factors_line_plot(figure_width,figure_height,xmin,xmax,ymin,ymax,ax,bx,ay,by)
      
      IMPLICIT NONE
      
      REAL(rp), INTENT(IN) :: figure_width
      REAL(rp), INTENT(IN) :: figure_height
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
               
      lr_margin = (rmax_page - figure_width)/2d0 
      axes_width = figure_width
      
      cscale_width = axes_width*.05d0
      
      dash = 0.2d0*cscale_width
      xticklabel_pad = 2d0*dash
      yticklabel_pad = 2d0*dash
      xlabel_pad = 7d0*dash
      ylabel_pad = 12d0*dash        
      

      rmin_axes = rmin_page + lr_margin + ylabel_pad
      rmax_axes = rmax_page - lr_margin 
      smin_axes = smin_page + lr_margin + xlabel_pad
      smax_axes = smax_page      
      
    
            
      ax = ((rmax_axes-rmin_axes)/(xmax-xmin))
      bx = (rmin_axes*xmax-rmax_axes*xmin)/(xmax-xmin)
      
      IF (figure_height < 0d0) THEN
        ay = ax     ! axis equal            
      ELSE
        axes_height = figure_height
        smax_axes = smin_axes+axes_height
        ay = ((smax_axes-smin_axes)/(ymax-ymin))      
      ENDIF
      by = smin_axes-ay*ymin              
      
      smax_axes = ay*ymax + by
      
      IF (smax_axes > smax_page) THEN
        PRINT*, "Error: Equal axis scaling puts y axis off page"
        STOP
      ENDIF      
      
      RETURN
      END SUBROUTINE scale_factors_line_plot      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
      
      
      SUBROUTINE write_psheader(file_name,file_unit)
      
      USE labels_mod, ONLY: write_latex_ps_header
      
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
      
      WRITE(file_unit,"(A)") "/draw-tri-element-color {"
      WRITE(file_unit,"(A)") "newpath"        
      WRITE(file_unit,"(A)") "moveto"
      WRITE(file_unit,"(A)") "lineto"
      WRITE(file_unit,"(A)") "lineto"
      WRITE(file_unit,"(A)") "closepath"
      WRITE(file_unit,"(A)") ".5 setlinewidth 2 setlinejoin"
      WRITE(file_unit,"(A)") "gsave setrgbcolor stroke grestore"           
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
      
      WRITE(file_unit,"(A)") "/draw-quad-element-color {"
      WRITE(file_unit,"(A)") "newpath"        
      WRITE(file_unit,"(A)") "moveto"
      WRITE(file_unit,"(A)") "lineto"
      WRITE(file_unit,"(A)") "lineto"
      WRITE(file_unit,"(A)") "lineto"      
      WRITE(file_unit,"(A)") "closepath"
      WRITE(file_unit,"(A)") ".5 setlinewidth 2 setlinejoin"
      WRITE(file_unit,"(A)") "gsave setrgbcolor stroke grestore"
      WRITE(file_unit,"(A)") "} def"         
      
      WRITE(file_unit,"(A)") "/fill-tri-element {"
      WRITE(file_unit,"(A)") "newpath"        
      WRITE(file_unit,"(A)") "moveto"
      WRITE(file_unit,"(A)") "lineto"
      WRITE(file_unit,"(A)") "lineto"
      WRITE(file_unit,"(A)") "closepath"
      WRITE(file_unit,"(A)") "3 copy"            
      WRITE(file_unit,"(A)") "gsave setrgbcolor fill grestore"     
      WRITE(file_unit,"(A)") "gsave setrgbcolor 0 setlinewidth stroke grestore"        
      WRITE(file_unit,"(A)") "} def" 
      
      WRITE(file_unit,"(A)") "/fill-quad-element {"
      WRITE(file_unit,"(A)") "newpath"        
      WRITE(file_unit,"(A)") "moveto"
      WRITE(file_unit,"(A)") "lineto"
      WRITE(file_unit,"(A)") "lineto"
      WRITE(file_unit,"(A)") "lineto"      
      WRITE(file_unit,"(A)") "closepath"   
      WRITE(file_unit,"(A)") "3 copy"            
      WRITE(file_unit,"(A)") "gsave setrgbcolor fill grestore"         
      WRITE(file_unit,"(A)") "gsave setrgbcolor 0 setlinewidth stroke grestore"         
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
      
      WRITE(file_unit,"(A)") "/draw-dot {"
      WRITE(file_unit,"(A)") "newpath"      
      WRITE(file_unit,"(A)") "0 360 arc closepath"
!       WRITE(file_unit,"(A)") "stroke"
!       WRITE(file_unit,"(A)") "gsave 0 0 0 setrgbcolor fill grestore"
      WRITE(file_unit,"(A)") "gsave setrgbcolor fill grestore"
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
      
      WRITE(file_unit,"(A)") "/fill-pixel {"
      WRITE(file_unit,"(A)") "newpath"        
      WRITE(file_unit,"(A)") "moveto"
      WRITE(file_unit,"(A)") "lineto"
      WRITE(file_unit,"(A)") "lineto"
      WRITE(file_unit,"(A)") "lineto"      
      WRITE(file_unit,"(A)") "closepath"
      WRITE(file_unit,"(A)") "3 copy"      
      WRITE(file_unit,"(A)") "gsave setrgbcolor fill grestore"  
      WRITE(file_unit,"(A)") "gsave setrgbcolor 0 setlinewidth stroke grestore"              
      WRITE(file_unit,"(A)") "} def"         
      
      WRITE(file_unit,"(A)") "/region-box{"      
      WRITE(file_unit,"(A)") "newpath" 
      WRITE(file_unit,"(A)") "moveto" 
      WRITE(file_unit,"(A)") "lineto" 
      WRITE(file_unit,"(A)") "lineto" 
      WRITE(file_unit,"(A)") "lineto"
      WRITE(file_unit,"(A)") "closepath"   
      WRITE(file_unit,"(A)") "gsave 3 setlinewidth 2 setlinejoin 0 0 0 setrgbcolor stroke grestore"         
!       WRITE(file_unit,"(A)") "gsave 3 setlinewidth 2 setlinejoin 0.81 0.1 0.11 setrgbcolor stroke grestore"      
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
      
            
!       CALL write_latex_ps_header(file_unit) 

     
 
  
      RETURN
      END SUBROUTINE write_psheader
      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      

      SUBROUTINE plot_background(file_unit,r,b,g)
     
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: file_unit
      REAL(rp), INTENT(IN) :: r,b,g
     
      WRITE(file_unit,"(A)") "gsave newpath"        
      WRITE(file_unit,"(2(F9.5,1x),A)") rmin_axes,smin_axes,"moveto"
      WRITE(file_unit,"(2(F9.5,1x),A)") rmax_axes,smin_axes,"lineto"
      WRITE(file_unit,"(2(F9.5,1x),A)") rmax_axes,smax_axes,"lineto"      
      WRITE(file_unit,"(2(F9.5,1x),A)") rmin_axes,smax_axes,"lineto"      
      WRITE(file_unit,"(A)") "closepath"     
      WRITE(file_unit,"(3(F9.5,1x),A)") r,g,b," setrgbcolor fill grestore"          
     
      RETURN
      END SUBROUTINE plot_background


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      
      SUBROUTINE close_ps(filename,unit_number)
      
      IMPLICIT NONE
      
      CHARACTER(*), INTENT(IN) :: filename
      INTEGER, INTENT(IN) :: unit_number                    
          
!       WRITE(unit_number,"(A)") "showpage"
      CLOSE(unit_number)            
!       CALL SYSTEM("ps2eps -B -C < "//filename//".ps > "//filename//".eps")
      
      RETURN
      END SUBROUTINE close_ps        
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE plot_filled_contours(file_unit,snap,fig,sol1,sol2)
     
      USE plot_globals, ONLY: xyplt
      USE evaluate_mod, ONLY: evaluate_solution         
     
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: file_unit
      INTEGER, INTENT(IN) :: snap
      TYPE(plot_type), INTENT(INOUT) :: fig
      TYPE(solution_type), INTENT(IN) :: sol1      
      TYPE(solution_type), INTENT(IN) :: sol2
      
      INTEGER :: el,pt
      INTEGER :: et,i,nnd
      INTEGER :: elin
      INTEGER :: sol_snap
      REAL(rp) :: xy(2)
      REAL(rp) :: rs(2),r(1),s(1)
      REAL(rp) :: sol2_val(mnpp),val(1)      
      
      
      
      
      ! Allocate and initialize max array (if needed)
      IF (fig%plot_max_option > 0) THEN
        IF (.not. ALLOCATED(fig%sol_maxval)) THEN
          ALLOCATE(fig%sol_maxval(mnpp,sol1%ne))
          fig%sol_maxval = max_init
        ENDIF
        fig%max_maxval = max_init
        fig%min_maxval = min_init        
      ENDIF 
      
      ! Find the solution 2 elements each plotting node lies in and calculate local r,s coordinates
      IF (sol_diff_option == 1) THEN      
        IF (.not. ALLOCATED(elpt_diff)) THEN
          CALL find_diff_nodes(fig,sol1,sol2)
        ENDIF        
      ENDIF
           
      
      
 elem:DO el = 1,sol1%ne
      
        ! Skip elements not in the specified zoom box
        IF (sol1%el_in(el) == 0) THEN
          CYCLE elem
        ENDIF                             
        
        ! Determine element plotting node order
        et = sol1%el_type(el)            
                              
        IF (snap <= snap_end) THEN
        
          ! Evaluate solution 1 at plotting nodes, using pre-calculated basis functions         
          CALL evaluate_solution(el,fig%type_flag,snap,sol1,fig%sol_val(:,el),sol1%npplt(et),phi_sol=phi_sol(:,:,et))         

          
          
        
          ! Evaluate solution 1 at plotting nodes, truncating high-order terms (if needed)
          IF (fig%ho_diff_option == 1) THEN
            CALL evaluate_solution(el,fig%type_flag,snap,sol1,sol2_val,sol1%npplt(et),phi_sol=phi_sol(:,:,et),plim=fig%plim)        
          ENDIF
        
          ! Evaluate solution 2 at plotting nodes (if needed)
          IF (fig%sol_diff_option == 1) THEN
            DO pt = 1,sol1%npplt(et)              
              elin = elpt_diff(pt,el)    
              r(1) = rspt_diff(1,pt,el)
              s(1) = rspt_diff(2,pt,el) 
          
              CALL evaluate_solution(elin,fig%type_flag,snap,sol2,val,1,r,s)
              sol2_val(pt) = val(1)
            ENDDO                                     
          ENDIF
        
          ! Calculate solution difference (if_needed)
          IF (fig%sol_diff_option == 1 .or. fig%ho_diff_option == 1) THEN
            DO pt = 1,sol1%npplt(et)
              fig%sol_val(pt,el) = abs(fig%sol_val(pt,el) - sol2_val(pt))        
            ENDDO          
          ENDIF
          
          ! Determine maximums (if needed)       
          IF (fig%plot_max_option > 0) THEN
           DO pt = 1,sol1%npplt(et)
              IF (fig%sol_val(pt,el) > fig%sol_maxval(pt,el)) THEN
                fig%sol_maxval(pt,el) = fig%sol_val(pt,el)
              ENDIF   
              IF (fig%sol_maxval(pt,el) > fig%max_maxval) THEN
                fig%max_maxval = fig%sol_maxval(pt,el)
              ENDIF
              IF (fig%sol_maxval(pt,el) < fig%min_maxval .and. fig%sol_val(pt,el) > max_init) THEN
                fig%min_maxval = fig%sol_maxval(pt,el)
              ENDIF
            ENDDO
          ENDIF          
          
        ELSE 
        
          ! Set the plotting array to max values (when called with: snap = snap_end+1)
          IF (fig%plot_max_option > 0) THEN
              DO pt = 1,sol1%npplt(et)
                fig%sol_val(pt,el) = fig%sol_maxval(pt,el)
              ENDDO
          ENDIF        
          
        ENDIF
       
        
      
        ! Write the PS code to contour fill the element
        IF (file_unit > 0) THEN
          CALL contour_fill_element(file_unit,sol1%nptri(et),sol1%rect(:,:,et),fig%sol_min,fig%sol_max,fig%sol_val(:,el),xyplt(:,el,1),xyplt(:,el,2))
        ENDIF
        
      ENDDO elem     
     
      RETURN
      END SUBROUTINE plot_filled_contours
      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE plot_filled_contours_adapt(file_unit,sol,xyplt,snap,fig)
      
      USE transformation, ONLY: init_vandermonde,element_transformation,xy2rs
      USE basis, ONLY: element_basis,linear_basis              
      USE area_qpts_mod, ONLY: tri_cubature,area_qpts
      USE shape_functions_mod, ONLY: shape_functions_area_eval
      USE evaluate_mod, ONLY: evaluate_solution,evaluate_plotting_nodes
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: file_unit     
      TYPE(solution_type), INTENT(IN) :: sol      
      REAL(rp), DIMENSION(:,:,:), INTENT(INOUT) :: xyplt 
      INTEGER, INTENT(IN) :: snap      
      TYPE(plot_type), INTENT(INOUT) :: fig

      
      INTEGER :: i,j,v,k
      INTEGER :: el,nd,dof,lev,tri,ord,pt
      INTEGER :: et,nv,pl    
      INTEGER :: nqpt,nnd,mnnds,mnp,ndf,mnqpta,nqpta(nel_type)
      INTEGER :: qpt_order
      INTEGER :: nptri_total
      INTEGER :: ne_total
      INTEGER :: pplt_max
      INTEGER :: nord

      INTEGER, DIMENSION(:), ALLOCATABLE :: pplt 
      INTEGER, DIMENSION(:), ALLOCATABLE :: nptri
      INTEGER, DIMENSION(:), ALLOCATABLE :: npplt      
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: rect 
      REAL(rp) :: sol_lev    
      REAL(rp), DIMENSION(:,:), ALLOCATABLE :: qpt
      REAL(rp), DIMENSION(:), ALLOCATABLE :: wpt
      REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: qpta
      REAL(rp), DIMENSION(:,:), ALLOCATABLE ::  wpta
      REAL(rp), DIMENSION(:,:), ALLOCATABLE :: l
      REAL(rp), DIMENSION(:), ALLOCATABLE :: xpt,ypt
      REAL(rp), DIMENSION(:), ALLOCATABLE :: rpt,spt
      REAL(rp), DIMENSION(:), ALLOCATABLE :: sol_lin,sol_el  
      REAL(rp), DIMENSION(:,:), ALLOCATABLE :: r,s
      REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: psic
      REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: phia
      REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: psi,dpsidr,dpsids      
      REAL(rp) :: xpta,ypta
      REAL(rp) :: drdx,drdy,dsdx,dsdy      
      REAL(rp) :: qpts(13**2,3)
      REAL(rp) :: xv(3),yv(3)
      REAL(rp) :: detJ
      REAL(rp) :: H,u_vel,v_vel
      REAL(rp) :: sol_avg
      REAL(rp) :: error_total
   
      REAL(rp) :: color_val(3)
      REAL(rp) :: err,max_err
      
      INTEGER :: elcnt,ndcnt
      REAL(rp) :: el_max
      INTEGER, DIMENSION(:), ALLOCATABLE :: el_list,nd_list,nd_flag
      
      err = 0d0
                  
      ! Get plotting nodes and triangularization for all nodal set orders
      CALL evaluate_plotting_nodes(sol%nel_type,p_low,p_skip,nord,mnpp,sol%np,sol%mnnds, &
                                           pplt,npplt,r,s,psic,nptri,rect)
      
      ! quadrature points for elemental solution average
      CALL area_qpts(1,sol%p,sol%ctp,nel_type,nqpta,mnqpta,wpta,qpta)      
      ALLOCATE(phia((sol%p+1)**2,mnqpta,nel_type))

      ! quadrature points for error integral
      qpt_order = 2
      CALL tri_cubature(qpt_order,nqpt,qpts)      
      mnqpta = max(mnqpta,nqpt)
         
      ALLOCATE(qpt(mnqpta,2),wpt(mnqpta))
      ALLOCATE(l(3,mnqpta))
      ALLOCATE(xpt(mnqpta),ypt(mnqpta))
      ALLOCATE(rpt(mnqpta),spt(mnqpta))
      ALLOCATE(sol_lin(mnqpta),sol_el(mnqpta))    
      
      DO pt = 1,nqpt
        qpt(pt,1) = qpts(pt,1)
        qpt(pt,2) = qpts(pt,2)
        wpt(pt) = qpts(pt,3)
      ENDDO
      
      
      
      ! shape functions for elemental solution average
      mnp = maxval(sol%np)      
      mnnds = (mnp+1)**2     
      ALLOCATE(psi(mnnds,mnqpta,nel_type),dpsidr(mnnds,mnqpta,nel_type),dpsids(mnnds,mnqpta,nel_type))
      
      DO et = 1,nel_type         
        CALL shape_functions_area_eval(et,sol%np(et),nnd,nqpta(et),qpta(:,1,et),qpta(:,2,et), &
                                       psi(:,:,et),dpsidr(:,:,et),dpsids(:,:,et))
      ENDDO      
      
      
      
       
      
      CALL init_vandermonde(nel_type,sol%np)      
            

            
      nptri_total = 0
      error_total = 0d0
      ne_total = 0
      pplt_max = 0
      
 elem:DO el = 1,sol%ne

        IF (sol%el_in(el) == 0) THEN
          CYCLE elem
        ENDIF
        
        et = sol%el_type(el)       
        nnd = sol%nnds(et)        
          
        IF (et > 2) THEN
         pl = pc
        ELSE
          pl = 1
        ENDIF
        
        IF (adapt_option == 1) THEN           
          ! Find element average value for solution
          CALL evaluate_solution(el,fig%type_flag,snap,sol,sol_el,nqpta(et),qpta(:,1,et),qpta(:,2,et))
        
          sol_avg = 0d0
          DO pt = 1,nqpta(et)
            CALL element_transformation(nnd,sol%elxy(:,el,1),sol%elxy(:,el,2),psi(:,pt,et),xpta,ypta, &
                                        dpsidr(:,pt,et),dpsids(:,pt,et),drdx,drdy,dsdx,dsdy,detJ)                
            sol_avg = sol_avg + wpta(pt,et)*sol_el(pt)*detJ
          ENDDO
          sol_avg = sol_avg/sol%el_area(el)
             
          PRINT*, ""        
          PRINT "(A,I9,A,I9,A,F15.7)", "ELEMENT: ", el, "    et: ", et, "    sol_avg = ", sol_avg        
        ENDIF
         
 order: DO ord = pl,nord
 
          i = (et-1)*nord+ord
          

          DO pt = 1,npplt(i)              
            CALL element_transformation(nnd,sol%elxy(:,el,1),sol%elxy(:,el,2),psic(:,pt,i),xyplt(pt,el,1),xyplt(pt,el,2))           
          ENDDO                  
     
          ! Evaluate solution at plotting nodes              
          CALL evaluate_solution(el,fig%type_flag,snap,sol,fig%sol_val(:,el),npplt(i),phi_sol=phi_sol(:,:,i))           

          IF (adapt_option == 0) THEN
            EXIT order
          ENDIF
    
          
          CALL linear_basis(ndf,nqpt,qpt(:,1),qpt(:,2),l)          
          
          err = 0d0
          max_err = 0d0
  sub_tri:DO tri = 1,nptri(i)
          
            ! Evaluate linear solution at plotting element quadrature points  
            sol_lin = 0d0
            DO pt = 1,nqpt
              DO v = 1,3              
                 nd = rect(v,tri,i)
                sol_lin(pt) = sol_lin(pt) + l(v,pt)*fig%sol_val(nd,el)
              ENDDO
            ENDDO
            
            DO pt = 1,nqpt
              xpt(pt) = 0d0
              ypt(pt) = 0d0
              DO v = 1,3
                nd = rect(v,tri,i)
                
                xv(v) = xyplt(nd,el,1)
                yv(v) = xyplt(nd,el,2)
                
                xpt(pt) = xpt(pt) + l(v,pt)*xv(v)
                ypt(pt) = ypt(pt) + l(v,pt)*yv(v)
              ENDDO
            ENDDO
            
            detJ = .25d0*((xv(2)-xv(1))*(yv(3)-yv(1))-(xv(3)-xv(1))*(yv(2)-yv(1)))
            
            CALL xy2rs(et,sol%np(et),sol%elxy(:,el,1),sol%elxy(:,el,2),nqpt,xpt,ypt,rpt,spt)   
            
            ! Evaluate DG solution at plotting element quadrature points
            CALL evaluate_solution(el,fig%type_flag,snap,sol,sol_el,nqpt,rpt,spt)

            
            
            DO pt = 1,nqpt                        
              err = err + detJ*wpt(pt)*(sol_el(pt)-sol_lin(pt))**2
            ENDDO
 
              max_err = sqrt(err)

!             DO pt = 1,nqpt               
!               err = abs(sol_el(pt)-sol_lin(pt))
!               IF (err > max_err) THEN
!                 max_err = err
!               ENDIF
!             ENDDO
            
            IF (abs(max_err/sol_avg) > fig%rel_tol .and. max_err > fig%abs_tol) THEN
              PRINT "(A,I4,A,I4)", "    Early sub_tri exit: ", tri, "/", nptri(i)              
              EXIT sub_tri
            ENDIF
            
          ENDDO sub_tri
          
          err = sqrt(err)
!           err = max_err
          
          PRINT "(A,I9,A,F15.7,A,F15.7)", "  ORDER: ", pplt(ord), "  Abs. Error = ", err, "  Rel. Error = ", abs(err/sol_avg)         
          
          IF (abs(err/sol_avg) < fig%rel_tol .or. err < fig%abs_tol) THEN
            fig%el_plt(el) = i       
            EXIT order     
          ENDIF          
          
        ENDDO order                 
          
                
        
        ord = fig%el_plt(el)
        
        
         
        CALL contour_fill_element(file_unit,nptri(ord),rect(:,:,ord),fig%sol_min,fig%sol_max,fig%sol_val(:,el),xyplt(:,el,1),xyplt(:,el,2))

        
        
        error_total = error_total + err**2
        nptri_total = nptri_total + nptri(ord)
        ne_total = ne_total + 1
        
        IF (pplt(ord) > pplt_max) THEN
          pplt_max = pplt(ord)
        ENDIF

      ENDDO elem
      
      IF (adapt_option == 1) THEN
        WRITE(998,"(A25,I25,E24.17,1x,I25,I25,I25)") fig%name, snap, sqrt(error_total), nptri_total, pplt_max, ne_total
      ENDIF
      
      
!       ALLOCATE(el_list(sol%ne),nd_list(sol%ne),nd_flag(sol%ne))
!       nd_flag = 0
!       elcnt = 0
!       ndcnt = 0
!       DO el = 1,sol%ne
!       
!         IF (sol%el_in(el) == 0) THEN
!           CYCLE
!         ENDIF      
!       
!         i = fig%el_plt(el)
!         et = sol%el_type(el)
!         nv = sol%nverts(et)
!         
!         el_max = -1d99
!         DO pt = 1,npplt(i)
!           IF (fig%sol_val(pt,el) > el_max) THEN
!             el_max = fig%sol_val(pt,el)
!           ENDIF          
!         ENDDO
!         
!         IF (el_max < 0.0004) THEN
!           elcnt = elcnt + 1
!           el_list(elcnt) = el
!           DO j = 1,nv
!             nd = sol%ect(j,el)
!             IF (nd_flag(nd) == 0)THEN
!               ndcnt = ndcnt + 1
!               nd_list(ndcnt) = nd
!               nd_flag(nd) = 1
!             ENDIF
!           ENDDO
!         ENDIF      
!         
!         
!       ENDDO
!       
!       OPEN(unit=123,file="element.fill")
!       WRITE(123,*) elcnt
!       DO i = 1,elcnt
!         WRITE(123,*) el_list(i)
!       ENDDO
!       CLOSE(123)
!       
!       OPEN(unit=123,file="node.list")
!       WRITE(123,*) ndcnt
!       DO i = 1,ndcnt
!         WRITE(123,*) nd_list(i)
!       ENDDO
!       CLOSE(123)      
              

      END SUBROUTINE plot_filled_contours_adapt
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE find_diff_nodes(fig,sol1,sol2)
            
      USE plot_globals, ONLY: xyplt   
      USE find_element, ONLY: in_element,find_element_init,find_element_final      
            
      IMPLICIT NONE      
      
      TYPE(plot_type), INTENT(IN) :: fig
      TYPE(solution_type), INTENT(IN) :: sol1
      TYPE(solution_type), INTENT(IN) :: sol2   
      
      INTEGER :: el,pt
      INTEGER :: et,i,nnd
      INTEGER :: elin
      REAL(rp) :: xy(2)
      REAL(rp) :: rs(2)    
      
      LOGICAL :: file_exists1,file_exists2          
      
      CALL find_element_init(sol2%nel_type,sol2%nverts,sol2%np,sol2%nnds,sol2%nn,sol2%xy,sol2%nepn,sol2%epn)
      
      ALLOCATE(elpt_diff(mnpp,sol1%ne),rspt_diff(2,mnpp,sol1%ne))
      
 elem:DO el = 1,sol1%ne
      
        IF (sol1%el_in(el) == 0) THEN
          CYCLE elem
        ENDIF
      
        et = sol1%el_type(el)
        i = fig%el_plt(el)    
        
        DO pt = 1,sol1%npplt(i)          
        
          xy(1) = xyplt(pt,el,1)
          xy(2) = xyplt(pt,el,2)
          CALL in_element(xy,sol2%el_type,sol2%elxy,elin,rs) 
          
          elpt_diff(pt,el) = elin
          rspt_diff(1,pt,el) = rs(1)
          rspt_diff(2,pt,el) = rs(2)

        ENDDO                       
        
      ENDDO elem
      
      CALL find_element_final()      
      
      RETURN
      END SUBROUTINE find_diff_nodes      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   

      SUBROUTINE contour_fill_element(file_unit,nptri,rect,sol_min,sol_max,sol_val,xplt,yplt)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: file_unit
      INTEGER, INTENT(IN) :: nptri
      INTEGER, DIMENSION(:,:), INTENT(IN) :: rect
      REAL(rp), INTENT(IN) :: sol_min
      REAL(rp), INTENT(IN) :: sol_max
      REAL(rp), DIMENSION(:), INTENT(IN) :: sol_val
      REAL(rp), DIMENSION(:), INTENT(IN) :: xplt
      REAL(rp), DIMENSION(:), INTENT(IN) :: yplt
      
      INTEGER :: tri,v,i
      INTEGER :: nd,lev
      REAL(rp) :: dc
      REAL(rp) :: sol_lev
      REAL(rp) :: color_val(3)        
      
      dc = (sol_max-sol_min)/real(ncolors-1,rp)             
      
      
      DO tri = 1,nptri

        DO v = 1,3  
   
          nd = rect(v,tri)
        
          sol_lev = sol_min
          lev = ncolors
  levels: DO i = 1,ncolors-1
            IF ((sol_val(nd) >= sol_lev) .and. (sol_val(nd) < sol_lev+dc)) THEN
              lev = i
              CALL interp_colors(lev,sol_lev,dc,colors,sol_val(nd),color_val)
              EXIT levels
            ENDIF
            sol_lev = sol_lev + dc
          ENDDO levels
          
          IF (sol_val(nd) <= sol_min) THEN
            lev = 1
            color_val(1) = colors(lev,1)
            color_val(2) = colors(lev,2)
            color_val(3) = colors(lev,3)
          ELSE IF (sol_val(nd) > sol_max) THEN
            lev = ncolors
            color_val(1) = colors(lev,1)
            color_val(2) = colors(lev,2)
            color_val(3) = colors(lev,3)           
          ENDIF          
          
          WRITE(file_unit,"(A,5(F9.5,1x),A)") "[",ax*xplt(nd)+bx,ay*yplt(nd)+by,color_val(1),color_val(2),color_val(3),"]"  
        ENDDO        

        WRITE(file_unit,"(A)") "trifill"        
        
      ENDDO       
      
      RETURN
      END SUBROUTINE contour_fill_element

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE plot_line_contours(file_unit,sol,nptri,rect,xyplt,rre,sre,snap,fig)

      USE evaluate_mod, ONLY: contour_line_newton,contour_line_linear
      USE transformation, ONLY: element_transformation
      USE shape_functions_mod, ONLY: shape_functions_area_eval, shape_functions_edge_eval
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: file_unit
      TYPE(solution_type), INTENT(IN) :: sol
      INTEGER, DIMENSION(:), INTENT(IN) :: nptri
      INTEGER, DIMENSION(:,:,:), INTENT(IN) :: rect
      REAL(rp), DIMENSION(:,:,:), INTENT(IN) :: xyplt  
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: rre
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: sre
      INTEGER, INTENT(IN) :: snap
      TYPE(plot_type), INTENT(IN) :: fig

      
      INTEGER :: i,j,vrt,it
      INTEGER :: el,nd,dof,lev,tri,ig
      INTEGER :: et,nv,nnd,ord    
      REAL(rp) :: sol_lev
      REAL(rp) :: dc
      INTEGER :: above(3)
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: el_below
      INTEGER :: n1,n2
      INTEGER :: nd1,nd2
      REAL(rp) :: r1,r2,s1,s2,f1,f2
      REAL(rp) :: x(2),y(2)
      REAL(rp) :: re(1),se(1)
      INTEGER :: fail_flag
      INTEGER :: sol_snap
      REAL(rp) :: l(sol%mnnds,1)
      REAL(rp) :: xv(3),yv(3)
      REAL(rp) :: rv(3),sv(3),fv(3)
      REAL(rp) :: hbe(2),Ze(2),ve(2),fe(2)
      REAL(rp) :: e,a,b,c,w,u,v
      
      

      ! Account for initial condition in dgswe output      
      IF (sol%output_type == "dgswe") THEN
        sol_snap = snap + 1
      ELSE 
        sol_snap = snap
      ENDIF      
      
      

      ALLOCATE(el_below(maxval(nptri),sol%ne))
      el_below = 0
      
       dc = (fig%sol_max-fig%sol_min)/real(nctick-1,rp)                  
       sol_lev = fig%sol_min      
       
levels:DO lev = 1,nctick
      
    elem:DO el = 1,sol%ne
           et = sol%el_type(el)
           IF (sol%el_in(el) == 0) THEN
             CYCLE elem
           ENDIF
     
           ord = fig%el_plt(el)
           
    subtri:DO tri = 1,nptri(ord)
    
             IF (el_below(tri,el) == 1) THEN
               CYCLE subtri
             ENDIF 

!              PRINT("(3(A,I9))"), "lev = ", lev, " el = ", el, " tri = ",tri
             
             DO vrt = 1,3  
   
               nd = rect(vrt,tri,ord)
          
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
               nd = rect(i,tri,ord)
               rv(i) = rre(nd,ord)
               sv(i) = sre(nd,ord)
               fv(i) = fig%sol_val(nd,el)
             ENDDO

             
             i = 0             
             
       edge: DO vrt = 1,3
               n1 = mod(vrt+0,3)+1
               n2 = mod(vrt+1,3)+1
                              
               IF ((above(n1) == 1 .and. above(n2) == 1) .or. &
                   (above(n1) == 0 .and. above(n2) == 0)) THEN    
                   CYCLE edge
               ENDIF
                   
               r1 = rv(n1)
               r2 = rv(n2)
               s1 = sv(n1)
               s2 = sv(n2)   
               f1 = fv(n1)
               f2 = fv(n2)
               
!                CALL contour_line_newton(sol,fig%type_flag,sol_lev,el,et,sol_snap,r1,r2,s1,s2,re,se)
               
               CALL contour_line_linear(f1,f2,sol_lev,r1,r2,s1,s2,re,se)
                 
               i = i + 1                 
 
               CALL shape_functions_area_eval(et,sol%np(et),nnd,1,re,se,l)
               CALL element_transformation(nnd,sol%elxy(:,el,1),sol%elxy(:,el,2),l(:,1),x(i),y(i))  
                 
 
               !!!! Edge midpoint !!!!
!                  i = i + 1
!                  nd1 = rect(n1,tri,et)
!                  nd2 = rect(n2,tri,et)                 
!                  x(i) = 0.5d0*(xyplt(nd1,el,1) + xyplt(nd2,el,1))
!                  y(i) = 0.5d0*(xyplt(nd1,el,2) + xyplt(nd2,el,2))                   
                 
             ENDDO edge

             IF (i > 2) THEN
               PRINT*, "Error: "
               STOP
             ENDIF 
             
             WRITE(file_unit,"(2(F9.5,1x))") ax*x(1)+bx,ay*y(1)+by
             WRITE(file_unit,"(2(F9.5,1x))") ax*x(2)+bx,ay*y(2)+by     
             WRITE(file_unit,"(A)") "draw-line"  
             
           ENDDO subtri
         ENDDO elem
         sol_lev = sol_lev + dc
      ENDDO levels

      END SUBROUTINE plot_line_contours

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE plot_mesh(file_unit,ne,nverts,nnds,pplt,el_type,el_in,psic,elxy,line_color)
      
      USE transformation, ONLY: element_transformation
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: file_unit
      INTEGER, INTENT(IN) :: ne
      INTEGER, DIMENSION(:), INTENT(IN) :: nverts
      INTEGER, DIMENSION(:), INTENT(IN) :: nnds
      INTEGER, DIMENSION(:), INTENT(IN) :: pplt
      INTEGER, DIMENSION(:), INTENT(IN) :: el_type
      INTEGER, DIMENSION(:), INTENT(IN) :: el_in      
      REAL(rp), DIMENSION(:,:,:), INTENT(IN) :: psic
      REAL(rp), DIMENSION(:,:,:), INTENT(IN) :: elxy
      REAL(rp), DIMENSION(:), INTENT(IN) :: line_color

      INTEGER :: el,nd
      INTEGER :: i,et
      INTEGER :: nnd,npts
      REAL(rp), DIMENSION(:), ALLOCATABLE :: xpt,ypt            
      
      ALLOCATE(xpt(mnpp),ypt(mnpp))

      
 elem:DO el = 1,ne              
      
        IF (el_in(el) == 0) THEN
          CYCLE elem
        ENDIF
        
        et = el_type(el)
        
        IF (et == 1) THEN         
          WRITE(file_unit,"(3(F9.5,1x))")  line_color(1), line_color(2), line_color(3)      
          DO nd = 1,nverts(et)       
            WRITE(file_unit,"(2(F9.5,1x))") ax*elxy(nd,el,1)+bx,ay*elxy(nd,el,2)+by    
          ENDDO                 
          WRITE(file_unit,"(A)") "draw-tri-element-color"              
        ELSE IF (et == 2) THEN
          WRITE(file_unit,"(3(F9.5,1x))")  line_color(1), line_color(2), line_color(3)          
          DO nd = 1,nverts(et)
            WRITE(file_unit,"(2(F9.5,1x))") ax*elxy(nd,el,1)+bx,ay*elxy(nd,el,2)+by    
          ENDDO                
          WRITE(file_unit,"(A)") "draw-quad-element-color"           
        ELSE
          npts = nverts(et)*pplt(et)
          nnd = nnds(et)
          DO nd = 1,npts
            CALL element_transformation(nnd,elxy(:,el,1),elxy(:,el,2),psic(:,nd,et),xpt(nd),ypt(nd))
          ENDDO
          WRITE(file_unit,"(A)") "newpath"
          WRITE(file_unit,"(2(F9.5,1x),A)") ax*xpt(1)+bx,ay*ypt(1)+by,"moveto" 
          DO nd = 2,npts
            WRITE(file_unit,"(2(F9.5,1x),A)") ax*xpt(nd)+bx,ay*ypt(nd)+by, "lineto"
          ENDDO
          WRITE(file_unit,"(A)") "closepath"
          WRITE(file_unit,"(A)") ".5 setlinewidth 2 setlinejoin"          
          WRITE(file_unit,"(3(F9.5,1x))")  line_color(1), line_color(2), line_color(3)               
          WRITE(file_unit,"(A)") "gsave setrgbcolor stroke grestore"
        ENDIF
      ENDDO elem

      END SUBROUTINE plot_mesh
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE plot_vis_mesh(file_unit,ne,el_type,el_in,xyplt,nptri,rect,fig)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: file_unit
      INTEGER, INTENT(IN) :: ne
      INTEGER, DIMENSION(:), INTENT(IN) :: el_type
      INTEGER, DIMENSION(:), INTENT(IN) :: el_in      
      REAL(rp), DIMENSION(:,:,:), INTENT(IN) :: xyplt
      INTEGER, DIMENSION(:), INTENT(IN) :: nptri
      INTEGER, DIMENSION(:,:,:), INTENT(IN) :: rect
      TYPE(plot_type), INTENT(IN) :: fig      

      INTEGER :: el,nd,tri,v
      INTEGER :: et,ord      

      
 elem:DO el = 1,ne
        
        IF (el_in(el) == 0) THEN
          CYCLE elem
        ENDIF
        
        et = el_type(el)
        ord = fig%el_plt(el)        

        DO tri = 1,nptri(ord)        
          DO v = 1,3        
            nd = rect(v,tri,ord)
            WRITE(file_unit,"(2(F9.5,1x))") ax*xyplt(nd,el,1)+bx,ay*xyplt(nd,el,2)+by
          ENDDO                
          WRITE(file_unit,"(A)") "draw-tri-element"            
        ENDDO
      ENDDO elem

      END SUBROUTINE plot_vis_mesh  
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE plot_decomp_mesh(file_unit,nverts)
      
      USE grid_file_mod, ONLY: read_header,read_coords,read_connectivity
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: file_unit
      INTEGER, DIMENSION(:), INTENT(IN) :: nverts 

      INTEGER :: el,nd,pe,i
      INTEGER :: et,nv
      
      INTEGER :: npe
      CHARACTER(6) :: dirname
      INTEGER, PARAMETER :: lname = 6
      CHARACTER(100) :: grd_name
      INTEGER :: ne_pe,nn_pe
      REAL(rp), DIMENSION(:,:), ALLOCATABLE :: xy_pe
      REAL(rp), DIMENSION(:), ALLOCATABLE :: depth_pe
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: ect_pe
      INTEGER, DIMENSION(:), ALLOCATABLE :: el_type_pe
      INTEGER :: color
      REAL(rp) :: xpt,ypt
      INTEGER :: outside
          
      
      dirname = "PE0000"
      
      OPEN(unit=8080,file='PE0000/fort.80')
      READ(8080,*) npe
      CLOSE(8080)
      
      i = 0      
      
  pes:DO pe = 1,npe
  
        WRITE(dirname(3:lname),"(I4.4)") pe-1     
        CALL read_header(0,dirname(1:lname)//'/'//'fort.14',grd_name,ne_pe,nn_pe)
        CALL read_coords(nn_pe,xy_pe,depth_pe)
        CALL read_connectivity(ne_pe,ect_pe,el_type_pe)        
        CLOSE(14)
      
        i = i + 1
        color = mod(i,ncolors)+1      

   elem:DO el = 1,ne_pe        
        
          et = el_type_pe(el)
          nv = nverts(et)
          outside = 0
          DO nd = 1,nv
            xpt = xy_pe(1,ect_pe(nd,el))
            ypt = xy_pe(2,ect_pe(nd,el))
            
            IF (xpt < xbox_min .or. xpt > xbox_max) THEN
              outside = 1
            ENDIF
           
            IF (ypt < ybox_min .or. ypt > ybox_max) THEN
              outside = 1
            ENDIF            
          ENDDO
          
          IF (outside == 1) THEN
            CYCLE elem
          ENDIF
        

        
          IF (et == 1) THEN
            WRITE(file_unit,"(3(F9.5,1x))") colors(color,1), colors(color,2), colors(color,3)          
            DO nd = 1,nverts(et)
              WRITE(file_unit,"(2(F9.5,1x))") ax*xy_pe(1,ect_pe(nd,el))+bx,ay*xy_pe(2,ect_pe(nd,el))+by    
            ENDDO                
            WRITE(file_unit,"(A)") "draw-tri-element-color"         
          ELSE IF (et == 2) THEN
            WRITE(file_unit,"(3(F9.5,1x))") colors(color,1), colors(color,2), colors(color,3)             
            DO nd = 1,nverts(et)         
              WRITE(file_unit,"(2(F9.5,1x))") ax*xy_pe(1,ect_pe(nd,el))+bx,ay*xy_pe(2,ect_pe(nd,el))+by    
            ENDDO                
            WRITE(file_unit,"(A)") "draw-quad-element-color"           
          ENDIF
        ENDDO elem
      
      ENDDO pes

      END SUBROUTINE plot_decomp_mesh      

      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   

     SUBROUTINE plot_boundaries(file_unit,nverts,nnds,pplt,nbed,bedn,ged2el,ged2led,el_type,el_in,psic,elxy)
     
     USE transformation, ONLY: element_transformation

     IMPLICIT NONE
     
     INTEGER, INTENT(IN) :: file_unit
     INTEGER, DIMENSION(:), INTENT(IN) :: nverts
     INTEGER, DIMENSION(:), INTENT(IN) :: nnds
     INTEGER, DIMENSION(:), INTENT(IN) :: pplt
     INTEGER, INTENT(IN) :: nbed
     INTEGER, DIMENSION(:), INTENT(IN) :: bedn
     INTEGER, DIMENSION(:,:), INTENT(IN) :: ged2el
     INTEGER, DIMENSION(:,:), INTENT(IN) :: ged2led
     INTEGER, DIMENSION(:), INTENT(IN) :: el_type
     INTEGER, DIMENSION(:), INTENT(IN) :: el_in
     REAL(rp), DIMENSION(:,:,:), INTENT(IN) :: psic
     REAL(rp), DIMENSION(:,:,:), INTENT(IN) :: elxy
     
     INTEGER :: ed,nd,j
     INTEGER :: el,et,nv
     INTEGER :: ged,led
     INTEGER :: n1,n2
     INTEGER :: nnd,npts
     REAL(rp), DIMENSION(:), ALLOCATABLE :: xpt,ypt            
      
     ALLOCATE(xpt(mnpp),ypt(mnpp))
     
     
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
          
          WRITE(file_unit,"(2(F9.5,1x))") ax*elxy(n1,el,1)+bx,ay*elxy(n1,el,2)+by
          WRITE(file_unit,"(2(F9.5,1x))") ax*elxy(n2,el,1)+bx,ay*elxy(n2,el,2)+by
          WRITE(file_unit,"(A)") "draw-line"          
        ELSE

          npts = nverts(et)*pplt(et)
          nnd = nnds(et)
          DO nd = 1,npts
            CALL element_transformation(nnd,elxy(:,el,1),elxy(:,el,2),psic(:,nd,et),xpt(nd),ypt(nd))
          ENDDO
          
          n1 = mod(led,nv)*pplt(et) + 1
          n2 = n1 + pplt(et)
          WRITE(file_unit,"(A)") "newpath"
          WRITE(file_unit,"(2(F9.5,1x),A)") ax*xpt(n1)+bx,ay*ypt(n1)+by,"moveto" 
          DO j = 1, pplt(et)
            nd = n1 + j
            IF (nd == nv*pplt(et)+1) THEN
              nd = 1
            ENDIF
            WRITE(file_unit,"(2(F9.5,1x),A)") ax*xpt(nd)+bx,ay*ypt(nd)+by, "lineto"
          ENDDO
          WRITE(file_unit,"(A)") ".5 setlinewidth 2 setlinejoin"          
          WRITE(file_unit,"(A)") "stroke"              
        ENDIF
       
     ENDDO edge
     
     END SUBROUTINE plot_boundaries
     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE plot_nodestring(file_unit,nbou,nvel,fbseg,fbnds,xy)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: file_unit
      INTEGER, INTENT(IN) :: nbou
      INTEGER, INTENT(IN) :: nvel
      INTEGER, DIMENSION(:,:), INTENT(IN) :: fbseg
      INTEGER, DIMENSION(:,:), INTENT(IN) :: fbnds
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: xy
      
      INTEGER :: bou,j
      INTEGER :: nd
      INTEGER :: segtype,nbnds
      INTEGER :: started
      INTEGER, DIMENSION(:), ALLOCATABLE :: ndflag
      
      ALLOCATE(ndflag(nvel))
      
      DO bou = 1,nbou
        segtype = fbseg(2,bou)
        
        IF (segtype /= 30) THEN
          CYCLE
        ENDIF
                
        nbnds = fbseg(1,bou)        
        
        ndflag = 0
        DO j = 1,nbnds
        
          nd = fbnds(j,bou)        
            
          IF (xy(1,nd) < xbox_min .or. xy(1,nd) > xbox_max) THEN
            ndflag(j) = 1
          ENDIF
           
          IF (xy(2,nd) < ybox_min .or. xy(2,nd) > ybox_max) THEN
            ndflag(j) = 1
          ENDIF     
        ENDDO
        
        
        started = 0
        DO j = 1,nbnds
          nd = fbnds(j,bou)
          IF (started == 0 .and. ndflag(j) == 0) THEN
            WRITE(file_unit,"(A)") "newpath"                  
            WRITE(file_unit,"(2(F9.5,1x),A)") ax*xy(1,nd)+bx,ay*xy(2,nd)+by,"moveto"   
            started = 1
          ELSE IF (started == 1 .and. ndflag(j) == 0) THEN
            WRITE(file_unit,"(2(F9.5,1x),A)") ax*xy(1,nd)+bx,ay*xy(2,nd)+by, "lineto"          
          ENDIF
        ENDDO
        
        IF (started /= 0) THEN
          WRITE(file_unit,"(A)") "2 setlinewidth 2 setlinejoin"    
          WRITE(file_unit,"(A)") "gsave 1 .5 0 setrgbcolor stroke grestore"        
        ENDIF
      ENDDO
      
      RETURN
      END SUBROUTINE plot_nodestring


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE plot_cb_nodes(file_unit,ctp,nbou,fbseg,fbnds,xy,bndxy,nepn,epn,el_in)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: file_unit
      INTEGER, INTENT(IN) :: ctp
      INTEGER, INTENT(IN) :: nbou
      INTEGER, DIMENSION(:,:), INTENT(IN) :: fbseg
      INTEGER, DIMENSION(:,:), INTENT(IN) :: fbnds
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: xy
      REAL(rp), DIMENSION(:,:,:,:), INTENT(IN) :: bndxy
      INTEGER, DIMENSION(:), INTENT(IN) :: nepn
      INTEGER, DIMENSION(:,:), INTENT(IN) :: epn
      INTEGER, DIMENSION(:), INTENT(IN) :: el_in
      
      INTEGER :: i,j,k
      INTEGER :: nd,el
      INTEGER :: nbnds
      INTEGER :: btype
      INTEGER :: in_flag
      
      DO i = 1,nbou
        nbnds = fbseg(1,i)
        btype = fbseg(2,i)
        
        IF( btype == 0 .OR. btype == 10 .OR. btype == 20  .OR. &   ! land boundaries
            btype == 1 .OR. btype == 11 .OR. btype == 21 ) THEN    ! island boundaries
            
          DO j = 1,nbnds-1
            
            
            nd = fbnds(j,i) 
            
            in_flag = 0
            DO k = 1,nepn(nd)
              el = epn(k,nd)
              IF (el_in(el) == 1) THEN
                in_flag = 1
              ENDIF
            ENDDO
            
            IF (in_flag /= 1) THEN
              CYCLE
            ENDIF
            
            
            WRITE(file_unit,"(3(I5,1x))") 1,0,0                    
            WRITE(file_unit,"(2(F9.5,1x))") ax*xy(1,nd)+bx, ay*xy(2,nd)+by
            WRITE(file_unit,"(I5)") 2                 
            WRITE(file_unit,"(A)") "draw-dot"
            
            DO k = 1,ctp-1
              WRITE(file_unit,"(3(I5,1x))") 1,0,0             
              WRITE(file_unit,"(2(F9.5,1x))") ax*bndxy(1,k,j,i)+bx, ay*bndxy(2,k,j,i)+by
              WRITE(file_unit,"(I5)") 2                   
              WRITE(file_unit,"(A)") "draw-dot"              
            ENDDO
          
          ENDDO
          
          nd = fbnds(nbnds,i) 
            
          in_flag = 0
          DO k = 1,nepn(nd)
            el = epn(k,nd)
            IF (el_in(el) == 1) THEN
              in_flag = 1
            ENDIF
          ENDDO
          
          IF (in_flag == 1) THEN
            WRITE(file_unit,"(3(I5,1x))") 1,0,0           
            WRITE(file_unit,"(2(F9.5,1x))") ax*xy(1,nbnds)+bx, ay*xy(2,nbnds)+by
            WRITE(file_unit,"(I5)") 2                 
            WRITE(file_unit,"(A)") "draw-dot"         
          ENDIF  
        ENDIF    
      ENDDO
      
      RETURN
      END SUBROUTINE plot_cb_nodes
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE plot_elxy_nodes(file_unit,ne,el_type,el_in,nnds,elxy)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: file_unit    
      INTEGER, INTENT(IN) :: ne
      INTEGER, DIMENSION(:), INTENT(IN) :: el_type
      INTEGER, DIMENSION(:), INTENT(IN) :: el_in
      INTEGER, DIMENSION(:), INTENT(IN) :: nnds
      REAL(rp), DIMENSION(:,:,:), INTENT(IN) :: elxy
      
      INTEGER :: el,nd
      INTEGER :: et
      INTEGER :: nnd

      DO el = 1,ne
        et = el_type(el)
        nnd = nnds(et)
        
        IF (el_in(el) == 0) THEN
          CYCLE
        ENDIF
        
        DO nd = 1,nnd
          WRITE(file_unit,"(3(I5,1x))") 0,0,0         
          WRITE(file_unit,"(2(F9.5,1x))") ax*elxy(nd,el,1)+bx, ay*elxy(nd,el,2)+by
          WRITE(file_unit,"(I5)") 2               
          WRITE(file_unit,"(A)") "draw-dot"         
        ENDDO
      ENDDO

      
      RETURN
      END SUBROUTINE plot_elxy_nodes     
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE plot_qpts(file_unit,ne,p,ctp,np,el_type,el_in,elxy,nverts,ned,ed_type,ged2el,ged2led)
      
      USE transformation, ONLY: element_transformation               
      USE area_qpts_mod, ONLY: area_qpts
      USE edge_qpts_mod, ONLY: edge_qpts
      USE shape_functions_mod, ONLY: shape_functions_area_eval 
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: file_unit
      INTEGER, INTENT(IN) :: ne
      INTEGER, INTENT(IN) :: p
      INTEGER, INTENT(IN) :: ctp
      INTEGER, DIMENSION(:), INTENT(IN) :: np
      INTEGER, DIMENSION(:), INTENT(IN) :: el_type
      INTEGER, DIMENSION(:), INTENT(IN) :: el_in
      REAL(rp), DIMENSION(:,:,:), INTENT(IN) :: elxy
      INTEGER, DIMENSION(:), INTENT(IN) :: nverts
      INTEGER, INTENT(IN) :: ned
      INTEGER, DIMENSION(:), INTENT(IN) :: ed_type
      INTEGER, DIMENSION(:,:), INTENT(IN) :: ged2el
      INTEGER, DIMENSION(:,:), INTENT(IN) :: ged2led
      
      INTEGER :: el,pt,nd,i,ged
      INTEGER :: et,nnd,edt,led
      INTEGER :: mnqpta,nqpta(nel_type)   
      INTEGER :: mnqpte,nqpte(nel_type)    
      INTEGER :: nnds(nel_type)
      INTEGER :: mnnds,mnp
      REAL(rp) :: xpt,ypt
      REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: qpta
      REAL(rp), DIMENSION(:,:), ALLOCATABLE ::  wpta 
      REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: qpte
      REAL(rp), DIMENSION(:,:), ALLOCATABLE ::  wpte         
      REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: psi     
      
      CALL area_qpts(1,p,ctp,nel_type,nqpta,mnqpta,wpta,qpta)           
      
      mnp = maxval(np)      
      mnnds = (mnp+1)**2           
      ALLOCATE(psi(mnnds,mnqpta,nel_type))      
      DO et = 1,nel_type         
        CALL shape_functions_area_eval(et,np(et),nnds(et),nqpta(et),qpta(:,1,et),qpta(:,2,et), &
                                       psi(:,:,et))
      ENDDO            
      
      DO el = 1,ne
        et = el_type(el)
        nnd = nnds(et)
        
        IF (el_in(el) == 0) THEN
          CYCLE
        ENDIF
        
        DO pt = 1,nqpta(et)             
          CALL element_transformation(nnd,elxy(:,el,1),elxy(:,el,2),psi(:,pt,et),xpt,ypt)           
          WRITE(file_unit,"(3(I5,1x))") 0,0,0         
          WRITE(file_unit,"(2(F9.5,1x))") ax*xpt+bx, ay*ypt+by
          WRITE(file_unit,"(I5)") 2               
          WRITE(file_unit,"(A)") "draw-dot"         
        ENDDO                
      ENDDO  
      DEALLOCATE(psi)
      
      
      
      
      CALL edge_qpts(1,p,ctp,nel_type,nqpte,mnqpte,wpte,qpte)           
      
      mnp = maxval(np)      
      mnnds = (mnp+1)**2           
      ALLOCATE(psi(mnnds,4*mnqpte,nel_type))      
      
      
      DO ged = 1,ned
        edt = ed_type(ged)
        el = ged2el(1,ged)
        led = ged2led(1,ged)
        et = el_type(el)
        
        nnd = nnds(et)          
       
        IF (el_in(el) == 0) THEN
          CYCLE
        ENDIF       
       
        IF (edt /= 10 .and. et > 2) THEN
          edt = 1
        ELSE 
          edt = et
        ENDIF                 
               
        CALL shape_functions_area_eval(et,np(et),nnds(et),nverts(et)*nqpte(edt),qpte(:,1,edt),qpte(:,2,edt), &
                                       psi(:,:,et))          
        
        DO i = 1,nqpte(edt)   
          pt = (led-1)*nqpte(edt)+i
          CALL element_transformation(nnd,elxy(:,el,1),elxy(:,el,2),psi(:,pt,et),xpt,ypt)           
          WRITE(file_unit,"(3(I5,1x))") 0,0,0         
          WRITE(file_unit,"(2(F9.5,1x))") ax*xpt+bx, ay*ypt+by
          WRITE(file_unit,"(I5)") 2               
          WRITE(file_unit,"(A)") "draw-dot"         
        ENDDO          
        
      ENDDO
      
      
      RETURN
      END SUBROUTINE plot_qpts

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE fill_elements(file_unit,ne,nverts,nnds,pplt,el_type,el_in,psic,elxy)

      USE transformation, ONLY: element_transformation
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: file_unit
      INTEGER, INTENT(IN) :: ne
      INTEGER, DIMENSION(:), INTENT(IN) :: nverts
      INTEGER, DIMENSION(:), INTENT(IN) :: nnds
      INTEGER, DIMENSION(:), INTENT(IN) :: pplt
      INTEGER, DIMENSION(:), INTENT(IN) :: el_type
      INTEGER, DIMENSION(:), INTENT(IN) :: el_in   
      REAL(rp), DIMENSION(:,:,:), INTENT(IN) :: psic
      REAL(rp), DIMENSION(:,:,:), INTENT(IN) :: elxy

      INTEGER :: el,nd
      INTEGER :: npts
      INTEGER :: nnd
      INTEGER :: et
      INTEGER :: i,n
      INTEGER :: fill
      INTEGER, ALLOCATABLE, DIMENSION(:) :: fill_list
      LOGICAL :: file_exists
      REAL(rp), DIMENSION(:), ALLOCATABLE :: xpt,ypt            
      
      ALLOCATE(xpt(mnpp),ypt(mnpp))
      
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

!       fill_list = 1
      
 elem:DO el = 1,ne
        
        IF (el_in(el) == 0) THEN
          CYCLE elem
        ENDIF
        
        et = el_type(el)
        fill = fill_list(el)        
        
        IF (et == 1) THEN
          IF (fill == 1) THEN
            WRITE(file_unit,"(3(I5,1x))") 0,0,1          
          ELSE
            WRITE(file_unit,"(3(I5,1x))") 1,1,1  
          ENDIF
          
          DO nd = 1,nverts(et)
            WRITE(file_unit,"(2(F9.5,1x))") ax*elxy(nd,el,1)+bx,ay*elxy(nd,el,2)+by
          ENDDO                
          WRITE(file_unit,"(A)") "fill-tri-element"            
        ELSE IF (et == 2) THEN        
          IF (fill == 1) THEN
            WRITE(file_unit,"(3(I5,1x))") 0,0,1          
          ELSE
            WRITE(file_unit,"(3(I5,1x))") 1,1,1  
          ENDIF        

          DO nd = 1,nverts(et)
            WRITE(file_unit,"(2(F9.5,1x))") ax*elxy(nd,el,1)+bx,ay*elxy(nd,el,2)+by
          ENDDO                
          WRITE(file_unit,"(A)") "fill-quad-element"            
       ELSE   
          npts = nverts(et)*pplt(et)
          nnd = nnds(et)
          DO nd = 1,npts
            CALL element_transformation(nnd,elxy(:,el,1),elxy(:,el,2),psic(:,nd,et),xpt(nd),ypt(nd))
          ENDDO
          WRITE(file_unit,"(A)") "newpath"
          WRITE(file_unit,"(2(F9.5,1x),A)") ax*xpt(1)+bx,ay*ypt(1)+by,"moveto" 
          DO nd = 2,npts
            WRITE(file_unit,"(2(F9.5,1x),A)") ax*xpt(nd)+bx,ay*ypt(nd)+by, "lineto"
          ENDDO
          
          IF (fill == 1) THEN
            WRITE(file_unit,"(A)") "gsave 0 0 1 setrgbcolor fill grestore"
          ELSE
            WRITE(file_unit,"(A)") "gsave 1 1 1 setrgbcolor fill grestore"      
          ENDIF
        ENDIF
      ENDDO elem

      END SUBROUTINE fill_elements     

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   

      SUBROUTINE plot_stations(file_unit,sta_start,sta_end,xysta,sta_pick)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: file_unit
      INTEGER, INTENT(IN) :: sta_start
      INTEGER, INTENT(IN) :: sta_end
      REAL(rp), DIMENSION(:,:), INTENT(IN) :: xysta
      INTEGER, INTENT(IN), OPTIONAL :: sta_pick
      
      INTEGER :: sta
      REAL(rp) :: x,y
     
      DO sta = sta_start,sta_end
            
        x = xysta(1,sta)
        y = xysta(2,sta)
            
        IF (x > xbox_min .AND. x < xbox_max .AND. y > ybox_min .AND. y < ybox_max) THEN  
          WRITE(file_unit,"(3(F9.5,1x))") .89,0.1,0.11        
          WRITE(file_unit,"(2(F9.5,1x))") ax*x+bx,ay*y+by
          WRITE(file_unit,"(I5)") 2              
          WRITE(file_unit,"(A)") "draw-dot"        
        ENDIF
      
      ENDDO         

      
      IF (PRESENT(sta_pick)) THEN
        x = xysta(1,sta_pick)
        y = xysta(2,sta_pick)
            
        IF (x > xbox_min .AND. x < xbox_max .AND. y > ybox_min .AND. y < ybox_max) THEN  
          WRITE(file_unit,"(3(F9.5,1x))") 0.21,0.49,.72
          WRITE(file_unit,"(2(F9.5,1x))") ax*x+bx,ay*y+by
          WRITE(file_unit,"(I5)") 6            
          WRITE(file_unit,"(A)") "draw-dot"        
        ENDIF      
      ENDIF
      
      END SUBROUTINE plot_stations

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE plot_region_box(file_unit,region_box)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: file_unit
      REAL(rp), INTENT(IN) :: region_box(4)
      
      
      WRITE(file_unit,"(2(F9.5,1x))") ax*region_box(1)+bx,ay*region_box(3)+by        
      WRITE(file_unit,"(2(F9.5,1x))") ax*region_box(2)+bx,ay*region_box(3)+by                    
      WRITE(file_unit,"(2(F9.5,1x))") ax*region_box(2)+bx,ay*region_box(4)+by        
      WRITE(file_unit,"(2(F9.5,1x))") ax*region_box(1)+bx,ay*region_box(4)+by              
      WRITE(file_unit,"(A)") "region-box"       
      
      END SUBROUTINE plot_region_box
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE plot_xy(file_unit,npt,xvec,yvec,rgb,width)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: file_unit
      INTEGER, INTENT(IN) :: npt
      REAL(rp), DIMENSION(:), INTENT(IN) :: xvec
      REAL(rp), DIMENSION(:), INTENT(IN) :: yvec
      REAL(rp), DIMENSION(:), INTENT(IN) :: rgb
      REAL(rp), INTENT(IN) :: width
      
      INTEGER :: pt
!       PRINT*, rgb(1),rgb(2),rgb(3)
      WRITE(file_unit,"(3(F9.5,1x),A)") rgb(1),rgb(2),rgb(3)," setrgbcolor"      
      WRITE(file_unit,"(A)") "newpath"        
      WRITE(file_unit,"(2(F9.5,1x),A)") ax*xvec(1)+bx,ay*yvec(1)+by, " moveto"      
      DO pt = 2,npt                
        WRITE(file_unit,"(2(F9.5,1x),A)") ax*xvec(pt)+bx,ay*yvec(pt)+by, " lineto"      
      ENDDO
      WRITE(file_unit,"(F9.5,1x,A)") width," setlinewidth 2 setlinejoin"
      WRITE(file_unit,"(A)") "stroke"            
      WRITE(file_unit,"(A)")  "0 0 0 setrgbcolor"       
      
      RETURN
      END SUBROUTINE plot_xy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE plot_xy_scatter(file_unit,npt,xvec,yvec,rgb,width)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: file_unit
      INTEGER, INTENT(IN) :: npt
      REAL(rp), DIMENSION(:), INTENT(IN) :: xvec
      REAL(rp), DIMENSION(:), INTENT(IN) :: yvec
      REAL(rp), DIMENSION(:), INTENT(IN) :: rgb
      REAL(rp), INTENT(IN) :: width
      
      INTEGER :: pt
     
      DO pt = 1,npt                
        WRITE(file_unit,"(3(F9.5,1x))") rgb(1),rgb(2),rgb(3)        
        WRITE(file_unit,"(2(F9.5,1x))") ax*xvec(pt)+bx,ay*yvec(pt)+by
        WRITE(file_unit,"(F9.5)") width              
        WRITE(file_unit,"(A)") "draw-dot"       
      ENDDO
    
      RETURN
      END SUBROUTINE plot_xy_scatter

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
        PRINT*, "colormap file does not exist: using default2 colormap"
        ncolors = 12
        ALLOCATE(colors(ncolors,3))        
        colors(1,:)  = (/ 0d0  , 0d0  , 139d0 /)
        colors(2,:)  = (/ 0d0  , 0d0  , 255d0 /)
        colors(3,:)  = (/ 125d0, 158d0, 192d0 /)
        colors(4,:)  = (/ 98d0 , 221d0, 221d0 /)
        colors(5,:)  = (/ 0d0  , 210d0, 0d0   /)
        colors(6,:)  = (/ 255d0, 255d0, 0d0   /)
        colors(7,:)  = (/ 255d0, 215d0, 0d0   /)
        colors(8,:)  = (/ 255d0, 104d0, 32d0  /)
        colors(9,:)  = (/ 251d0, 57d0 , 30d0  /)
        colors(10,:) = (/ 232d0, 0d0  , 0d0   /)
        colors(11,:) = (/ 179d0, 0d0  , 0d0   /)
        colors(12,:) = (/ 221d0, 0d0  , 221d0 /)        
      ELSE      
        OPEN(UNIT=101,FILE=TRIM(ADJUSTL(cmap_file)))
        READ(101,*) ncolors
        ALLOCATE(colors(ncolors,3))
        DO lev = 1,ncolors
          READ(101,*) (colors(lev,j), j=1,3)
        ENDDO
        CLOSE(101)
      ENDIF
      
      cmax = MAXVAL(colors)
      IF (cmax > 1d0+1d-12) THEN
        colors = colors/255d0
      ENDIF      
      
      RETURN
      END SUBROUTINE read_colormap      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      SUBROUTINE setup_cbounds(fig,sol,snap_start,snap_end)
      
      USE evaluate_mod, ONLY: evaluate_solution
      
      IMPLICIT NONE
      
      TYPE(plot_type), INTENT(INOUT) :: fig
      TYPE(solution_type), INTENT(IN) :: sol
      INTEGER, INTENT(IN) :: snap_start
      INTEGER, INTENT(IN) :: snap_end         
      
      INTEGER :: i,snap,el,pt
      INTEGER :: et,npts
      INTEGER :: start_snap,end_snap
      INTEGER :: sol_snap
      CHARACTER(:), ALLOCATABLE :: filename

           
      ! Early exit if figure is not being plotted, or if it is the mesh
      IF (fig%plot_sol_option == 0 .or. fig%type_flag == 1) THEN
        RETURN
      ENDIF
      
     
            
  
      
      ALLOCATE(fig%cscale_vals(snap_end+1,3)) 
      ALLOCATE(fig%sol_val(mnpp,sol%ne))       
      
      fig%num_cscale_vals = snap_end - snap_start + 1      
      
      ! Setup color scale bounds output file
      IF (fig%type_flag > 2) THEN   
        filename = TRIM(ADJUSTL(fig%name))//".cscale"          
        OPEN(unit=fig%cscale_unit,file=filename//".out")  
        WRITE(fig%cscale_unit,"(2I5)") snap_start-1,snap_end-1
      ENDIF      
      
      ! Read in minimum and maximum solution values from file
      IF (fig%cscale_option == "file") THEN
        OPEN(unit=22,file=filename)
        READ(22,*) start_snap,end_snap
        fig%num_cscale_vals = end_snap - start_snap + 1           
        DO snap = snap_start,snap_end
          READ(22,*) fig%cscale_vals(snap,1),fig%cscale_vals(snap,2),fig%cscale_vals(snap,3)
        ENDDO
        CLOSE(22)
      ENDIF
      
      ! Set specified minimum and maximum values for all timesnaps
      IF (fig%cscale_option == "spec") THEN
        DO snap = 1,snap_end+1
          fig%cscale_vals(snap,2) = fig%cscale_min
          fig%cscale_vals(snap,3) = fig%cscale_max
        ENDDO             
      ENDIF
      
      ! Find minimum and maximum solution values for each timesnap
      IF (fig%cscale_option == "auto-snap" .or. fig%cscale_option == "auto-all") THEN
      
        DO snap = snap_start,snap_end 
    
          fig%cscale_min = min_init
          fig%cscale_max = max_init  
                 
    
    elem: DO el = 1,sol%ne
    
            IF (sol%el_in(el) == 0) THEN
              CYCLE elem
            ENDIF
            
            et = sol%el_type(el)           
            npts = sol%npplt(et)       
            
            CALL evaluate_solution(el,fig%type_flag,snap,sol,fig%sol_val(:,el),npts,phi_sol=phi_sol(:,:,et))  
            
            DO pt = 1,npts
              IF (fig%sol_val(pt,el) > fig%cscale_max) THEN
                fig%cscale_max = fig%sol_val(pt,el)
              ENDIF
              
              IF (fig%sol_val(pt,el) < fig%cscale_min) THEN
                fig%cscale_min = fig%sol_val(pt,el)
              ENDIF
            ENDDO            
          
          ENDDO elem
                
          fig%cscale_vals(snap,2) = fig%cscale_min
          fig%cscale_vals(snap,3) = fig%cscale_max
        ENDDO            
      ENDIF
      
      ! Find overall minimum and maximum solution values
      IF (fig%cscale_option == "auto-all") THEN
        fig%cscale_min = min_init
        fig%cscale_max = max_init
        DO snap = snap_start,snap_end
          IF (fig%cscale_vals(snap,3) > fig%cscale_max) THEN
            fig%cscale_max = fig%cscale_vals(snap,3)
          ENDIF
          
          IF (fig%cscale_vals(snap,2) < fig%cscale_min) THEN
            fig%cscale_min = fig%cscale_vals(snap,2)
          ENDIF
        ENDDO
        
        DO snap = 1,snap_end
          fig%cscale_vals(snap,2) = fig%cscale_min
          fig%cscale_vals(snap,3) = fig%cscale_max          
        ENDDO
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

      SUBROUTINE plot_ref_el(filename,figure_width,et,ctp,ntri,ect,r,s)
      
      USE basis, ONLY: element_nodes
      USE axes_mod, ONLY: write_xyaxis    
      USE labels_mod, ONLY: write_xyaxis_labels,write_texheader,run_latex, & 
                            write_latex_ps_body,remove_latex_files      
      
      IMPLICIT NONE
      
      CHARACTER(*), INTENT(IN) :: filename
      REAL(rp), INTENT(IN) :: figure_width
      INTEGER, INTENT(IN) :: et
      INTEGER, INTENT(IN) :: ctp
      INTEGER, INTENT(IN) :: ntri
      INTEGER, DIMENSION(:,:), INTENT(IN) :: ect
      REAL(rp), DIMENSION(:), INTENT(IN) :: r
      REAL(rp), DIMENSION(:), INTENT(IN) :: s
      LOGICAL :: tex_file_exists
      
      INTEGER :: el,nd
      INTEGER :: file_unit
      INTEGER :: nnd
      REAL(rp) :: r_min,r_max,s_min,s_max
      REAL(rp) :: a_x,b_x,a_y,b_y
      REAL(rp), DIMENSION(:), ALLOCATABLE :: rc,sc
      
      file_unit = 90
      
      r_max = 1d0
      r_min = -1d0
      s_max = 1d0
      s_min = -1d0
      
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
            
      ax = ((rmax_axes-rmin_axes)/(r_max-r_min))
      bx = (rmin_axes*r_max-rmax_axes*r_min)/(r_max-r_min)
      
      ay = ax    
      by = smin_axes-ax*s_min     
      
      smax_axes = ay*s_max + by      

      CALL write_texheader()              
      CALL write_xyaxis_labels("rs")  
      CALL run_latex(tex_file_exists)         

      CALL write_psheader(TRIM(filename),file_unit)     
      
      DO el = 1,ntri
        DO nd = 1,3
          WRITE(file_unit,"(2(F9.5,1x))") ax*r(ect(nd,el))+bx,ay*s(ect(nd,el))+by                       
        ENDDO
        WRITE(file_unit,"(A)") "draw-tri-element"           
      ENDDO
      
      ALLOCATE(rc((ctp+1)**2),sc((ctp+1)**2))
      CALL element_nodes(et,1,ctp,nnd,rc,sc)          
    
      CALL write_xyaxis(file_unit,1)  
      
      DO nd = 1,nnd
        WRITE(file_unit,"(3(I5,1x))") 0,0,0       
        WRITE(file_unit,"(2(F9.5,1x))") ax*rc(nd)+bx,ay*sc(nd)+by
        WRITE(file_unit,"(I5)") 2             
        WRITE(file_unit,"(A)") "draw-dot"            
      ENDDO         
      
      CALL write_latex_ps_body(file_unit)      
    
      WRITE(file_unit,"(A)") "showpage"           
      CLOSE(file_unit)   
      
      CALL remove_latex_files()      

      END SUBROUTINE plot_ref_el
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      END MODULE plot_mod