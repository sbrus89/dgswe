      SUBROUTINE read_plot_input()

      USE plot_globals, ONLY: rp,input_path,cmap_file,ps,pc,p_low,p_high,p_skip, &
                              frmt,density, &
                              xbox_min,xbox_max,ybox_min,ybox_max,figure_width, &
                              snap_start,snap_end, &
                              zeta,bathy,vel,mesh, &
                              fontsize,nxtick,nytick,nctick, &
                              nxdec,nydec,ncdec,ntdec, &
                              substitute_path,replace_path,sub_path

      IMPLICIT NONE
      
      INTEGER, PARAMETER :: ninp = 27
      INTEGER :: i,n
      INTEGER :: inp_read,skipped
      INTEGER :: rm_ps
      CHARACTER(100) :: temp      
      CHARACTER(100) :: zbox  
      CHARACTER(100) :: ytick   
      CHARACTER(100) :: cscale 
      REAL(rp) :: tmp
      LOGICAL :: file_exists      

      
      zeta%cscale_unit = 20
      vel%cscale_unit = 21
      
      INQUIRE(file='plot.inp',exist=file_exists)
      IF (file_exists == .FALSE.) THEN
        PRINT*, "plot.inp file does not exist"
        STOP
      ENDIF      
      OPEN(unit=15,file='plot.inp')
      
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", "             Input Information               "
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", " "
      
      inp_read = 0
      skipped = 0
      DO WHILE (inp_read < ninp)
      
        READ(15,"(A100)") temp
                    
        IF ( INDEX(temp,"!") == 1 .or. INDEX(temp,"          ") == 1) THEN
!           PRINT*, "input ignored"
            skipped = skipped + 1
        ELSE
          inp_read = inp_read + 1
          SELECT CASE (inp_read)
            CASE (1)
              input_path = TRIM(ADJUSTL(temp))
              PRINT("(A,A)"), "input file path = ", input_path 
              
            CASE(2)
              IF (TRIM(ADJUSTL(temp)) == "0") THEN
                substitute_path = 0
                PRINT("(A,I3)"), "substitute path: ", substitute_path
              ELSE 
                substitute_path = 1
                n = INDEX(temp,",")
                replace_path = temp(1:n-1)
                sub_path = temp(n+1:)
                PRINT("(A,A,A,A)"), "substitute path: ", TRIM(ADJUSTL(replace_path)),", ",TRIM(ADJUSTL(sub_path))                
              ENDIF
              
            CASE (3)
              READ(temp,*) zeta%plot_sol_option, zeta%plot_mesh_option, zeta%plot_lines_option 
              PRINT("(A,I3)"), "zeta plot option = ", zeta%plot_sol_option 
              PRINT("(A,I3)"), "plot zeta mesh = ", zeta%plot_mesh_option  
              PRINT("(A,I3)"), "plot zeta contour lines = ", zeta%plot_lines_option               
              
            CASE (4)
              READ(temp,*) vel%plot_sol_option, vel%plot_mesh_option, vel%plot_lines_option
              PRINT("(A,I3)"), "velocity plot option = ", vel%plot_sol_option
              PRINT("(A,I3)"), "plot velocity mesh = ", vel%plot_mesh_option
              PRINT("(A,I3)"), "plot velocity contour lines = ", vel%plot_lines_option     
              
            CASE (5)
              READ(temp,*) bathy%plot_sol_option, bathy%plot_mesh_option, bathy%plot_lines_option
              PRINT("(A,I3)"), "bathymetry plot option = ", bathy%plot_sol_option 
              PRINT("(A,I3)"), "plot bathymetry mesh = ", bathy%plot_mesh_option
              PRINT("(A,I3)"), "plot bathymetry contour lines = ", bathy%plot_lines_option              
              
            CASE (6)
              READ(temp,*) mesh%plot_mesh_option
              mesh%plot_sol_option = mesh%plot_mesh_option
              PRINT("(A,I3)"), "mesh plot option = ", mesh%plot_mesh_option  
              
            CASE (7)
              mesh%el_label_option = TRIM(ADJUSTL(temp)) 
              PRINT("(A,A)"), "mesh plot element labels = ", mesh%el_label_option     
              IF (mesh%el_label_option == "file") THEN
                INQUIRE(file="element.label",exist=file_exists)
                IF (file_exists /= .true.) THEN
                  PRINT("(A)"), "Error: Mesh plot element label file not found"
                  STOP
                ENDIF              
              ENDIF
              
            CASE (8)
              mesh%nd_label_option = TRIM(ADJUSTL(temp)) 
              PRINT("(A,A)"), "mesh plot node labels = ", mesh%nd_label_option   
              IF (mesh%nd_label_option == "file") THEN
                INQUIRE(file="node.label",exist=file_exists)
                IF (file_exists /= .true.) THEN
                  PRINT("(A)"), "Error: Mesh plot node label file not found"
                  STOP
                ENDIF              
              ENDIF
              
            CASE (9)
              READ(temp,*) p_low,p_high,p_skip
              PRINT("(A,I3)"), "straight element plot nodes order = ", p_low,p_high,p_skip
              
            CASE (10)
              READ(temp,*) pc
              PRINT("(A,I3)"), "curved element plot nodes order = ", pc 
              
            CASE (11)
              zbox = TRIM(ADJUSTL(temp))
              IF (TRIM(ADJUSTL(zbox)) == "all") THEN
                xbox_min = -1d10
                xbox_max = 1d10
                ybox_min = -1d10
                ybox_max = 1d10
              ELSE
                READ(temp,*) xbox_min,xbox_max,ybox_min,ybox_max
              ENDIF
              PRINT("(A,A)"), "zoom box = ", zbox    
              
            CASE (12)
              READ(temp,*) figure_width
              figure_width = figure_width*72d0  
              PRINT("(A,F9.5)"), "figure width = ", figure_width
              
            CASE (13)
              READ(temp,*) snap_start,snap_end
              PRINT("(A,I5)"), "start snap = ", snap_start
              PRINT("(A,I5)"), "snap_end = ", snap_end
              
            CASE (14)
              cmap_file = TRIM(ADJUSTL(temp))
              PRINT("(A,A)"), "color map file = ", cmap_file     
              
            CASE (15)
              zeta%cscale_option = TRIM(ADJUSTL(temp)) 
              IF (TRIM(ADJUSTL(zeta%cscale_option)) == "auto-snap" .OR. TRIM(ADJUSTL(zeta%cscale_option)) == "auto-all") THEN
              
              ELSEIF (TRIM(ADJUSTL(zeta%cscale_option)) == "file") THEN
                INQUIRE(file='zeta.cscale',exist=file_exists)
                IF(file_exists /= .true.) THEN
                  PRINT("(A)"), "Error: zeta color scale file not found"
                  STOP
                ENDIF
              ELSE
                READ(temp,*) zeta%cscale_min,zeta%cscale_max
                zeta%cscale_option = "spec"
              ENDIF
              PRINT("(A,2(1x,F15.7))"), "zeta color scale = ", zeta%cscale_min,zeta%cscale_max   
              
            CASE (16)
              vel%cscale_option = TRIM(ADJUSTL(temp)) 
              IF (TRIM(ADJUSTL(vel%cscale_option)) == "auto-snap" .OR. TRIM(ADJUSTL(vel%cscale_option)) == "auto-all") THEN
              
              ELSEIF (TRIM(ADJUSTL(vel%cscale_option)) == "file") THEN
                INQUIRE(file='vel.cscale',exist=file_exists)
                IF(file_exists /= .true.) THEN                     
                  PRINT("(A)"), "Error: velocity color scale file not found"
                  STOP
                ENDIF
              ELSE
                READ(temp,*) vel%cscale_min,vel%cscale_max
                vel%cscale_option = "spec"
              ENDIF
              PRINT("(A,2(1x,F15.7))"), "velocity color scale = ", vel%cscale_min,vel%cscale_max  
              
            CASE (17)
              READ(temp,*) fontsize
              PRINT("(A,I3)"), "font size = ", fontsize 
              
            CASE (18)  
              READ(temp,*) nxtick
              PRINT("(A,I3)"), "number of x ticks = ", nxtick
              
            CASE (19)
              READ(temp,*) ytick
              IF (TRIM(ADJUSTL(ytick)) == 'auto') THEN
                nytick = 10000
              ELSE
                READ(temp,*) nytick
              ENDIF
              PRINT("(A,A)"), "number of y ticks = ", ytick
              
            CASE (20)
              READ(temp,*) nctick
              PRINT("(A,I3)"), "number of c ticks = ", nctick
              
            CASE (21)
              READ(temp,*) nxdec
              PRINT("(A,I3)"), "number of x decimals = ", nxdec
              
            CASE (22)
              READ(temp,*) nydec
              PRINT("(A,I3)"), "number of y decimals = ", nydec
              
            CASE (23) 
              READ(temp,*) ncdec
              PRINT("(A,I3)"), "number of c decimals = ", ncdec
              
            CASE (24)
              READ(temp,*) ntdec
              PRINT("(A,I3)"), "number of t decimals = ", ntdec
              
            CASE (25)
              READ(temp,*) frmt,rm_ps
              PRINT("(A,A)"), "additional file format = ", frmt
              PRINT("(A,I3)"), "remove PostScript files = ", rm_ps
              mesh%rm_ps = 0
              bathy%rm_ps = rm_ps
              zeta%rm_ps = rm_ps
              vel%rm_ps = rm_ps
            CASE (26)
              READ(temp,*) density
              PRINT("(A,A)"), "density of raster format = ", density 
              
            CASE (27)
              READ(temp,*) zeta%movie_flag
              vel%movie_flag = zeta%movie_flag
              PRINT("(A,I5)"), "movie flag = ", zeta%movie_flag
              
          END SELECT
            
        ENDIF
      
      
      ENDDO
      
      PRINT*, " "
      PRINT("(A,I5)"), "Lines skipped: ", skipped
      PRINT*, " "  
      
      IF (zeta%cscale_max < zeta%cscale_min) THEN
        tmp = zeta%cscale_min
        zeta%cscale_min = zeta%cscale_max
        zeta%cscale_max = tmp
        PRINT("(A)"), "Warning: zeta color scale max and min have been switched"        
      ENDIF         
      
      IF (vel%cscale_max < vel%cscale_min) THEN
        tmp = vel%cscale_min
        vel%cscale_min = vel%cscale_max
        vel%cscale_max = tmp
        PRINT("(A)"), "Warning: velocity color scale max and min have been switched"
      ENDIF
      
   
      IF (xbox_min >= xbox_max) THEN
        PRINT("(A)"), "Error: the order of x-values in zoom box is incorrect"
        STOP
      ENDIF
      
      IF (ybox_min >= ybox_max) THEN
        PRINT("(A)"), "Error: the order of -values in zoom box is incorrect"
        STOP
      ENDIF

      CLOSE(15)
      


      END SUBROUTINE  read_plot_input
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
      
      SUBROUTINE substitute_partial_path(path,replace,substitute)      
      
      IMPLICIT NONE
      
      CHARACTER(*) :: path
      CHARACTER(*) :: replace
      CHARACTER(*) :: substitute 
      
      CHARACTER(100) :: tmp
      INTEGER :: rlen
      INTEGER :: n
            
      rlen = LEN(TRIM(ADJUSTL(replace)))
      
      tmp = path(rlen+1:)      
      path = " "
      
      path = TRIM(ADJUSTL(substitute)) // TRIM(ADJUSTL(tmp))
      
      END SUBROUTINE substitute_partial_path