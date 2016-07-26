      SUBROUTINE read_plot_input()

      USE plot_globals, ONLY: rp,input_path,cmap_file,ps,pc,frmt,density,rm_ps,make_movie, &
                              xbox_min,xbox_max,ybox_min,ybox_max,figure_width, &
                              plot_mesh_option,plot_zeta_option,plot_vel_option,plot_bathy_option, &
                              plot_zeta_mesh,plot_vel_mesh,plot_bathy_mesh, &
                              snap_start,snap_end, &
                              cscale_zeta,cscale_zeta_min,cscale_zeta_max, &
                              cscale_vel,cscale_vel_min,cscale_vel_max, &
                              cscale_zeta_unit,cscale_vel_unit, &
                              mesh_el_label,mesh_nd_label
                              
      USE plot_mod, ONLY: fontsize,nxtick,nytick,nctick, &
                          nxdec,nydec,ncdec,ntdec                              

      IMPLICIT NONE
      
      INTEGER, PARAMETER :: ninp = 26
      INTEGER :: i,n
      INTEGER :: inp_read,skipped
      CHARACTER(100) :: temp      
      CHARACTER(100) :: zbox  
      CHARACTER(100) :: ytick   
      CHARACTER(100) :: cscale 
      REAL(rp) :: tmp
      LOGICAL :: file_exists      
      
      cscale_zeta_unit = 20      
      cscale_vel_unit = 21      
      
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
              
            CASE (2)
              READ(temp,*) plot_zeta_option, plot_zeta_mesh
              PRINT("(A,I3)"), "zeta plot option = ", plot_zeta_option 
              PRINT("(A,I3)"), "plot zeta mesh = ", plot_zeta_option 
              
            CASE (3)
              READ(temp,*) plot_vel_option, plot_vel_mesh
              PRINT("(A,I3)"), "velocity plot option = ", plot_vel_option
              PRINT("(A,I3)"), "plot velocity mesh = ", plot_vel_mesh
              
            CASE (4)
              READ(temp,*) plot_bathy_option, plot_bathy_mesh
              PRINT("(A,I3)"), "bathymetry plot option = ", plot_bathy_option 
              PRINT("(A,I3)"), "plot bathymetry mesh = ", plot_bathy_mesh 
              
            CASE (5)
              READ(temp,*) plot_mesh_option
              PRINT("(A,I3)"), "mesh plot option = ", plot_mesh_option  
              
            CASE (6)
              mesh_el_label = TRIM(ADJUSTL(temp)) 
              PRINT("(A,A)"), "mesh plot element labels = ", mesh_el_label              
              
            CASE (7)
              mesh_nd_label = TRIM(ADJUSTL(temp)) 
              PRINT("(A,A)"), "mesh plot node labels = ", mesh_nd_label             
              
            CASE (8)
              READ(temp,*) ps
              PRINT("(A,I3)"), "straight element plot nodes order = ", ps 
              
            CASE (9)
              READ(temp,*) pc
              PRINT("(A,I3)"), "curved element plot nodes order = ", pc 
              
            CASE (10)
              zbox = TRIM(ADJUSTL(temp))
              IF (zbox == "all") THEN
                xbox_min = -1d10
                xbox_max = 1d10
                ybox_min = -1d10
                ybox_max = 1d10
              ELSE
                READ(temp,*) xbox_min,xbox_max,ybox_min,ybox_max
              ENDIF
              PRINT("(A,A)"), "zoom box = ", zbox    
              
            CASE (11)
              READ(temp,*) figure_width
              figure_width = figure_width*72d0  
              PRINT("(A,F9.5)"), "figure width = ", figure_width
              
            CASE (12)
              READ(temp,*) snap_start,snap_end
              PRINT("(A,I5)"), "start snap = ", snap_start
              PRINT("(A,I5)"), "snap_end = ", snap_end
              
            CASE (13)
              cmap_file = TRIM(ADJUSTL(temp))
              PRINT("(A,A)"), "color map file = ", cmap_file     
              
            CASE (14)
              cscale_zeta = TRIM(ADJUSTL(temp)) 
              IF (TRIM(ADJUSTL(cscale_zeta)) == "auto-snap" .OR. TRIM(ADJUSTL(cscale_zeta)) == "auto-all") THEN
              
              ELSEIF (TRIM(ADJUSTL(cscale_zeta)) == "file") THEN
                INQUIRE(file='zeta.cscale',exist=file_exists)
                IF(file_exists) THEN
                  OPEN(unit=cscale_zeta_unit,file="zeta.cscale")
                ELSE
                  PRINT("(A)"), "Error: zeta color scale file not found"
                ENDIF
              ELSE
                READ(temp,*) cscale_zeta_min,cscale_zeta_max
              ENDIF
              PRINT("(A,A)"), "zeta color scale = ", cscale_zeta   
              
            CASE (15)
              cscale_vel = TRIM(ADJUSTL(temp)) 
              IF (TRIM(ADJUSTL(cscale_vel)) == "auto-snap" .OR. TRIM(ADJUSTL(cscale_vel)) == "auto-all") THEN
              
              ELSEIF (TRIM(ADJUSTL(cscale_vel)) == "file") THEN
                INQUIRE(file='vel.cscale',exist=file_exists)
                IF(file_exists) THEN                     
                  OPEN(unit=cscale_vel_unit,file="vel.cscale")
                ELSE
                  PRINT("(A)"), "Error: velocity color scale file not found"
                  STOP
                ENDIF
              ELSE
                READ(temp,*) cscale_vel_min,cscale_vel_max
              ENDIF
              PRINT("(A,A)"), "velocity color scale = ", cscale_vel   
              
            CASE (16)
              READ(temp,*) fontsize
              PRINT("(A,I3)"), "font size = ", fontsize 
              
            CASE (17)  
              READ(temp,*) nxtick
              PRINT("(A,I3)"), "number of x ticks = ", nxtick
              
            CASE (18)
              READ(temp,*) ytick
              IF (TRIM(ADJUSTL(ytick)) == 'auto') THEN
                nytick = 10000
              ELSE
                READ(temp,*) nytick
              ENDIF
              PRINT("(A,A)"), "number of y ticks = ", ytick
              
            CASE (19)
              READ(temp,*) nctick
              PRINT("(A,I3)"), "number of c ticks = ", nctick
              
            CASE (20)
              READ(temp,*) nxdec
              PRINT("(A,I3)"), "number of x decimals = ", nxdec
              
            CASE (21)
              READ(temp,*) nydec
              PRINT("(A,I3)"), "number of y decimals = ", nydec
              
            CASE (22) 
              READ(temp,*) ncdec
              PRINT("(A,I3)"), "number of c decimals = ", ncdec
              
            CASE (23)
              READ(temp,*) ntdec
              PRINT("(A,I3)"), "number of t decimals = ", ntdec
              
            CASE (24)
              READ(temp,*) frmt,rm_ps
              PRINT("(A,A)"), "additional file format = ", frmt
              
            CASE (25)
              READ(temp,*) density
              PRINT("(A,A)"), "density of raster format = ", density 
              
            CASE (26)
              READ(temp,*) make_movie
              PRINT("(A,I5)"), "movie flag = ", make_movie     
              
          END SELECT
            
        ENDIF
      
      
      ENDDO
      
      PRINT*, " "
      PRINT("(A,I5)"), "Lines skipped: ", skipped
      PRINT*, " "  
      
      IF (cscale_zeta_max < cscale_zeta_min) THEN
        tmp = cscale_zeta_min
        cscale_zeta_min = cscale_zeta_max
        cscale_zeta_max = tmp
        PRINT("(A)"), "Warning: zeta color scale max and min have been switched"        
      ENDIF         
      
      IF (cscale_vel_max < cscale_vel_min) THEN
        tmp = cscale_vel_min
        cscale_vel_min = cscale_vel_max
        cscale_vel_max = tmp
        PRINT("(A)"), "Warning: velocity color scale max and min have been switched"
      ENDIF
      
   

      CLOSE(15)
      


      END SUBROUTINE  read_plot_input