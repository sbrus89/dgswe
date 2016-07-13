      SUBROUTINE read_plot_input()

      USE plot_globals, ONLY: rp,input_path,cmap_file,ps,pc,frmt,density,rm_ps,make_movie, &
                              xbox_min,xbox_max,ybox_min,ybox_max,figure_width, &
                              plot_mesh_option,plot_zeta_option,plot_vel_option,plot_bathy_option, &
                              snap_start,snap_end, &
                              cscale_zeta,cscale_zeta_min,cscale_zeta_max, &
                              cscale_vel,cscale_vel_min,cscale_vel_max
      USE plot_mod, ONLY: fontsize,nxtick,nytick,nctick, &
                          nxdec,nydec,ncdec,ntdec                              

      IMPLICIT NONE
      
      INTEGER, PARAMETER :: ninp = 24
      INTEGER :: i,n
      INTEGER :: inp_read,skipped
      CHARACTER(100) :: temp      
      CHARACTER(100) :: zbox  
      CHARACTER(100) :: ytick   
      CHARACTER(100) :: cscale 
      REAL(rp) :: tmp
      LOGICAL :: file_exists      
      
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
              READ(temp,*) plot_zeta_option
              PRINT("(A,I3)"), "zeta plot option = ", plot_zeta_option 
            CASE (3)
              READ(temp,*) plot_vel_option
              PRINT("(A,I3)"), "velocity plot option = ", plot_vel_option
            CASE (4)
              READ(temp,*) plot_bathy_option
              PRINT("(A,I3)"), "bathymetry plot option = ", plot_bathy_option              
            CASE (5)
              READ(temp,*) plot_mesh_option
              PRINT("(A,I3)"), "mesh plot option = ", plot_mesh_option              
            CASE (6)
              READ(temp,*) ps
              PRINT("(A,I3)"), "straight element plot nodes order = ", ps 
            CASE (7)
              READ(temp,*) pc
              PRINT("(A,I3)"), "curved element plot nodes order = ", pc               
            CASE (8)
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
            CASE (9)
              READ(temp,*) figure_width
              figure_width = figure_width*72d0  
              PRINT("(A,F9.5)"), "figure width = ", figure_width              
            CASE (10)
              READ(temp,*) snap_start,snap_end
            CASE (11)
              cmap_file = TRIM(ADJUSTL(temp))
              PRINT("(A,A)"), "color map file = ", cmap_file                
            CASE (12)
              cscale_zeta = TRIM(ADJUSTL(temp)) 
              IF (TRIM(ADJUSTL(cscale_zeta)) == "auto-snap" .OR. TRIM(ADJUSTL(cscale_zeta)) == "auto-all") THEN
              
              ELSE
                READ(temp,*) cscale_zeta_min,cscale_zeta_max
              ENDIF
              PRINT("(A,A)"), "zeta color scale = ", cscale_zeta      
            CASE (13)
              cscale_vel = TRIM(ADJUSTL(temp)) 
              IF (TRIM(ADJUSTL(cscale_vel)) == "auto-snap" .OR. TRIM(ADJUSTL(cscale_vel)) == "auto-all") THEN
              
              ELSE
                READ(temp,*) cscale_vel_min,cscale_vel_max
              ENDIF
              PRINT("(A,A)"), "velocity color scale = ", cscale_vel                  
            CASE (14)
              READ(temp,*) fontsize
              PRINT("(A,I3)"), "font size = ", fontsize 
            CASE (15)  
              READ(temp,*) nxtick
              PRINT("(A,I3)"), "number of x ticks = ", nxtick
            CASE (16)
              READ(temp,*) ytick
              IF (TRIM(ADJUSTL(ytick)) == 'auto') THEN
                nytick = 10000
              ELSE
                READ(temp,*) nytick
              ENDIF
              PRINT("(A,A)"), "number of y ticks = ", ytick
            CASE (17)
              READ(temp,*) nctick
              PRINT("(A,I3)"), "number of c ticks = ", nctick
            CASE (18)
              READ(temp,*) nxdec
              PRINT("(A,I3)"), "number of x decimals = ", nxdec
            CASE (19)
              READ(temp,*) nydec
              PRINT("(A,I3)"), "number of y decimals = ", nydec
            CASE (20) 
              READ(temp,*) ncdec
              PRINT("(A,I3)"), "number of c decimals = ", ncdec
            CASE (21)
              READ(temp,*) ntdec
              PRINT("(A,I3)"), "number of t decimals = ", ntdec
            CASE (22)
              READ(temp,*) frmt,rm_ps
              PRINT("(A,A)"), "additional file format = ", frmt
            CASE (23)
              READ(temp,*) density
              PRINT("(A,A)"), "density of raster format = ", density     
            CASE (24)
              READ(temp,*) make_movie
              PRINT("(A,A)"), "movie flag = ", make_movie                 
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