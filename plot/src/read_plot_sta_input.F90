      SUBROUTINE read_plot_sta_input()
      
      USE plot_globals, ONLY: zeta,vel, &
                              nadc,adc_sol,ndg,dg_sol, &
                              p_high,p_low,p_skip, &
                              xsta_min,xsta_max, &
                              xtime_min,xtime_max, &
                              snap_start,snap_end, &
                              sta_start,sta_end, &
                              xc_snap_opt,ts_sta_opt, &
                              plot_ts_sta_opt,plot_xc_snap_opt,plot_sta_loc_opt, &
                              line_width,figure_width,figure_height, &
                              nxtick,nytick,nxdec,nydec, &
                              fontsize,font, &
                              frmt,density, &
                              cmap_file, &
                              plot_error_opt,plot_scatter_opt
      USE read_dginp, ONLY: stations_file,grid_file,curve_file,ctp       

      IMPLICIT NONE
      
      INTEGER :: i
      LOGICAL :: file_exists  
      CHARACTER(100) :: temp      
      INTEGER :: rm_ps

      INQUIRE(file='plot_sta.inp',exist=file_exists)
      IF (file_exists == .FALSE.) THEN
        PRINT*, "plot_sta.inp file does not exist"
        STOP
      ENDIF  
      
      OPEN(unit=15,file='plot_sta.inp')
      
      READ(15,*) nadc
      PRINT*, nadc
      ALLOCATE(adc_sol(nadc))
      DO i = 1,nadc
        READ(15,"(A100)") temp
        adc_sol(i)%line = TRIM(ADJUSTL(temp))
        PRINT*, adc_sol(i)%line
      ENDDO
      
      READ(15,*) ndg
      PRINT*, ndg
      ALLOCATE(dg_sol(ndg))
      DO i = 1,ndg
        READ(15,"(A100)") temp
        dg_sol(i)%line = TRIM(ADJUSTL(temp))
        PRINT*, dg_sol(i)%line
      ENDDO
      
      READ(15,"(A100)") stations_file
      PRINT*, stations_file
      
      READ(15,*) plot_xc_snap_opt
      READ(15,"(A100)") temp
      IF (TRIM(ADJUSTL(temp)) == "all") THEN
        xc_snap_opt = TRIM(ADJUSTL(temp))
      ELSE IF (TRIM(ADJUSTL(temp)) == "file") THEN
        xc_snap_opt = TRIM(ADJUSTL(temp))      
        INQUIRE(file="xc_snap.list",exist=file_exists)
        IF (file_exists == .false.) THEN
          PRINT*, "xc_snap.list file not found"
          STOP
        ENDIF
      ELSE
        READ(temp,*) snap_start,snap_end
      ENDIF
      
      READ(15,*) xsta_min,xsta_max
      
      READ(15,*) plot_ts_sta_opt
      READ(15,"(A100)") temp
      IF (TRIM(ADJUSTL(temp)) == "all") THEN
        ts_sta_opt = TRIM(ADJUSTL(temp))
      ELSE IF (TRIM(ADJUSTL(temp)) == "file") THEN
        ts_sta_opt = TRIM(ADJUSTL(temp))      
        INQUIRE(file="ts_sta.list",exist=file_exists)
        IF (file_exists == .false.) THEN
          PRINT*, "ts_sta.list file not found"
          STOP
        ENDIF        
      ELSE
        READ(temp,*) sta_start,sta_end
      ENDIF      
      
      READ(15,*) xtime_min, xtime_max
      
      
      READ(15,*) plot_sta_loc_opt
      READ(15,"(A100)") grid_file
      READ(15,*) ctp
      READ(15,"(A100)") curve_file
      READ(15,*) p_high
        p_low = p_high
        p_skip = 1
      READ(15,*) plot_error_opt
      READ(15,*) plot_scatter_opt
      READ(15,"(A100)") cmap_file  
      READ(15,*) line_width
      READ(15,*) figure_width
         figure_width = figure_width*72d0
      READ(15,*) figure_height
         figure_height = figure_height*72d0
      READ(15,*) nxtick
      READ(15,*) nytick
      READ(15,*) nxdec
      READ(15,*) nydec
      READ(15,*) fontsize
      READ(15,*) font
      READ(15,*) frmt
      READ(15,*) density
      READ(15,*) rm_ps
        zeta%rm_ps = rm_ps
        vel%rm_ps = rm_ps
      
      
      CLOSE(15)

      
      zeta%plot_sol_option = 1
      vel%plot_sol_option = 1

      zeta%name = "zeta"
      vel%name = "vel"    
      
    
      
      RETURN
      END SUBROUTINE