      SUBROUTINE read_plot_input()

      USE plot_globals, ONLY: input_path,cmap_file,zbox,ps,pc, &
                              xbox_min,xbox_max,ybox_min,ybox_max,figure_width, &
                              plot_mesh_option,plot_zeta_option,plot_vel_option,plot_bathy_option, &
                              snap_start,snap_end

      IMPLICIT NONE
      
      INTEGER, PARAMETER :: ninp = 11
      INTEGER :: i,n
      INTEGER :: inp_read,skipped
      CHARACTER(100) :: temp      
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
!             CASE ()
!               cscale = TRIM(ADJUSTL(temp) 
!               IF (cscale = "auto") THEN
!                 PRINT("(A,A)"), "color scale = ", cscale  
!               ELSE
!                 READ(temp,*) cscale_min,cscale_max
!                 PRINT("(A,2(e10.3))"), "color scale = " cscale_min,cscale_max
!               ENDIF

          END SELECT
            
        ENDIF
      
      
      ENDDO
      
      PRINT*, " "
      PRINT("(A,I5)"), "Lines skipped: ", skipped
      PRINT*, " "  
      

      CLOSE(15)
      


      END SUBROUTINE  read_plot_input