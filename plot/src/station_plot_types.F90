      MODULE station_plot_types

      USE plot_globals, ONLY: rp,plot_type,colors,line_width
      USE plot_mod, ONLY: prep_station_plot,plot_xy,finish_station_plot,plot_xy_scatter
      USE stats, ONLY: r_squared

      IMPLICIT NONE
      
      CHARACTER(:), ALLOCATABLE :: filename
      CHARACTER(4) :: snap_char 
      CHARACTER(2) :: axis_label_type        

      CONTAINS
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       

      SUBROUTINE cross_section_snap(fig,nrun,nsnap,xc_snaps,xsta_min,xsta_max)            
      
      IMPLICIT NONE
      
      TYPE(plot_type), INTENT(INOUT) :: fig
      INTEGER, INTENT(IN) :: nrun
      INTEGER, INTENT(IN) :: nsnap
      INTEGER, DIMENSION(:), INTENT(IN) :: xc_snaps
      INTEGER, INTENT(IN) :: xsta_min
      INTEGER, INTENT(IN) :: xsta_max
      
      INTEGER :: i,j,run,pt,snap      
      INTEGER :: nsta_plot
      REAL(rp) :: x_min,x_max,y_min,y_max
      REAL(rp), DIMENSION(:), ALLOCATABLE :: xvec,yvec
      
      nsta_plot = xsta_max-xsta_min+1
      
      x_min = real(xsta_min,rp)
      x_max = real(xsta_max,rp)
      
      ALLOCATE(xvec(nsta_plot),yvec(nsta_plot))
      
      DO i = 1,nsnap
      
        snap = xc_snaps(i)
        PRINT*, fig%name//" snap",snap
      
        y_min = MINVAL(fig%sta_val(xsta_min:xsta_max,snap,1:nrun))  
        y_max = MAXVAL(fig%sta_val(xsta_min:xsta_max,snap,1:nrun)) 
        
        IF (fig%name == "vel") THEN
          y_min = 0d0
        ELSE
          y_min = y_min - abs((y_max-y_min)*.05d0)
        ENDIF
        y_max = y_max + abs((y_max-y_min)*.05d0)    
        
        
        IF (abs(y_min) < 1d-8 .and. abs(y_max) < 1d-8) THEN
          CYCLE 
        ENDIF
        
        WRITE(snap_char,"(I4.4)") snap           
        filename = TRIM(ADJUSTL(fig%name))//"_xs_"//snap_char   
        
        IF (fig%name == "vel") THEN
          axis_label_type = "vs"
        ELSE
          axis_label_type = "zs"
        ENDIF
        
        CALL prep_station_plot(fig,axis_label_type,snap,x_min,x_max,y_min,y_max,filename)        
        
        DO run = 1,nrun
          j = 1
          DO pt = xsta_min,xsta_max
            xvec(j) = pt
            yvec(j) = fig%sta_val(pt,snap,run)
            j = j + 1
          ENDDO
        
          CALL plot_xy(fig%ps_unit,nsta_plot,xvec,yvec,colors(run,:),line_width)
        ENDDO
        
        CALL finish_station_plot(fig,filename)
      ENDDO      
      
      
      RETURN
      END SUBROUTINE cross_section_snap
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        

      SUBROUTINE time_series_station(fig,nrun,nsta,ts_stas,xtime_min,xtime_max,nsnap_tot)
      
      IMPLICIT NONE
      
      TYPE(plot_type), INTENT(INOUT) :: fig
      INTEGER, INTENT(IN) :: nrun
      INTEGER, INTENT(IN) :: nsta
      INTEGER, DIMENSION(:), INTENT(IN) :: ts_stas      
      REAL(rp), INTENT(IN) :: xtime_min
      REAL(rp), INTENT(IN) :: xtime_max
      INTEGER :: nsnap_tot
      
      INTEGER :: i,j,run,pt,sta
      REAL(rp) :: x_min,x_max,y_min,y_max
      REAL(rp) :: tf
      INTEGER :: snap_start,snap_end
      INTEGER :: nsnap
      REAL(rp), DIMENSION(:), ALLOCATABLE :: xvec,yvec
      
      
      x_min = xtime_min
      x_max = xtime_max
      
      tf = fig%t_sta(nsnap_tot,1)
      
      snap_start = INT(((nsnap_tot-1)/tf)*xtime_min)+1
      snap_end = INT(((nsnap_tot-1)/tf)*xtime_max)
      
      nsnap = snap_end-snap_start+1       
      
      ALLOCATE(xvec(nsnap),yvec(nsnap))
      
      DO i = 1,nsta
      
        sta = ts_stas(i)
        PRINT*, fig%name//" station",sta      
      
        y_min = MINVAL(fig%sta_val(sta,snap_start:snap_end,1:nrun))  
        y_max = MAXVAL(fig%sta_val(sta,snap_start:snap_end,1:nrun))         
        
        IF (fig%name == "vel") THEN
          y_min = 0d0
        ELSE
          y_min = y_min - abs((y_max-y_min)*.05d0)
        ENDIF
        y_max = y_max + abs((y_max-y_min)*.05d0)            

        
        WRITE(snap_char,"(I4.4)") sta           
        filename = TRIM(ADJUSTL(fig%name))//"_ts_"//snap_char
        
        IF (fig%name == "vel") THEN
          axis_label_type = "vt"
        ELSE
          axis_label_type = "zt"
        ENDIF        
        
        CALL prep_station_plot(fig,axis_label_type,sta,x_min,x_max,y_min,y_max,filename)        
        
        DO run = 1,nrun
          j = 1
          DO pt = snap_start,snap_end
            xvec(j) = fig%t_sta(pt,run)
            yvec(j) = fig%sta_val(sta,pt,run)
            j = j+1
          ENDDO
        
          CALL plot_xy(fig%ps_unit,nsnap,xvec,yvec,colors(run,:),line_width)
        ENDDO
        
        CALL finish_station_plot(fig,filename)
      ENDDO      
      
      RETURN
      END SUBROUTINE time_series_station

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      SUBROUTINE error_station(fig,nrun,nsta,ts_stas,xsta_min,xsta_max,nsnap,xc_snaps,maxavg_opt,absrel_opt)
      
      IMPLICIT NONE
      
      TYPE(plot_type), INTENT(INOUT) :: fig
      INTEGER, INTENT(IN) :: nrun
      INTEGER, INTENT(IN) :: nsta
      INTEGER, DIMENSION(:), INTENT(IN) :: ts_stas
      INTEGER, INTENT(IN) :: xsta_min
      INTEGER, INTENT(IN) :: xsta_max
      INTEGER, INTENT(IN) :: nsnap
      INTEGER, DIMENSION(:), INTENT(IN) :: xc_snaps
      CHARACTER(3), INTENT(IN) :: maxavg_opt
      CHARACTER(3), INTENT(IN) :: absrel_opt
      
      INTEGER :: i,j,run,sta,snap
      INTEGER :: nsta_plot
      REAL(rp), DIMENSION(:,:), ALLOCATABLE :: xvec_err,yvec_err
      REAL(rp) :: max_err,err
      REAL(rp) :: x_min,x_max,y_min,y_max
      
      nsta_plot = xsta_max-xsta_min+1               
      
      ALLOCATE(xvec_err(nsta_plot,nrun-1),yvec_err(nsta_plot,nrun-1))
    
      DO run = 1,nrun-1
      DO i = 1,nsta_plot
        sta = ts_stas(i)
        
        max_err = -1d10        
        err = 0d0
        DO j = 1,nsnap
          snap = xc_snaps(j)
          
          IF (maxavg_opt == "avg") THEN
            err = err + abs(fig%sta_val(sta,snap,1)-fig%sta_val(sta,snap,run+1)) 
          ELSE IF (maxavg_opt == "max") THEn
            err = abs(fig%sta_val(sta,snap,1)-fig%sta_val(sta,snap,run+1))                        
            IF (err > max_err) THEN
              max_err = err
            ENDIF          
          ENDIF
        ENDDO
        
        IF (maxavg_opt == "avg") THEN
          err = err/real(nsnap,rp)
        ELSE IF (maxavg_opt == "max") THEN
          err = max_err    
        ENDIF

        
        
        IF (absrel_opt == "rel") THEN
          y_max = -1d10
          y_min = 1d10
          DO j = 1,nsnap
            snap = xc_snaps(j)          
          
            IF (fig%sta_val(sta,snap,1) > y_max) THEN
              y_max = fig%sta_val(sta,snap,1)
            ENDIF
          
            IF (fig%sta_val(sta,snap,1) < y_min) THEN
              y_min = fig%sta_val(sta,snap,1)
            ENDIF          
          ENDDO

        
          err = err/(y_max - y_min)
        ENDIF
        
        xvec_err(i,run) = sta
        yvec_err(i,run) = err              
        
      ENDDO
      ENDDO
     
      
      y_min = 0d0  
      y_max = MAXVAL(yvec_err(1:nsta,1:nrun-1)) 
      y_max = y_max + y_max*.05
        
      x_min = real(xsta_min,rp)
      x_max = real(xsta_max,rp)        
      filename = TRIM(fig%name)//"_"//maxavg_opt//"_"//absrel_opt//"_error"
      
      IF (absrel_opt == "rel") THEN
        IF (fig%name == "vel") THEN
          axis_label_type = "vr"      
        ELSE IF (fig%name == "zeta") THEN
          axis_label_type = "zr"          
        ENDIF
      ELSE IF (absrel_opt == "abs") THEN
        IF (fig%name == "vel") THEN
          axis_label_type = "va"      
        ELSE IF (fig%name == "zeta") THEN
          axis_label_type = "za"          
        ENDIF
      ENDIF
      

      CALL prep_station_plot(fig,axis_label_type,nsnap,x_min,x_max,y_min,y_max,filename) 
      DO run = 1,nrun-1
        CALL plot_xy(fig%ps_unit,nsta_plot,xvec_err(:,run),yvec_err(:,run),colors(run,:),line_width)      
      ENDDO
      CALL finish_station_plot(fig,filename)           
      
      RETURN
      END SUBROUTINE error_station

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 


      SUBROUTINE error_scatter(fig,nrun,nsta,ts_stas,nsnap,xc_snaps)
      
      IMPLICIT NONE
      
      TYPE(plot_type), INTENT(INOUT) :: fig
      INTEGER, INTENT(IN) :: nrun
      INTEGER, INTENT(IN) :: nsta
      INTEGER, DIMENSION(:), INTENT(IN) :: ts_stas
      INTEGER, INTENT(IN) :: nsnap
      INTEGER, DIMENSION(:), INTENT(IN) :: xc_snaps
      
      
      INTEGER :: i,j,k,run,sta,snap
      REAL(rp), DIMENSION(:,:), ALLOCATABLE :: xvec_err,yvec_err
      REAL(rp) :: x_min,x_max,y_min,y_max      
      REAL(rp) :: r2
      
      ALLOCATE(xvec_err(nsta*nsnap,nrun-1),yvec_err(nsta*nsnap,nrun-1))
    
      DO run = 1,nrun-1
      k = 0
      DO i = 1,nsta
        sta = ts_stas(i)
        
        DO j = 1,nsnap
          snap = xc_snaps(j)
          k = k+1
          xvec_err(k,run) = fig%sta_val(sta,snap,1)
          yvec_err(k,run) = fig%sta_val(sta,snap,run+1) 
        ENDDO

      ENDDO
      ENDDO
     
      IF (fig%name == "vel") THEN
        y_min = 0d0   
      ELSE IF (fig%name == "zeta") THEN
        y_min = MINVAL(yvec_err(1:nsta*nsnap,1:nrun-1))  
      ENDIF
      y_max = MAXVAL(yvec_err(1:nsta*nsnap,1:nrun-1))       
        
      IF (fig%name == "vel") THEN  
        x_min = 0d0
      ELSE IF (fig%name == "zeta") THEN
        x_min = MINVAL(xvec_err(1:nsta*nsnap,1:nrun-1)) 
      ENDIF
      x_max = MAXVAL(xvec_err(1:nsta*nsnap,1:nrun-1))       
      
      IF (y_max > x_max) THEN
        x_max = y_max
      ELSE IF(x_max > y_max) THEN
       y_max = x_max
      ENDIF
      
      filename = TRIM(fig%name)//"_scatter_error"
      
      IF (fig%name == "vel") THEN
        axis_label_type = "vc"
      ELSE IF (fig%name == "zeta") THEN
        axis_label_type = "zc"
      ENDIF
      
      
      CALL prep_station_plot(fig,axis_label_type,snap,x_min,x_max,y_min,y_max,filename) 
      DO run = 1,nrun-1
        CALL plot_xy_scatter(fig%ps_unit,nsta*nsnap,xvec_err(:,run),yvec_err(:,run),colors(run,:),line_width)      
        CALL r_squared(nsta*nsnap,yvec_err(:,run),xvec_err(:,run),r2)    
        PRINT*,r2
      ENDDO
      CALL plot_xy_scatter(fig%ps_unit,nsta*nsnap,xvec_err(:,1),xvec_err(:,1),colors(nrun,:),line_width)       
      CALL finish_station_plot(fig,filename)          
      
      RETURN
      END SUBROUTINE error_scatter

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 



      END MODULE station_plot_types