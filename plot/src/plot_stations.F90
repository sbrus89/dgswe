      PROGRAM plot_sta


      USE plot_globals
      USE globals, ONLY: ne,nn,nverts,ndof,mndof,np,mnp,nnds,mnnds,nel_type,ect,xy,depth, &
                         nope,neta,obseg,obnds,nvel,nbou,fbseg,fbnds,bndxy,grid_name, &
                         el_type,elxy, &
                         nepn,epn,mnepn,ned,ged2el,ged2el,ged2led,ged2nn,ed_type,recv_edge, &
                         nied,iedn,nobed,obedn,nfbed,fbedn,nnfbed,nfbedn,nfbednn, &
                         psic,psiv, &
                         nsta,xysta
      USE grid_file_mod
      USE read_write_output
      USE plot_mod
      USE read_dginp, ONLY: stations_file,grid_file,curve_file,cb_file_exists,ctp     
      USE edge_connectivity_mod
      USE curvilinear_nodes_mod
      USE transformation
      USE shape_functions_mod   
      USE triangulation, ONLY: reference_element_delaunay  
      USE basis, ONLY: element_nodes,element_basis      
      USE station_plot_types, ONLY: cross_section_snap,time_series_station,error_station,error_scatter

      IMPLICIT NONE
      
      INTEGER :: k
      CHARACTER(100) :: out_direc
      CHARACTER(:), ALLOCATABLE :: filename
      INTEGER :: sta,run,ord
      REAL(rp) :: h,uu,vv
      REAL(rp) :: x_max,x_min,y_max,y_min
      INTEGER :: npt
      REAL(rp), DIMENSION(:), ALLOCATABLE :: xvec,yvec
      REAL(rp), DIMENSION(:,:), ALLOCATABLE :: xvec_err,yvec_err      
      INTEGER :: nrun 
      INTEGER :: nsta_plot
      INTEGER :: nsnap_max
      REAL(rp) :: tf
      INTEGER, DIMENSION(:), ALLOCATABLE :: xc_snaps
      INTEGER, DIMENSION(:), ALLOCATABLE :: ts_stas
      REAL(rp) :: err,max_err
      REAL(rp) :: yval
      REAL(rp) :: r2
      CHARACTER(3) :: maxavg_opt,absrel_opt
      CHARACTER(2) :: axis_label_option

      
      

!       nadc = 2
!       ALLOCATE(adc_sol(nadc))     
!       adc_sol(1)%line = "/home/sbrus/data-drive/galveston_spline_flux_fix/galveston_tri/adcirc/ESL0/"
!       adc_sol(2)%line = "/home/sbrus/data-drive/galveston_spline_flux_fix/galveston_tri_x64/adcirc/ESL0/"
!       ndg = 1
!       ALLOCATE(dg_sol(ndg))
!       dg_sol(1)%line = "/home/sbrus/data-drive/galveston_spline_flux_fix/galveston_tri/p3/ctp3/hbp3/"

!       nadc = 1
!       ALLOCATE(adc_sol(nadc))
!       adc_sol(1)%line = "/home/sbrus/data-drive/galveston_spline_flux_fix/galveston_tri_x64/adcirc/ESL0/"
!       ndg = 2
!       ALLOCATE(dg_sol(ndg))
!       dg_sol(1)%line = "/home/sbrus/data-drive/galveston_spline_flux_fix/galveston_tri/p3/ctp3/hbp3/"
!       dg_sol(2)%line = "/home/sbrus/data-drive/galveston_spline_flux_fix/galveston_tri/p3/ctp1/hbp1/"      

      
      CALL read_plot_sta_input()
                 
      nsnap_max = 2000  
       
      CALL read_stations(0,stations_file,1,nsta,xysta)   
      nsta_plot = nsta      
      

      
      CALL read_colormap(cmap_file)
      
      ALLOCATE(zeta%sta_val(nsta,nsnap_max+1,ndg+nadc),zeta%t_sta(nsnap_max+1,ndg+nadc))         
      ALLOCATE(vel%sta_val(nsta,nsnap_max+1,ndg+nadc),vel%t_sta(nsnap_max+1,ndg+nadc))           
      
      
      nrun = 0            
      
      DO i = 1,nadc
      
        out_direc = adc_sol(i)%line    
        
        nrun = nrun + 1        
      
        IF (zeta%plot_sol_option == 1 .or. vel%plot_sol_option == 1) THEN
          PRINT("(A)"), "Reading zeta solution..."          
          CALL read_fort6163(out_direc,"61",t,eta,nsnap_Z) 
        ENDIF     
        IF (vel%plot_sol_option == 1) THEN
          PRINT("(A)"), "Reading u and v solutions..."       
          CALL read_fort6264(out_direc,"62",t,uu2,vv2,nsnap_Qx)                   
        ELSE 
          ALLOCATE(Qx(1,1,1),Qy(1,1,1))
        ENDIF   
      
        DO sta = 1,nsta
          zeta%sta_val(sta,1,nrun) = 0d0
        ENDDO            
        zeta%t_sta(1,nrun) = 0d0
        DO snap  = 1,nsnap_Z
          DO sta = 1,nsta
            zeta%sta_val(sta,snap+1,nrun) = eta(sta,snap)
          ENDDO
          zeta%t_sta(snap+1,nrun) = t(snap)/86400d0
        ENDDO
      
      
        DO sta = 1,nsta
          vel%sta_val(sta,1,nrun) = 0d0
        ENDDO      
        vel%t_sta(1,nrun) = 0d0
        DO snap = 1,nsnap_Qx
          DO sta = 1,nsta
            vel%sta_val(sta,snap+1,nrun) = sqrt(uu2(sta,snap)**2 + vv2(sta,snap)**2)
          ENDDO
          vel%t_sta(snap+1,nrun) = t(snap)/86400d0
        ENDDO      
  
      ENDDO 
 
      DEALLOCATE(t)
 
 
      DO i = 1,ndg
      
        out_direc = dg_sol(i)%line          
      
        nrun = nrun + 1
        
        IF (zeta%plot_sol_option == 1 .or. vel%plot_sol_option == 1) THEN
          PRINT("(A)"), "Reading zeta solution..."          
          CALL read_solution_full(out_direc,"Z.sta","N",t,Z,nsnap_Z) 
        ENDIF
        IF (vel%plot_sol_option == 1) THEN
          PRINT("(A)"), "Reading Qx and Qy solutions..."       
          CALL read_solution_full(out_direc,"Qx.sta","N",t,Qx,nsnap_Qx)        
          CALL read_solution_full(out_direc,"Qy.sta","N",t,Qy,nsnap_Qy) 
        ELSE 
          ALLOCATE(Qx(1,1,1),Qy(1,1,1))
        ENDIF  
        IF (vel%plot_sol_option == 1) THEN
          PRINT("(A)"), "Reading bathymetry solution..."          
          CALL read_solution_full(out_direc,"hb.sta","N",t,hb,nsnap_hb)  
        ENDIF   
 
   
        DO snap  = 1,nsnap_Z
          DO sta = 1,nsta
            zeta%sta_val(sta,snap,nrun) = Z(sta,1,snap)
          ENDDO
          zeta%t_sta(snap,nrun) = t(snap)/86400d0
        ENDDO
      

        DO snap = 1,nsnap_Qx
          DO sta = 1,nsta
            h = Z(sta,1,snap) + hb(sta,1,1)
            uu = Qx(sta,1,snap)/h
            vv = Qy(sta,1,snap)/h
            vel%sta_val(sta,snap,nrun) = sqrt(uu**2 + vv**2)
          ENDDO
          vel%t_sta(snap,nrun) = t(snap)/86400d0
        ENDDO
      
      ENDDO
      
      
      
      
      
      nsta_plot = xsta_max-xsta_min+1
      
      IF (xc_snap_opt == "all") THEN
        nsnap = nsnap_Z
        ALLOCATE(xc_snaps(nsnap))
        DO i = 1,nsnap
          xc_snaps(i) = i
        ENDDO
      ELSE IF (xc_snap_opt == "file") THEN
        OPEN(unit=10,file="xc_snap.list")
        READ(10,*) nsnap
        ALLOCATE(xc_snaps(nsnap))
        DO i = 1,nsnap
          READ(10,*) xc_snaps(i)
        ENDDO
        CLOSE(10)
      ELSE
        nsnap = snap_end-snap_start+1
        ALLOCATE(xc_snaps(nsnap))
        DO i = 1,nsnap
          xc_snaps(i) = snap_start + (i-1)
        ENDDO
      ENDIF
      
      
          
      
      IF (ts_sta_opt == "all") THEN
        ALLOCATE(ts_stas(nsta))
        DO i = 1,nsta
          ts_stas(i) = i
        ENDDO
      ELSE IF (ts_sta_opt == "file") THEN
        OPEN(unit=10,file="ts_sta.list")
        READ(10,*) nsta
        ALLOCATE(ts_stas(nsta))
        DO i = 1,nsta
          READ(10,*) ts_stas(i)
        ENDDO
        CLOSE(10)
      ELSE
        nsta = sta_end-sta_start+1
        ALLOCATE(ts_stas(nsta))
        DO i = 1,nsta
          ts_stas(i) = sta_start + (i-1)
        ENDDO
      ENDIF      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      IF (plot_xc_snap_opt > 0) THEN
      
        CALL cross_section_snap(zeta,nrun,nsnap,xc_snaps,xsta_min,xsta_max)      
        CALL cross_section_snap(vel,nrun,nsnap,xc_snaps,xsta_min,xsta_max)

      ENDIF
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!            
      
      IF (plot_error_opt == 1) THEN

        CALL error_station(vel,nrun,nsta,ts_stas,xsta_min,xsta_max,nsnap,xc_snaps,"max","abs")   
        CALL error_station(vel,nrun,nsta,ts_stas,xsta_min,xsta_max,nsnap,xc_snaps,"max","rel") 
        CALL error_station(vel,nrun,nsta,ts_stas,xsta_min,xsta_max,nsnap,xc_snaps,"avg","abs")  
        CALL error_station(vel,nrun,nsta,ts_stas,xsta_min,xsta_max,nsnap,xc_snaps,"avg","rel")   
      
        CALL error_station(zeta,nrun,nsta,ts_stas,xsta_min,xsta_max,nsnap,xc_snaps,"max","abs")   
        CALL error_station(zeta,nrun,nsta,ts_stas,xsta_min,xsta_max,nsnap,xc_snaps,"max","rel") 
        CALL error_station(zeta,nrun,nsta,ts_stas,xsta_min,xsta_max,nsnap,xc_snaps,"avg","abs")  
        CALL error_station(zeta,nrun,nsta,ts_stas,xsta_min,xsta_max,nsnap,xc_snaps,"avg","rel")                            
      
      ENDIF
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      IF (plot_scatter_opt == 1) THEN      
      
        CALL error_scatter(vel,nrun,nsta,ts_stas,nsnap,xc_snaps)
        CALL error_scatter(zeta,nrun,nsta,ts_stas,nsnap,xc_snaps)          

      ENDIF
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
      
      IF (plot_ts_sta_opt > 0) THEN
      
        CALL time_series_station(zeta,nrun,nsta,ts_stas,xtime_min,xtime_max,nsnap_Z)  
        CALL time_series_station(vel,nrun,nsta,ts_stas,xtime_min,xtime_max,nsnap_Z)              
      
      ENDIF
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IF (plot_sta_loc_opt > 0) THEN
      CALL sizes()

      CALL read_header(0,TRIM(grid_file),grid_name,ne,nn)
      CALL read_coords(nn,xy,depth)
      CALL read_connectivity(ne,ect,el_type)
      CALL init_element_coordinates(ne,ctp,el_type,nverts,xy,ect,elxy)                  
      CALL read_open_boundaries(nope,neta,obseg,obnds)            
      CALL read_flow_boundaries(nbou,nvel,fbseg,fbnds)                       
      CALL read_curve_file(0,curve_file,ctp,nbou,xy,bndxy,cb_file_exists)   
      
      CALL elements_per_node(ne,nn,nverts,el_type,ect,nepn,mnepn,epn)       
      CALL find_edge_pairs(ne,nverts,el_type,ect,nepn,epn,ned,ged2el,ged2nn,ged2led)      
      CALL find_interior_edges(ned,ged2el,nied,iedn,ed_type,recv_edge)      
      CALL find_open_edges(nope,obseg,obnds,ged2nn,nobed,obedn,ed_type,recv_edge)            
      CALL find_flow_edges(nbou,fbseg,fbnds,ged2nn,nnfbed,nfbedn,nfbednn,nfbed,fbedn,recv_edge,ed_type)          
      
      CALL shape_functions_linear_at_ctp(nel_type,np,psiv)                   
      CALL eval_coordinates_curved(ctp,nnds,nverts,el_type,xy,ect,fbseg,fbnds, &
                                   nnfbed,nfbedn,nfbednn,ged2el,ged2led, &
                                   psiv,bndxy,elxy)  
                                   
      nord = (p_high-p_low+1)/p_skip
      
      ALLOCATE(r(mnpp,nel_type*nord),s(mnpp,nel_type*nord))      
      ALLOCATE(psic(mnnds,mnpp,nel_type*nord))
      ALLOCATE(rect(3,3*mnpp,nel_type*nord))      
      ALLOCATE(nptri(nel_type*nord),npplt(nel_type*nord),pplt(nel_type*nord))
      DO et = 1,nel_type 
              
        DO ord = 1,nord
          i = (et-1)*nord+ord
          pplt(i) = (ord-1)*p_skip+p_low
          CALL element_nodes(et,space,pplt(i),npplt(i),r(:,i),s(:,i))                  
          CALL shape_functions_area_eval(et,np(et),nnd,npplt(i),r(:,i),s(:,i),psic(:,:,i))  
!           CALL reference_element_delaunay(et,pplt(i),npplt(i),r(:,i),s(:,i),nptri(i),rect(:,:,i))        
        ENDDO         
        
      ENDDO                                    
           
      ALLOCATE(xyplt(mnpp,ne,2))
      DO el = 1,ne      
        et = el_type(el)                          
        nnd = nnds(et)
        i = (et-1)*nord+nord
        npts = npplt(i)
        DO pt = 1,npts              
          CALL element_transformation(nnd,elxy(:,el,1),elxy(:,el,2),psic(:,pt,i),xpt,ypt)           
          xyplt(pt,el,1) = xpt
          xyplt(pt,el,2) = ypt
        ENDDO
      ENDDO                                          
      
      x_min = MINVAL(xy(1,:))
      x_max = MAXVAL(xy(1,:))
      y_min = MINVAL(xy(2,:))
      y_max = MAXVAL(xy(2,:))
      
      xbox_min = x_min      
      xbox_max = x_max
      ybox_min = y_min      
      ybox_max = y_max      
      
      figure_height = -1d0
      CALL scale_factors(figure_width,figure_height,x_min,x_max,y_min,y_max,ax,bx,ay,by)   

      ALLOCATE(vel%el_plt(ne))  
      ALLOCATE(el_in(el))
      DO el = 1,ne
        et = el_type(el)
        vel%el_plt(el) = (et-1)*nord + nord
        el_in(el) = 1
      ENDDO             

      DO i = 1,nsta
      
        sta = ts_stas(i)      
        WRITE(snap_char,"(I4.4)") sta           
        filename = "sta_loc_"//snap_char      

        CALL write_psheader(filename//".ps",vel%ps_unit)  
        CALL plot_background(vel%ps_unit,1d0,1d0,1d0)        
!         CALL plot_mesh(vel%ps_unit,ne,nverts,vel%el_plt,pplt,el_type,el_in,xy,ect,xyplt)   
        CALL plot_boundaries(vel%ps_unit,nverts,vel%el_plt,pplt,nnfbed,nfbedn,ged2el,ged2led,el_type,el_in,ect,xy,xyplt)        
        CALL plot_boundaries(vel%ps_unit,nverts,vel%el_plt,pplt,nfbed,fbedn,ged2el,ged2led,el_type,el_in,ect,xy,xyplt)         
        CALL plot_stations(vel%ps_unit,1,nsta_plot,xysta,sta)
        CALL close_ps(filename,vel%ps_unit)
        CALL convert_ps(filename,frmt,density,vel%rm_ps)           
                
      ENDDO
      
      ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 


      END PROGRAM plot_sta