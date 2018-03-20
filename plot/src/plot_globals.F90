      MODULE plot_globals
      
      USE globals, ONLY: rp
      
      INTEGER :: sol_diff_option
      INTEGER :: ho_diff_option
      CHARACTER(100) :: input_path
      CHARACTER(100) :: input_path2      
      CHARACTER(100) :: cmap_file
      CHARACTER(100) :: tsnap_spec
      CHARACTER(3) :: frmt
      CHARACTER(3) :: density
      CHARACTER(5) :: font
      
      INTEGER :: plot_google_map
      
      REAL(rp), PARAMETER :: max_init = -1d10
      REAL(rp), PARAMETER :: min_init =  1d10
      
      INTEGER :: substitute_path
      CHARACTER(100) :: replace_path
      CHARACTER(100) :: sub_path
      
!       INTEGER :: coord_sys
!       REAL(rp) :: slam0,sphi0      
      REAL(rp) :: xmin,xmax
      REAL(rp) :: ymin,ymax
      REAL(rp) :: rmin,rmax
      REAL(rp) :: smin,smax
      
      REAL(rp) :: lr_margin 
      REAL(rp) :: top_margin
      REAL(rp) :: cscale_width
      REAL(rp) :: axes_width
      REAL(rp) :: axes_height
      REAL(rp) :: rmin_axes,rmax_axes
      REAL(rp) :: smin_axes,smax_axes    
      REAL(rp) :: rmin_cbar,rmax_cbar
      REAL(rp) :: smin_cbar,smax_cbar   
      REAL(rp) :: rmin_tbar,rmax_tbar
      REAL(rp) :: smin_tbar,smax_tbar   
      REAL(rp) :: rmin_scale,rmax_scale
      REAL(rp) :: smin_scale,smax_scale
      REAL(rp) :: ax,bx  
      REAL(rp) :: ay,by      
      
      INTEGER :: ncolors      
      REAL(rp), DIMENSION(:,:), ALLOCATABLE :: colors        
      
      REAL(rp) :: dash      
      REAL(rp) :: xticklabel_pad,yticklabel_pad
      REAL(rp) :: xlabel_pad,ylabel_pad  
      REAL(rp) :: cticklabel_pad
      REAL(rp) :: clabel_pad     
      REAL(rp) :: scale_pad
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
      
      INTEGER :: scale_flag
      INTEGER :: scale_label
      CHARACTER(2) :: scale_loc  
      
      INTEGER :: region_box_option
      REAL(rp) :: region_box(4)
      
      
      CHARACTER(100) :: main_font = "/Times-Roman"
      CHARACTER(100) :: math_font = "/Times-Italic"    
!       CHARACTER(100) :: main_font = "(/usr/share/fonts/type1/gsfonts/cmr10.pfb)"
      INTEGER :: fontsize          
      
      INTEGER :: el,nd,pt
      INTEGER :: et,nv

      
      INTEGER :: adapt_option
      INTEGER :: ps
      INTEGER :: pc
      INTEGER :: p_low
      INTEGER :: p_high
      INTEGER :: p_skip
      INTEGER :: nord
      INTEGER, DIMENSION(:), ALLOCATABLE :: pplt
      INTEGER, DIMENSION(:), ALLOCATABLE :: npplt
      INTEGER, DIMENSION(:), ALLOCATABLE :: nptri
      INTEGER :: mnpp
      INTEGER :: nred 
      
      
      INTEGER, DIMENSION(:), ALLOCATABLE :: ndof_sol      
      REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: phi_sol 

      
      INTEGER :: i,j
      INTEGER :: nsnap
      INTEGER :: snap
      INTEGER :: space
      INTEGER :: n,ndf,npts,nnd,dof
      REAL(rp) :: xpt,ypt
      REAL(rp), TARGET :: t_snap
      REAL(rp) :: t_start,t_end
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: t      
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: r,s  
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: rect        
      REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: xyplt  
      REAL(rp), DIMENSION(:), ALLOCATABLE :: el_area
      REAL(rp), DIMENSION(:), ALLOCATABLE :: el_size          
      INTEGER :: outside
      INTEGER, DIMENSION(:), ALLOCATABLE :: el_in
      REAL(rp) :: xbox_min,xbox_max,ybox_min,ybox_max
      REAL(rp) :: figure_width
      REAL(rp) :: figure_height
      REAL(rp) :: line_width      

      
!       TYPE :: viz
!         INTEGER, DIMENSION(:), ALLOCATABLE :: p
!         INTEGER, DIMENSION(:), ALLOCATABLE :: nnds
!         REAL(rp), DIMENSION(:,:), ALLOCATABLE :: r,s
!         INTEGER, DIMENSION(:), ALLOCATABLE :: ntri        
!         INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: tri
!         REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: xy
!         REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: phi_sol 
!         REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: phi_hb
!       END TYPE
!       
!       TYPE(viz) :: plt
      
      INTEGER :: snap_start
      INTEGER :: snap_end
      INTEGER :: sta_start
      INTEGER :: sta_end
      
      INTEGER :: nsnap_Z
      INTEGER :: nsnap_Qx
      INTEGER :: nsnap_Qy
      INTEGER :: nsnap_hb
      
      INTEGER :: tex_unit   
               
      CHARACTER(3) :: start_num
      CHARACTER(3) :: nframes
      
      TYPE :: char_array
        CHARACTER(:),ALLOCATABLE :: line
      END TYPE
            
      INTEGER :: spherical_flag     
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: map
      INTEGER :: map_width
      INTEGER :: map_height
      REAL(rp) :: map_res
      REAL(rp) :: lamc
      REAL(rp) :: phic

      TYPE :: solution_type
        CHARACTER(100) :: input_path      
        CHARACTER(:),ALLOCATABLE :: out_direc  
        CHARACTER(:),ALLOCATABLE :: grid_file 
        CHARACTER(:),ALLOCATABLE :: curve_file 
        CHARACTER(:),ALLOCATABLE :: bathy_file    
        CHARACTER(:),ALLOCATABLE :: stations_file
        CHARACTER(:),ALLOCATABLE :: output_type
        INTEGER :: sta_opt              
      
        INTEGER :: p
        INTEGER :: ctp
        INTEGER :: hbp      
        
        CHARACTER(100) :: grid_name          
        INTEGER :: ne
        INTEGER :: nn       
        INTEGER :: nel_type   
        INTEGER :: ndof(4)
        INTEGER :: mndof         
        INTEGER :: nverts(4)
        INTEGER :: np(4)
        INTEGER :: mnp
        INTEGER :: nnds(4)
        INTEGER :: mnnds

        INTEGER :: nord
        INTEGER :: mnpp
        INTEGER, DIMENSION(:), ALLOCATABLE :: pplt
        INTEGER, DIMENSION(:), ALLOCATABLE :: npplt
        INTEGER, DIMENSION(:), ALLOCATABLE :: nptri
        INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: rect
        
        
        REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: r
        REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: s
        REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: psic   
        REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: psiv         
        
        REAL(rp) :: sphi0
        REAL(rp) :: slam0
        REAL(rp) :: h0
        REAL(rp) :: g = 9.81d0          
        REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: xy
        REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: elxy
        REAL(rp), ALLOCATABLE, DIMENSION(:) :: depth
        REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: elhb
        REAL(rp), ALLOCATABLE, DIMENSION(:,:,:,:) :: bndxy   
        REAL(rp) :: mesh_line_color(3)        
               
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ect
        INTEGER, ALLOCATABLE, DIMENSION(:) :: el_type
        REAL(rp), DIMENSION(:), ALLOCATABLE :: el_area
        REAL(rp), DIMENSION(:), ALLOCATABLE :: el_size   
        INTEGER, DIMENSION(:), ALLOCATABLE :: el_in        
        
        INTEGER :: nope 
        INTEGER, ALLOCATABLE, DIMENSION(:) :: obseg 
        INTEGER :: neta
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: obnds 

        INTEGER :: nbou  
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: fbseg 
        INTEGER :: nvel  
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: fbnds 
      
        INTEGER :: ned       
        INTEGER :: mnepn
        INTEGER, ALLOCATABLE, DIMENSION(:) :: nepn 
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: epn   
     
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ged2nn 
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ged2el 
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ged2led 
        INTEGER, ALLOCATABLE, DIMENSION(:) :: recv_edge    
      
        INTEGER :: nied 
        INTEGER, ALLOCATABLE, DIMENSION(:) :: iedn 
        INTEGER :: nbed
        INTEGER, ALLOCATABLE, DIMENSION(:) :: bedn
        INTEGER :: nobed 
        INTEGER, ALLOCATABLE, DIMENSION(:) :: obedn  
        INTEGER :: nfbed 
        INTEGER, ALLOCATABLE, DIMENSION(:) :: fbedn 
        INTEGER :: nnfbed 
        INTEGER, ALLOCATABLE, DIMENSION(:) :: nfbedn 
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: nfbednn       
        INTEGER, ALLOCATABLE, DIMENSION(:) :: ed_type        
        
        REAL(rp) :: tf
        REAL(rp), ALLOCATABLE, DIMENSION(:) :: t         
        REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: Z
        REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: hb
        REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: Qx
        REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: Qy 
        
!         REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: hbm            
!         REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: eta
!         REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: uu2
!         REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: vv2          
        
        INTEGER :: nsnap_Z
        INTEGER :: nsnap_Qx
        INTEGER :: nsnap_Qy
        INTEGER :: nsnap_hb
      END TYPE
      
      TYPE(solution_type) :: sol1
      TYPE(solution_type) :: sol2
          
      
      TYPE :: plot_type
      
        INTEGER :: plot_sol_option
        INTEGER :: plot_mesh_option
        INTEGER :: plot_lines_option
        INTEGER :: plot_sta_option
        INTEGER :: plot_max_option
        INTEGER :: sol_diff_option
        INTEGER :: ho_diff_option
        INTEGER :: sta_start
        INTEGER :: sta_end
        INTEGER :: plim
        
        INTEGER :: cbar_flag
        INTEGER :: tbar_flag
        INTEGER :: axis_label_flag
        LOGICAL :: tex_file_exists
        
        INTEGER :: rm_ps
        INTEGER :: movie_flag        
        
        CHARACTER(100) :: sol_label
        CHARACTER(100) :: name
        
        INTEGER :: ps_unit
        REAL(rp) :: sol_min,sol_max
        REAL(rp) :: snap_min,snap_max
        INTEGER :: nsnap
        INTEGER, DIMENSION(:), ALLOCATABLE :: el_plt
        REAL(rp), DIMENSION(:,:), ALLOCATABLE :: sol_val
        REAL(rp), DIMENSION(:,:), ALLOCATABLE :: sol_maxval
        REAL(rp) :: max_maxval
        REAL(rp) :: min_maxval
        REAL(rp), DIMENSION(:), ALLOCATABLE :: t
        REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: phi     
        
        REAL(rp) :: rel_tol
        REAL(rp) :: abs_tol
      
        INTEGER :: cscale_unit
        CHARACTER(100) :: cscale_option
        INTEGER :: num_cscale_vals
        REAL(rp), DIMENSION(:,:), ALLOCATABLE :: cscale_vals
        REAL(rp) :: cscale_min,cscale_max
        
        CHARACTER(100) :: el_label_option    
        CHARACTER(100) :: nd_label_option 
        
        INTEGER :: type_flag 
        
        REAL(rp), POINTER :: t_snap
 
        INTEGER :: nline_header
        TYPE(char_array), DIMENSION(:), ALLOCATABLE :: latex_header
        INTEGER :: nline_body
        TYPE(char_array), DIMENSION(:), ALLOCATABLE :: latex_body        
      
      
        REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: sta_val
        REAL(rp), DIMENSION(:,:), ALLOCATABLE :: t_sta        
      END TYPE
      
      TYPE(plot_type) :: zeta
      TYPE(plot_type) :: vel
      TYPE(plot_type) :: bathy
      TYPE(plot_type) :: mesh  
      TYPE(plot_type) :: cfl
      
!       TYPE(plot_type), DIMENSION(:), ALLOCATABLE :: zeta_sta
!       TYPE(plot_type), DIMENSION(:), ALLOCATABLE :: vel_sta
      
      INTEGER :: nadc
      INTEGER :: ndg
      
      TYPE(char_array), DIMENSION(:), ALLOCATABLE :: adc_sol
      TYPE(char_array), DIMENSION(:), ALLOCATABLE :: dg_sol
      
      INTEGER :: xsta_min,xsta_max
      REAL(rp) :: xtime_min,xtime_max
      CHARACTER(4) :: xc_snap_opt
      CHARACTER(4) :: ts_sta_opt
      
      INTEGER :: plot_ts_sta_opt
      INTEGER :: plot_xc_snap_opt
      INTEGER :: plot_sta_loc_opt
      INTEGER :: plot_error_opt
      INTEGER :: plot_scatter_opt

      
      CONTAINS
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      SUBROUTINE setup_plot_types()
      
      IMPLICIT NONE           
      
      mesh%type_flag = 1     
      bathy%type_flag = 2    
      zeta%type_flag = 3     
      vel%type_flag = 4
      cfl%type_flag = 5
      
      mesh%axis_label_flag = 1
      mesh%cbar_flag = 0
      mesh%tbar_flag = 0
      bathy%axis_label_flag = 1
      bathy%cbar_flag = 1
      bathy%tbar_flag = 0
      zeta%axis_label_flag = 1
      zeta%cbar_flag = 1
      zeta%tbar_flag = 1
      vel%axis_label_flag = 1
      vel%cbar_flag = 1
      vel%tbar_flag = 1
      cfl%axis_label_flag = 1
      cfl%cbar_flag = 1
      cfl%tbar_flag = 0
      
      IF (nxdec < 0 .or. nydec < 0) THEN
        mesh%axis_label_flag = 0
        bathy%axis_label_flag = 0
        zeta%axis_label_flag = 0
        vel%axis_label_flag = 0   
        cfl%axis_label_flag = 0        
      ENDIF
      
      IF (ncdec < 0) THEN
        bathy%cbar_flag = 0
        zeta%cbar_flag = 0
        vel%cbar_flag = 0
        cfl%cbar_flag = 0
        ncdec = abs(ncdec)
      ENDIF
      
      IF (ntdec < 0) THEN
        zeta%tbar_flag = 0
        vel%tbar_flag = 0
      ENDIF       
      
      mesh%sol_diff_option = 0
      bathy%sol_diff_option = 0
      zeta%sol_diff_option = 0
      vel%sol_diff_option = 0
      cfl%sol_diff_option = 0
      IF (sol_diff_option == 1) THEN
        bathy%sol_diff_option = 1
        zeta%sol_diff_option = 1
        vel%sol_diff_option = 1
        mesh%sol_diff_option = 1
      ENDIF
      
      
      mesh%ho_diff_option = 0
      bathy%ho_diff_option = 0 
      zeta%ho_diff_option = 0
      vel%ho_diff_option = 0
      cfl%ho_diff_option = 0
      IF (ho_diff_option == 1) THEN
        bathy%ho_diff_option = 1
        zeta%ho_diff_option = 1
        vel%ho_diff_option = 1
        
        bathy%plim = 1
        zeta%plim = 1
        vel%plim = 1
      ENDIF
      
      bathy%el_label_option = "off"
      bathy%nd_label_option = "off"
      zeta%el_label_option = "off"
      zeta%nd_label_option = "off"
      vel%el_label_option = "off"
      vel%nd_label_option = "off"
      cfl%el_label_option = "off"
      cfl%nd_label_option = "off"
      
      mesh%sol_label = "mesh"
      mesh%name = "mesh"
      zeta%sol_label = "surface elevation (m)"
      zeta%name = "zeta"      
      vel%sol_label = "velocity (m/s)"
      vel%name = "vel"      
      bathy%sol_label = "bathymetry (m)"
      bathy%name = "bathy"       
      cfl%sol_label = "Courant number"
      cfl%name = "cfl"      
      
      mesh%t_snap => t_snap
      bathy%t_snap => t_snap
      zeta%t_snap => t_snap
      vel%t_snap => t_snap 
      cfl%t_snap => t_snap
      
      zeta%cscale_unit = 30
      bathy%cscale_unit = 31
      vel%cscale_unit = 32      
      
      cfl%cscale_option = "auto-snap"
!       cfl%cscale_option = "spec"      
!       cfl%cscale_min = 0.4d0
!       cfl%cscale_max = 7d0

      sol1%mesh_line_color = (/  0d0,  0d0,  0d0 /)
!       sol1%mesh_line_color = (/  .35d0,  .35d0,  .35d0 /)      
      sol2%mesh_line_color = (/ .5d0, .5d0, .5d0 /)
      
      END SUBROUTINE setup_plot_types

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      END MODULE plot_globals
