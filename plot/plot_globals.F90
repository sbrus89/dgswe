      MODULE plot_globals
      
      USE globals, ONLY: rp
      
      CHARACTER(100) :: input_path
      CHARACTER(100) :: cmap_file
      CHARACTER(100) :: zbox
      CHARACTER(100) :: tsnap_spec

      INTEGER :: coord_sys
      REAL(rp) :: slam0,sphi0      
      REAL(rp) :: xmin,xmax
      REAL(rp) :: ymin,ymax
      REAL(rp) :: rmin,rmax
      REAL(rp) :: smin,smax
      
      INTEGER :: el,nd,pt
      INTEGER :: et,nv
      REAL(rp) :: ax,bx
      REAL(rp) :: ay,by
      
      INTEGER :: ps
      INTEGER :: pc
      INTEGER :: pplt(4)
      INTEGER :: nplt(4)
      INTEGER :: mnpp
      INTEGER :: ntri(4)
      INTEGER :: nred 
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: rect
      INTEGER :: ndof_hb(4)
      
      INTEGER :: i,j
      INTEGER :: nsnap
      INTEGER :: snap
      INTEGER :: space
      INTEGER :: n,ndf,npts,nnd,dof
      REAL(rp) :: xpt,ypt
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: t      
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: Z,Qx,Qy,hb 
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: r,s
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: phi 
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: phi_hb        
      REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: xyplt   
      REAL(rp), DIMENSION(:,:), ALLOCATABLE :: Z_val 
      REAL(rp), DIMENSION(:,:), ALLOCATABLE :: hb_val           
      REAL(rp), DIMENSION(:,:), ALLOCATABLE :: vel_val   
      REAL(rp) :: Qx_val,Qy_val,H_val
      INTEGER :: outside
      INTEGER, DIMENSION(:), ALLOCATABLE :: el_in
      REAL(rp) :: xbox_min,xbox_max,ybox_min,ybox_max
      REAL(rp) :: figure_width
      
      INTEGER :: plot_mesh_option
      INTEGER :: plot_vel_option
      INTEGER :: plot_zeta_option
      INTEGER :: plot_bathy_option
      
      INTEGER :: snap_start
      INTEGER :: snap_end
      
      INTEGER :: nsnap_Z
      INTEGER :: nsnap_Qx
      INTEGER :: nsnap_Qy
      INTEGER :: nsnap_hb
      
      INTEGER :: Z_unit
      INTEGER :: hb_unit
      INTEGER :: vel_unit
      INTEGER :: mesh_unit
      
      CHARACTER(4) :: snap_char
      
      CONTAINS
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         
      
      SUBROUTINE sizes()
      
      USE globals, ONLY: nverts,ndof,mndof,np,mnp,nnds,mnnds
      USE read_dginp, ONLY: p,ctp,hbp
      
      IMPLICIT NONE
            
      ndof(1) = (p+1)*(p+2)/2
      ndof(2) = (p+1)**2
      ndof(3) = ndof(1)
      ndof(4) = ndof(2)      
      mndof = maxval(ndof)
      
      ndof_hb(1) = (hbp+1)*(hbp+2)/2
      ndof_hb(2) = (hbp+1)*(hbp+1)
      ndof_hb(3) = (hbp+1)*(hbp+2)/2
      ndof_hb(4) = (hbp+1)*(hbp+1)      
      
      nverts(1) = 3
      nverts(2) = 4
      nverts(3) = 3
      nverts(4) = 4
      
      np(1) = 1
      np(2) = 1
      np(3) = ctp
      np(4) = ctp  
      mnp = maxval(np)+1

      nnds(1) = 3
      nnds(2) = 4
      nnds(3) = (ctp+1)*(ctp+2)/2
      nnds(4) = (ctp+1)*(ctp+1)      
      mnnds = maxval(nnds)   
      
      pplt(1) = ps
      pplt(2) = ps
      pplt(3) = pc
      pplt(4) = pc

      nplt(1) = (ps+1)*(ps+2)/2
      nplt(2) = (ps+1)*(ps+1)
      nplt(3) = (pc+1)*(pc+2)/2
      nplt(4) = (pc+1)*(pc+1)
      mnpp = maxval(nplt)      
      
      END SUBROUTINE sizes
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      
      END MODULE plot_globals