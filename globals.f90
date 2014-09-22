      MODULE globals
      IMPLICIT NONE

      INTEGER, PARAMETER :: pres = kind(1d0) ! precision 
      INTEGER :: p ! polynomial order      
      INTEGER :: ne ! number of elements
      INTEGER :: nn ! number of nodes
      INTEGER :: npart ! number of element partitions      
      REAL(pres) :: cf ! bottom friction parameter  
      REAL(pres) :: dt ! time step 
      REAL(pres) :: tf ! final time      
      REAL(pres) :: dramp ! numer of ramp days      
      REAL(pres) :: lines ! number of lines in output files      
      CHARACTER(100) :: grid_name ! name of the grid
      CHARACTER(100) :: grid_file ! name of fort.14 file
      CHARACTER(100) :: forcing_file ! name of fort.15 file      
      CHARACTER(100) :: out_direc           
      
      INTEGER, PARAMETER :: nel_type = 4 !(type #s: 1 -> triangles, 2 -> quads, 3 -> curved triangles, 4-> curved quads)   
      INTEGER :: ctp
      INTEGER :: nverts(nel_type)
      INTEGER :: np(nel_type)
      INTEGER :: nnds(nel_type)
      INTEGER :: mnnds      
      INTEGER :: ndof(nel_type)
      INTEGER :: mndof
      
      
      INTEGER :: curved_grid
      INTEGER, ALLOCATABLE, DIMENSION(:) :: el_type      
      
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ect,vct ! element connectivity table
      INTEGER, ALLOCATABLE, DIMENSION(:) :: nelnds
      INTEGER :: mnelnds      
      REAL(pres), ALLOCATABLE, DIMENSION(:,:) :: xy ! x,y coordinates of nodes
      REAL(pres), ALLOCATABLE, DIMENSION(:,:,:) :: elxy    
      REAL(pres), ALLOCATABLE, DIMENSION(:) :: depth ! depth at each node
      REAL(pres), ALLOCATABLE, DIMENSION(:,:) :: elhb     
      
      INTEGER :: nsta
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ndsta
      INTEGER, ALLOCATABLE, DIMENSION(:) :: elsta
      REAL(pres), ALLOCATABLE, DIMENSION(:,:) :: xysta ! x,y coordinates of stations
      REAL(pres), ALLOCATABLE, DIMENSION(:,:) :: rssta

      INTEGER :: nope ! number of open boundary segents
      INTEGER, ALLOCATABLE, DIMENSION(:) :: obseg ! number of nodes in each open boundary segment
      INTEGER :: neta ! total elevation specified boundary nodes
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: obnds ! open boundary nodes

      INTEGER :: nbou  ! number of normal flow boundary segments
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: fbseg ! number of nodes and type of each normal flow boundary segment
      INTEGER :: nvel  ! total number of normal flow boundary nodes
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: fbnds ! normal flow boundary nodes
      
      INTEGER, ALLOCATABLE, DIMENSION(:) :: nepn ! number of elements per node
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: epn ! elements per node            
   
      REAL(pres), ALLOCATABLE, DIMENSION(:,:) :: H
      REAL(pres), ALLOCATABLE, DIMENSION(:,:) :: Qx
      REAL(pres), ALLOCATABLE, DIMENSION(:,:) :: Qy     
      
      REAL(pres), ALLOCATABLE, DIMENSION(:,:,:) :: V
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ipiv

      END MODULE globals
