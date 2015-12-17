      MODULE globals
      
      USE kdtree2_module
      
      IMPLICIT NONE

      INTEGER, PARAMETER :: rp = kind(1d0) ! precision 
      INTEGER :: p ! polynomial order      
      INTEGER :: npart ! number of element partitions      
      REAL(rp) :: cf ! bottom friction parameter  
      REAL(rp) :: dt ! time step 
      REAL(rp) :: tf ! final time      
      REAL(rp) :: dramp ! numer of ramp days      
      REAL(rp) :: lines ! number of lines in output files      

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
      INTEGER :: mninds
      
      LOGICAL :: refinement     
      
      TYPE :: grid
      
        CHARACTER(100) :: grid_name ! name of the grid
        CHARACTER(100) :: grid_file ! name of fort.14 file 
        
        INTEGER :: ne ! number of elements
        INTEGER :: nn ! number of nodes        
      
        INTEGER :: curved_grid
        INTEGER, ALLOCATABLE, DIMENSION(:) :: el_type      
      
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ect,vct ! element connectivity table
        INTEGER, ALLOCATABLE, DIMENSION(:) :: nelnds
        INTEGER :: mnelnds      
        REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: xy,vxy ! x,y coordinates of nodes
        INTEGER, ALLOCATABLE, DIMENSION(:) :: vxyn
        REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: elxy        
        REAL(rp), ALLOCATABLE, DIMENSION(:) :: depth ! depth at each node
        REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: elhb     
        REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: nhb
        REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: xyhc      
        REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: xyhe      
        REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: xyhi 
        REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: xyhv       
        REAL(rp), ALLOCATABLE, DIMENSION(:) :: h
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: bnd_flag 

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
        
        INTEGER, ALLOCATABLE, DIMENSION(:) :: nepe
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: el2el
      
        INTEGER :: ned
        INTEGER :: nied
        INTEGER :: nnfbed
        INTEGER, ALLOCATABLE, DIMENSION(:) :: bed_flag
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: bel2bed
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ged2nn,ged2el,ged2led  
        INTEGER, ALLOCATABLE, DIMENSION(:) :: nfbedn        
      
      END TYPE
      
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: xyhw          
      
      TYPE(grid) :: base
      TYPE(grid) :: eval
      
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: V
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: l
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: dldr
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: dlds      
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ipiv
      
      REAL(rp) :: r,sigma_n
      REAL(rp) :: hmin,hmax
      
      INTEGER :: srchdp
      REAL(rp) :: rsre(2,4,4)      
      TYPE(kdtree2), POINTER :: tree_xy,tree_c,tree_xy_rand
      TYPE(kdtree2_result), ALLOCATABLE, DIMENSION(:) :: kdresults  
      TYPE(kdtree2_result), ALLOCATABLE, DIMENSION(:) :: closest    
      
      INTEGER :: nrpt
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: xy_rand     
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: h_rand
      
      REAL(rp) :: Erad 
      REAL(rp) :: phi0,lambda0
      REAL(rp), PARAMETER :: pi = 3.141592653589793D0
      REAL(rp), PARAMETER :: deg2rad = pi/180d0
      
      CHARACTER(50) :: gitSHA
      CHARACTER(50) :: gitBranch
      CHARACTER(100) :: compiler_version
      CHARACTER(50) :: compiler_flags
      CHARACTER(1000) :: modified_files
      CHARACTER(50) :: compile_date
      CHARACTER(50) :: host      

      END MODULE globals
