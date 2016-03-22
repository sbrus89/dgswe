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
      INTEGER, PARAMETER :: norder = 6 ! # of different orders (straight sided elements for tri/quad = 1, curvilinear tri/quad = ctp, high-order bathymetry = hbp)       
      INTEGER :: order(2*nel_type)      
      INTEGER :: ctp
      INTEGER :: hbp
      INTEGER :: lsp
      INTEGER :: nverts(nel_type)
      INTEGER :: np(norder)
      INTEGER :: nnds(norder)
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
      
        INTEGER :: mnepn
        INTEGER, ALLOCATABLE, DIMENSION(:) :: nepn ! number of elements per node
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: epn ! elements per node            
      
        INTEGER :: ned
        INTEGER :: nied
        INTEGER, ALLOCATABLE, DIMENSION(:) :: iedn  
        INTEGER :: nobed    
        INTEGER, ALLOCATABLE, DIMENSION(:) :: obedn
        INTEGER :: nfbed
        INTEGER, ALLOCATABLE, DIMENSION(:) :: fbedn        
        INTEGER :: nnfbed
        INTEGER, ALLOCATABLE, DIMENSION(:) :: nfbedn 
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: nfbednn         
        INTEGER :: nbed
        INTEGER, ALLOCATABLE, DIMENSION(:) :: bedn   
        INTEGER :: nred
        INTEGER, ALLOCATABLE, DIMENSION(:) :: recv_edge
        INTEGER, ALLOCATABLE, DIMENSION(:) :: ed_type
        
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ged2nn,ged2el,ged2led  
        
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: el2el
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: el2ged            

               
      
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
      REAL(rp) :: eps
      
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
