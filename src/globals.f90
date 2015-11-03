      MODULE globals
      
      USE kdtree2_module
      
      IMPLICIT NONE

      INTEGER, PARAMETER :: pres = kind(1d0) ! precision 
      INTEGER :: p ! polynomial order      
      INTEGER :: npart ! number of element partitions      
      REAL(pres) :: cf ! bottom friction parameter  
      REAL(pres) :: dt ! time step 
      REAL(pres) :: tf ! final time      
      REAL(pres) :: dramp ! numer of ramp days      
      REAL(pres) :: lines ! number of lines in output files      

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
        REAL(pres), ALLOCATABLE, DIMENSION(:,:) :: xy,vxy ! x,y coordinates of nodes
        INTEGER, ALLOCATABLE, DIMENSION(:) :: vxyn
        REAL(pres), ALLOCATABLE, DIMENSION(:,:,:) :: elxy        
        REAL(pres), ALLOCATABLE, DIMENSION(:) :: depth ! depth at each node
        REAL(pres), ALLOCATABLE, DIMENSION(:,:) :: elhb     
        REAL(pres), ALLOCATABLE, DIMENSION(:,:) :: nhb
        REAL(pres), ALLOCATABLE, DIMENSION(:,:) :: xyhc      
        REAL(pres), ALLOCATABLE, DIMENSION(:,:,:) :: xyhe      
        REAL(pres), ALLOCATABLE, DIMENSION(:,:,:) :: xyhi 
        REAL(pres), ALLOCATABLE, DIMENSION(:,:,:) :: xyhv       
        REAL(pres), ALLOCATABLE, DIMENSION(:) :: h
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
      
      REAL(pres), ALLOCATABLE, DIMENSION(:,:,:) :: xyhw          
      
      TYPE(grid) :: base
      TYPE(grid) :: eval
      
      REAL(pres), ALLOCATABLE, DIMENSION(:,:,:) :: V
      REAL(pres), ALLOCATABLE, DIMENSION(:,:,:) :: l
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ipiv
      
      REAL(pres) :: r,sigma_n
      REAL(pres) :: hmin,hmax
      
      INTEGER :: srchdp
      REAL(pres) :: rsre(2,4,4)      
      TYPE(kdtree2), POINTER :: tree_xy,tree_c
      TYPE(kdtree2_result), ALLOCATABLE, DIMENSION(:) :: kdresults  
      TYPE(kdtree2_result), ALLOCATABLE, DIMENSION(:) :: closest       
      
      REAL(pres) :: Erad 
      REAL(pres) :: phi0,lambda0
      REAL(pres), PARAMETER :: pi = 3.141592653589793D0
      REAL(pres), PARAMETER :: deg2rad = pi/180d0

      END MODULE globals
