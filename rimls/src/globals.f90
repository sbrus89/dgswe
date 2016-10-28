      MODULE globals
      
      USE kdtree2_module
      
      IMPLICIT NONE

      INTEGER, PARAMETER :: rp = kind(1d0) ! precision 
    
      CHARACTER(100) :: out_direc           
      
      INTEGER, PARAMETER :: nel_type = 4 !(type #s: 1 -> triangles, 2 -> quads, 3 -> curved triangles, 4-> curved quads)   
      INTEGER, PARAMETER :: norder = 6 ! # of different orders (straight sided elements for tri/quad = 1, curvilinear tri/quad = ctp, high-order bathymetry = hbp)             
      INTEGER :: lsp
      INTEGER :: basis_opt
      INTEGER :: nverts(nel_type)    
      
      LOGICAL :: refinement     
      
      TYPE :: grid
      
        CHARACTER(100) :: grid_name ! name of the grid
        CHARACTER(100) :: grid_file ! name of fort.14 file 
        CHARACTER(100) :: curve_file          
        
        INTEGER :: ne ! number of elements
        INTEGER :: nn ! number of nodes        
        
        INTEGER :: ctp
        INTEGER :: hbp        
        INTEGER :: np(norder)
        INTEGER :: nnds(norder)
        INTEGER :: mnnds         
        INTEGER :: mninds
      
        INTEGER, ALLOCATABLE, DIMENSION(:) :: el_type      
      
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ect ! element connectivity table
        REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: xy ! x,y coordinates of nodes
        REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: elxy
        REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: elxyh           
        REAL(rp), ALLOCATABLE, DIMENSION(:) :: depth ! depth at each node
        REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: elhb 
        REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: dhbdx
        REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: dhbdy        
        REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: nhb  
        
        INTEGER :: tpts_edge
        INTEGER, ALLOCATABLE, DIMENSION(:) :: npts_edge
        REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: xyh_edge  
        INTEGER :: tpts_interior
        INTEGER, ALLOCATABLE, DIMENSION(:) :: npts_interior       
        REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: xyh_interior 
        INTEGER :: tpts_vertex
        INTEGER, ALLOCATABLE, DIMENSION(:) :: npts_vertex        
        REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: xyh_vertex          
        
        REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: elxy_center
        
        REAL(rp), ALLOCATABLE, DIMENSION(:) :: h
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: bnd_flag 
        REAL(rp), ALLOCATABLE, DIMENSION(:,:,:,:) :: bndxy          

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
        INTEGER :: nred
        INTEGER, ALLOCATABLE, DIMENSION(:) :: recv_edge
        INTEGER, ALLOCATABLE, DIMENSION(:) :: ed_type      
        
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ged2nn,ged2el,ged2led  
        
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: el2el
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: el2ged            

               
      
      END TYPE
         
      
      TYPE(grid) :: base
      TYPE(grid) :: eval

      
      REAL(rp) :: r,sigma_n
      REAL(rp) :: hmin,hmax      
      
      INTEGER :: nrpt
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: xy_rand     
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: h_rand
      REAL(rp) :: eps
      
      REAL(rp) :: Erad 
      REAL(rp) :: phi0,lambda0
      REAL(rp), PARAMETER :: g = 9.81d0      
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
