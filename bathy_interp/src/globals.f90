      MODULE globals
      
      USE kdtree2_module
      
      IMPLICIT NONE

      INTEGER, PARAMETER :: rp = kind(1d0) ! precision   
      REAL(rp), PARAMETER  ::  pi=3.141592653589793D0     
      REAL(rp), PARAMETER :: deg2rad = pi/180d0      
      REAL(rp), PARAMETER :: g = 9.81d0        
      
      CHARACTER(50) :: gitSHA
      CHARACTER(50) :: gitBranch
      CHARACTER(100) :: compiler_version
      CHARACTER(50) :: compiler_flags
      CHARACTER(1000) :: modified_files
      CHARACTER(50) :: compile_date
      CHARACTER(50) :: host       
      
      INTEGER, PARAMETER :: srchdp = 10
      TYPE(kdtree2), POINTER :: tree_xy
      TYPE(kdtree2_result), ALLOCATABLE, DIMENSION(:) :: closest         
      
      INTEGER, PARAMETER :: nel_type = 4 
      INTEGER :: nverts(nel_type)      

      INTEGER :: lines ! number of lines in output files  
      REAL(rp) :: tf ! final time   
      
      LOGICAL :: exclude_bndel
      
      TYPE :: grid 
        CHARACTER(10) :: sol_name
        CHARACTER(100) :: grid_file ! name of fort.14 file
        CHARACTER(100) :: grid_name ! name of the grid
        CHARACTER(100) :: out_direc ! name of output directory
        CHARACTER(100) :: bathy_file 
        CHARACTER(100) :: curve_file
        INTEGER :: hbp
        INTEGER :: ctp        
        INTEGER :: ne ! number of elements
        INTEGER :: nn ! number of nodes
        
        INTEGER, ALLOCATABLE, DIMENSION(:) :: el_type   
        INTEGER :: np(nel_type+2)
        INTEGER :: nnds(nel_type+2)
        INTEGER :: mnnds             
      
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ect ! element connectivity table     
        REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: xy ! x,y coordinates of nodes
        REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: elxy    
        REAL(rp), ALLOCATABLE, DIMENSION(:) :: depth ! depth at each node
        REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: elhb   
        REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: elhbxy        
        REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: hbxy 
        REAL(rp), ALLOCATABLE, DIMENSION(:,:,:,:) :: bndxy         
        INTEGER :: npts
 
        
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
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: epn ! elements associated with each node 
        
      
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
        
        INTEGER, ALLOCATABLE, DIMENSION(:) :: bnd_flag
        
        REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: V ! vandermonde matrix
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ipiv 
        
        REAL(rp) :: rsre(2,4,nel_type) ! reference element verticies 
        
        REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: l,dldr,dlds   
                 
        
      END TYPE
      
      TYPE(grid) :: eval    
      TYPE(grid) :: base

      INTEGER :: mnepn
      
      INTEGER :: nept(nel_type)
      INTEGER :: mnept
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: ept       
        


      

      


      END MODULE globals
