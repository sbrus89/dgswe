      MODULE globals
      
      USE kdtree2_module
      
      IMPLICIT NONE

      INTEGER, PARAMETER :: rp = kind(1d0) ! precision   
      
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

      INTEGER :: lines ! number of lines in output files  
      REAL(rp) :: tf ! final time   
      
      LOGICAL :: exclude_bndel
      
      TYPE :: grid 
        CHARACTER(10) :: sol_name
        CHARACTER(100) :: grid_file ! name of fort.14 file
        CHARACTER(100) :: grid_name ! name of the grid
        CHARACTER(100) :: out_direc ! name of output directory
        CHARACTER(100) :: bathy_file         
        INTEGER :: hbp
        INTEGER :: ctp        
        INTEGER :: ne ! number of elements
        INTEGER :: nn ! number of nodes
        INTEGER :: curved_grid
        
        INTEGER, ALLOCATABLE, DIMENSION(:) :: el_type   
        INTEGER :: nverts(nel_type)
        INTEGER :: np(nel_type+2)
        INTEGER :: nnds(nel_type+2)
        INTEGER :: mnnds             
      
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ect,vct ! element connectivity table
        INTEGER, ALLOCATABLE, DIMENSION(:) :: nelnds
        INTEGER :: mnelnds      
        REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: xy,vxy ! x,y coordinates of nodes
        INTEGER, ALLOCATABLE, DIMENSION(:) :: vxyn
        REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: elxy    
        REAL(rp), ALLOCATABLE, DIMENSION(:) :: depth ! depth at each node
        REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: elhb   
        REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: hbxy 
        INTEGER :: npts
 
        
        INTEGER :: nope ! number of open boundary segents
        INTEGER, ALLOCATABLE, DIMENSION(:) :: obseg ! number of nodes in each open boundary segment
        INTEGER :: neta ! total elevation specified boundary nodes
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: obnds ! open boundary nodes

        INTEGER :: nbou  ! number of normal flow boundary segments
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: fbseg ! number of nodes and type of each normal flow boundary segment
        INTEGER :: nvel  ! total number of normal flow boundary nodes
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: fbnds ! normal flow boundary nodes   
        
        INTEGER, ALLOCATABLE, DIMENSION(:) :: bndel ! array to flag land boundary elements
        
        INTEGER :: mnepn
        INTEGER, ALLOCATABLE, DIMENSION(:) :: nepn ! number of elements per node
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: epn ! elements associated with each node 
        
        INTEGER, ALLOCATABLE, DIMENSION(:) :: nepe
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: el2el
      
        INTEGER :: ned
        INTEGER :: nied
        INTEGER :: nnfbed
        INTEGER, ALLOCATABLE, DIMENSION(:) :: bed_flag
        INTEGER, ALLOCATABLE, DIMENSION(:) :: bnd_flag        
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: bel2bed
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ged2nn,ged2el,ged2led  
        INTEGER, ALLOCATABLE, DIMENSION(:) :: nfbedn              
        
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
