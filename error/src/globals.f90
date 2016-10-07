      MODULE globals
      IMPLICIT NONE

      INTEGER, PARAMETER :: rp = kind(1d0) ! precision   
      
      CHARACTER(50) :: gitSHA
      CHARACTER(50) :: gitBranch
      CHARACTER(100) :: compiler_version
      CHARACTER(50) :: compiler_flags
      CHARACTER(1000) :: modified_files
      CHARACTER(50) :: compile_date
      CHARACTER(50) :: host       
      
      REAL(rp), PARAMETER :: pi = 3.141592653589793D0
      REAL(rp), PARAMETER :: deg2rad = pi/180d0  
      REAL(rp), PARAMETER :: g = 9.81d0      
      
      INTEGER, PARAMETER :: nel_type = 4 !(type #s: 1 -> triangles, 2 -> quads, 3 -> curved triangles, 4-> curved quads)         
      INTEGER :: nverts(nel_type)   

      INTEGER :: lines ! number of lines in output files  
      REAL(rp) :: tf ! final time   
      
      LOGICAL :: exclude_bndel
         
      
      TYPE :: solution 
        CHARACTER(10) :: sol_name
        CHARACTER(100) :: grid_file ! name of fort.14 file
        CHARACTER(100) :: grid_name ! name of the grid
        CHARACTER(100) :: out_direc ! name of output directory
        CHARACTER(100) :: curve_file        
        INTEGER :: p ! polynomial order
        INTEGER :: ctp 
        REAL(rp) :: dt ! time step        
        INTEGER :: ne ! number of elements
        INTEGER :: nn ! number of nodes
        INTEGER :: curved_grid
        
        INTEGER, ALLOCATABLE, DIMENSION(:) :: el_type   
        INTEGER :: np(nel_type)
        INTEGER :: nnds(nel_type)
        INTEGER :: mnnds      
        INTEGER :: ndof(nel_type)
        INTEGER :: mndof         
      
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ect ! element connectivity table    
        REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: xy ! x,y coordinates of nodes
        REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: elxy    
        REAL(rp), ALLOCATABLE, DIMENSION(:) :: depth ! depth at each node
        REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: elhb   
        REAL(rp), ALLOCATABLE, DIMENSION(:,:,:,:) :: bndxy        
        
        INTEGER :: nope ! number of open boundary segents
        INTEGER, ALLOCATABLE, DIMENSION(:) :: obseg ! number of nodes in each open boundary segment
        INTEGER :: neta ! total elevation specified boundary nodes
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: obnds ! open boundary nodes

        INTEGER :: nbou  ! number of normal flow boundary segments
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: fbseg ! number of nodes and type of each normal flow boundary segment
        INTEGER :: nvel  ! total number of normal flow boundary nodes
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: fbnds ! normal flow boundary nodes   
        
        INTEGER, ALLOCATABLE, DIMENSION(:) :: bndel ! array to flag land boundary elements
        
        INTEGER, ALLOCATABLE, DIMENSION(:) :: nepn ! number of elements per node
        INTEGER :: mnepn
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
        
        
        REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: V ! vandermonde matrix
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ipiv 
        
        REAL(rp) :: rsre(2,4,nel_type) ! reference element verticies 
        
        REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: l,dldr,dlds
        REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: phi
        REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: psiv       
        REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: detJ      
        
        REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: H
        REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: Qx
        REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: Qy             
        
      END TYPE
      
      TYPE(solution) :: coarse
      TYPE(solution) :: fine    
      TYPE(solution) :: base
      
      INTEGER :: nqpta(nel_type)
      INTEGER :: mnqpta
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: wpta
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: qpta      
      
      INTEGER, ALLOCATABLE, DIMENSION(:) :: elf2elc
      INTEGER, ALLOCATABLE, DIMENSION(:) :: elf2elb
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: xysta ! x,y coordinates of stations
      
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: r,s,hb      
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: xf,yf      
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: Hc,Qxc,Qyc
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: Hf,Qxf,Qyf            
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: phi      


      

      


      END MODULE globals
