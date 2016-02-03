      MODULE globals
      
      USE kdtree2_module      
      
      IMPLICIT NONE

      INTEGER, PARAMETER :: rp = kind(1d0) ! precision 
      REAL(rp), PARAMETER  ::  pi=3.141592653589793D0     
   
      CHARACTER(100) :: out_direc           
      
      INTEGER, PARAMETER :: nel_type = 4 !(type #s: 1 -> triangles, 2 -> quads, 3 -> curved triangles, 4-> curved quads)   
      INTEGER :: ctp
      INTEGER :: nverts(nel_type)
      INTEGER :: np(nel_type)
      INTEGER :: nnds(nel_type)
      INTEGER :: mnnds         
      INTEGER :: mninds  
      
      REAL(rp) :: Erad
      REAL(rp) :: lambda0
      REAL(rp) :: phi0
                          
                          
      
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
        INTEGER, ALLOCATABLE, DIMENSION(:) :: bel_flag        
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ged2nn,ged2el,ged2led  
        INTEGER, ALLOCATABLE, DIMENSION(:) :: nfbedn        
        
        REAL(rp), ALLOCATABLE, DIMENSION(:) :: minedlen,edlen
      
      END TYPE      
      
      INTEGER :: nfbnds
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: fbnds_xy
      INTEGER , ALLOCATABLE, DIMENSION(:) :: fbnds
      
      
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: ax,bx,cx,dx
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: ay,by,cy,dy      
 
   
      TYPE(grid) :: base
      TYPE(grid) :: eval   
      
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: V
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: l
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: dldr
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: dlds      
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ipiv      
      
      INTEGER, PARAMETER :: srchdp = 20
      REAL(rp) :: rsre(2,4,4)        
      TYPE(kdtree2), POINTER :: tree_xy      
      TYPE(kdtree2_result), ALLOCATABLE, DIMENSION(:) :: closest       


      END MODULE globals
