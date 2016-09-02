      MODULE globals    
      
      IMPLICIT NONE

      SAVE
      
      INTEGER, PARAMETER :: rp = kind(1d0) ! precision 
      REAL(rp), PARAMETER  ::  pi=3.141592653589793D0     
      REAL(rp), PARAMETER :: deg2rad = pi/180d0      
      REAL(rp), PARAMETER :: g = 9.81d0            
   
      CHARACTER(100) :: out_direc           
      
      INTEGER, PARAMETER :: nel_type = 4 !(type #s: 1 -> triangles, 2 -> quads, 3 -> curved triangles, 4-> curved quads)   
      INTEGER :: nverts(nel_type)
      
   
      
      
      REAL(rp) :: Erad
      REAL(rp) :: lambda0
      REAL(rp) :: phi0
      
      REAL(rp) :: theta_tol
      REAL(rp) :: deform_tol
      REAL(rp) :: sig
                          
                          
      
      TYPE :: grid
      
        CHARACTER(100) :: grid_name ! name of the grid
        CHARACTER(100) :: grid_file ! name of fort.14 file 
        
        INTEGER :: ne ! number of elements
        INTEGER :: nn ! number of nodes        
        
        INTEGER :: ctp
        INTEGER :: np(nel_type)
        INTEGER :: nnds(nel_type)
        INTEGER :: mnnds           
        REAL(rp), ALLOCATABLE, DIMENSION(:) :: rpts  
        
        INTEGER :: nqpta(nel_type),mnqpta
        INTEGER :: nqpte(nel_type),mnqpte
        REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: wpta,wpte
        REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: qpta,qpte
        REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: psi,dpdr,dpds           
        REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: psiv
        
        INTEGER, ALLOCATABLE, DIMENSION(:) :: el_type,el_type_spline      
      
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ect ! element connectivity table     
        REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: xy ! x,y coordinates of nodes
        REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: elxy,elxy_spline        
        REAL(rp), ALLOCATABLE, DIMENSION(:) :: depth ! depth at each node      
        REAL(rp), ALLOCATABLE, DIMENSION(:) :: h
        REAL(rp), ALLOCATABLE, DIMENSION(:,:,:,:) :: bndxy
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: bnd_flag 

        INTEGER :: nope ! number of open boundaries
        INTEGER, ALLOCATABLE, DIMENSION(:) :: obseg ! number of nodes in each open boundary
        INTEGER :: neta ! total elevation specified boundary nodes
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: obnds ! open boundary nodes

        INTEGER :: nbou  ! number of normal flow boundaries
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: fbseg ! number of nodes and type of each normal flow boundary
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
        REAL(rp), ALLOCATABLE, DIMENSION(:) :: minedlen,edlen
        
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: el2el
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: el2ged          
      
      END TYPE      
      
      INTEGER :: nfbnds
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: fbnds_xy
      INTEGER , ALLOCATABLE, DIMENSION(:) :: fbnds
      INTEGER, ALLOCATABLE, DIMENSION(:) :: nfbnd2el
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: fbnd2el
      

      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: ax,bx,cx,dx
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: ay,by,cy,dy      
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: dt 

      

 
   
      TYPE(grid) :: base
      TYPE(grid) :: eval   
      
!       REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: V
!       REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: l
!       REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: dldr
!       REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: dlds      
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ipiv                
      
      CHARACTER(50) :: gitSHA
      CHARACTER(50) :: gitBranch
      CHARACTER(100) :: compiler_version
      CHARACTER(50) :: compiler_flags
      CHARACTER(1000) :: modified_files
      CHARACTER(50) :: compile_date
      CHARACTER(50) :: host           


      END MODULE globals
