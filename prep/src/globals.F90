      MODULE globals

      IMPLICIT NONE

      SAVE
      
      INTEGER, PARAMETER :: rp = kind(1d0)

      CHARACTER(100) :: grid_name
      
      INTEGER, PARAMETER :: nel_type = 4

      INTEGER :: ne
      INTEGER :: nn      
      
      INTEGER, DIMENSION(nel_type) :: nverts
      INTEGER, DIMENSION(nel_type) :: np
      INTEGER :: mnp
      INTEGER, DIMENSION(nel_type) :: nnds
      INTEGER :: mnnds
      INTEGER, DIMENSION(nel_type) :: ndof
      INTEGER :: mndof        
      INTEGER, DIMENSION(nel_type) :: hbnds
      INTEGER :: mhbnds
      INTEGER :: nlines
      
      INTEGER :: nproc

      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ect
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: vct      
      INTEGER, ALLOCATABLE, DIMENSION(:) :: nelnds
      INTEGER :: mnelnds
      INTEGER, ALLOCATABLE, DIMENSION(:) :: el_type    
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: xy
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: depth      
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: elxy
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: elhb   
      INTEGER :: curved_grid

      INTEGER :: nope
      INTEGER, ALLOCATABLE, DIMENSION(:) :: obseg
      INTEGER :: neta
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: obnds

      INTEGER :: nbou
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: fbseg
      INTEGER :: nvel
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: fbnds
      
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ged2nn
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ged2el
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ged2led
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: el2ged
      
      INTEGER :: mnepn ! maximum number of elements per node      
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
      INTEGER :: nred
      INTEGER, ALLOCATABLE, DIMENSION(:) :: redn
      INTEGER :: nbed
      INTEGER, ALLOCATABLE, DIMENSION(:) :: bedn
      
      INTEGER, ALLOCATABLE, DIMENSION(:) :: part

      INTEGER, ALLOCATABLE, DIMENSION(:) :: nresel
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: el_g2l,el_l2g

      INTEGER :: nsred
      INTEGER, ALLOCATABLE, DIMENSION(:) :: sredn
      
      INTEGER :: mned_sr
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ned_sr
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: pe_sr
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: el_sr
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: led_sr
      
      INTEGER, ALLOCATABLE, DIMENSION(:) :: nresnd
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: nd_l2g      
      INTEGER, ALLOCATABLE, DIMENSION(:) :: nd_g2l
       
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: lect
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: lnelnds
      
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: lobseg
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: lobnds
      INTEGER, ALLOCATABLE, DIMENSION(:) :: lneta
      INTEGER, ALLOCATABLE, DIMENSION(:) :: lnope      
      
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: lfbseg
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: lfbnds
      INTEGER, ALLOCATABLE, DIMENSION(:) :: lnvel
      INTEGER, ALLOCATABLE, DIMENSION(:) :: lnbou     
      
      INTEGER :: nobfr ! number of periodic forcings on elevation boundaries
      CHARACTER(10), ALLOCATABLE, DIMENSION(:) :: obtag,obtag2 ! constituent name
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: obfreq ! open bounday constituent frequency
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: obnfact ! open boundary constituent nodal factor
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: obeq ! open bounday constituent equilibrium argument
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: obamp ! open boundary node amplitute
      REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: obph ! open boundary node phase
      
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: lobamp      
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: lobph

      INTEGER :: nfbfr ! number of periodic forcings on flow boundaries
      CHARACTER(10), ALLOCATABLE, DIMENSION(:) :: fbtag,fbtag2 ! constituent name
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: fbfreq ! flow bounday constituent frequency
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: fbnfact ! flow boundary constituent nodal factor
      REAL(rp), ALLOCATABLE, DIMENSION(:) :: fbeq ! flow bounday constituent equilibrium argument
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: fbamp ! flow boundary node amplitute  
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:) :: fbph ! flow boundary node phase

      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:,:) :: lfbamp
      REAL(rp), ALLOCATABLE, DIMENSION(:,:,:,:) :: lfbph 
      INTEGER, ALLOCATABLE, DIMENSION(:) :: lnbouf
      
      END MODULE globals