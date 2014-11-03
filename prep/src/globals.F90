      MODULE globals

      IMPLICIT NONE

      SAVE
      
      INTEGER, PARAMETER :: pres = kind(1d0)

      CHARACTER(50) :: grid_file
      CHARACTER(50) :: grid_name
      CHARACTER(50) :: forcing_file
      CHARACTER(50) :: out_direc

      INTEGER :: ne
      INTEGER :: nn
      
      INTEGER :: p
      INTEGER :: ctp
      REAL(pres) :: dt
      REAL(pres) :: tf
      REAL(pres) :: dramp
      REAL(pres) :: cf
      REAL(pres) :: lines
      INTEGER :: nblk
      INTEGER :: npart
      INTEGER :: mndof
      
      INTEGER :: nproc

      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ect
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: xy
      INTEGER, ALLOCATABLE, DIMENSION(:) :: depth

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
      
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: lobseg
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: lobnds
      INTEGER, ALLOCATABLE, DIMENSION(:) :: lneta
      INTEGER, ALLOCATABLE, DIMENSION(:) :: lnope      
      
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: lfbseg
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: lfbnds
      INTEGER, ALLOCATABLE, DIMENSION(:) :: lnvel
      INTEGER, ALLOCATABLE, DIMENSION(:) :: lnbou     
      
      INTEGER :: nobfr ! number of periodic forcings on elevation boundaries
      CHARACTER(10), ALLOCATABLE, DIMENSION(:) :: obtag,obtag2 ! constituent name
      REAL(pres), ALLOCATABLE, DIMENSION(:) :: obfreq ! open bounday constituent frequency
      REAL(pres), ALLOCATABLE, DIMENSION(:) :: obnfact ! open boundary constituent nodal factor
      REAL(pres), ALLOCATABLE, DIMENSION(:) :: obeq ! open bounday constituent equilibrium argument
      REAL(pres), ALLOCATABLE, DIMENSION(:,:) :: obamp ! open boundary node amplitute
      REAL(pres), ALLOCATABLE, DIMENSION(:,:) :: obph ! open boundary node phase
      
      REAL(pres), ALLOCATABLE, DIMENSION(:,:,:) :: lobamp      
      REAL(pres), ALLOCATABLE, DIMENSION(:,:,:) :: lobph

      INTEGER :: nfbfr ! number of periodic forcings on flow boundaries
      CHARACTER(10), ALLOCATABLE, DIMENSION(:) :: fbtag,fbtag2 ! constituent name
      REAL(pres), ALLOCATABLE, DIMENSION(:) :: fbfreq ! flow bounday constituent frequency
      REAL(pres), ALLOCATABLE, DIMENSION(:) :: fbnfact ! flow boundary constituent nodal factor
      REAL(pres), ALLOCATABLE, DIMENSION(:) :: fbeq ! flow bounday constituent equilibrium argument
      REAL(pres), ALLOCATABLE, DIMENSION(:,:,:) :: fbamp ! flow boundary node amplitute  
      REAL(pres), ALLOCATABLE, DIMENSION(:,:,:) :: fbph ! flow boundary node phase

      REAL(pres), ALLOCATABLE, DIMENSION(:,:,:,:) :: lfbamp
      REAL(pres), ALLOCATABLE, DIMENSION(:,:,:,:) :: lfbph 
      INTEGER, ALLOCATABLE, DIMENSION(:) :: lnbouf
      
      END MODULE globals