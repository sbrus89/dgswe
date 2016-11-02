      PROGRAM convert_bou2sta

      USE globals, ONLY: rp,ne,nn,nverts,ect,xy,depth,el_type, &
                         nope,neta,obseg,obnds,nvel,nbou,fbseg,fbnds,grid_name
      USE grid_file_mod

      IMPLICIT NONE
      
      INTEGER :: i,j
      INTEGER :: bou,n,nbnds
      INTEGER :: navg,btype
      INTEGER :: nbou_mod,nvel_mod
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: fbseg_mod
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: fbnds_mod
      CHARACTER(100) :: grid_file_in,grid_file_out
      REAL(rp) :: sum_w
      
      REAL(rp), DIMENSION(:), ALLOCATABLE :: w
      REAL(rp), DIMENSION(:,:), ALLOCATABLE :: sta_xy
      
      grid_file_in = "/home/sbrus/data-drive/galveston_SL18/grid_dev/v18_cart/galveston_SL18_cart.grd"
!       grid_file_out = "/home/sbrus/data-drive/galveston_spline_flux/grids/galveston_tri_x16.grd"
      
      navg = 1
      btype = 30
      
      nverts(1) = 3
      nverts(2) = 4
      nverts(3) = 3
      nverts(4) = 4            

      CALL read_header(0,grid_file_in,grid_name,ne,nn)        
      CALL read_coords(nn,xy,depth)
      CALL read_connectivity(ne,ect,el_type)                 
      CALL read_open_boundaries(nope,neta,obseg,obnds)            
      CALL read_flow_boundaries(nbou,nvel,fbseg,fbnds)
      
      
      bou = 0
      DO i = 1,nbou
        IF (fbseg(2,i) == btype) THEN
          bou = i
          EXIT
        ENDIF
      ENDDO
      
      nbnds = fbseg(1,bou)
      
      ALLOCATE(sta_xy(2,nbnds))
      ALLOCATE(w(2*navg+1))
      
      w(1) = 1d0
      w(2) = 1d0
      w(3) = 1d0
      
      sum_w = 0d0
      DO i = 1,2*navg+1
        sum_w = sum_w + w(i)
      ENDDO
      
      DO i = 1,navg
        n = fbnds(i,bou)
        sta_xy(1,i) = xy(1,n)
        sta_xy(2,i) = xy(2,n)
      ENDDO
      DO i = navg+1,nbnds-navg      
        sta_xy(1,i) = 0d0
        sta_xy(2,i) = 0d0
        DO j = 1,2*navg+1
          n = fbnds(i-navg+(j-1),bou)        
          sta_xy(1,i) = sta_xy(1,i) + w(j)*xy(1,n)
          sta_xy(2,i) = sta_xy(2,i) + w(j)*xy(2,n)
        ENDDO
        sta_xy(1,i) = sta_xy(1,i)/sum_w
        sta_xy(2,i) = sta_xy(2,i)/sum_w
      ENDDO
      DO i = 0,navg-1
        n = fbnds(nbnds-i,bou)
        sta_xy(1,nbnds-i) = xy(1,n)
        sta_xy(2,nbnds-i) = xy(2,n)
      ENDDO
      
      OPEN(unit=15,file='stations_avg.d')
      WRITE(15,*) nbnds
      DO i = 1,nbnds
!         n = fbnds(i,bou)
!         WRITE(15,"(2(E25.17,1x))") xy(1,n),xy(2,n)
        WRITE(15,"(2(E25.17,1x))") sta_xy(1,i),sta_xy(2,i)
      ENDDO

!       CALL write_header(grid_file_out,grid_name,ne,nn)
!       CALL write_coords(nn,xy,depth)
!       CALL write_connectivity(ne,ect,el_type,nverts)
!       CALL write_open_boundaries(nope,neta,obseg,obnds)
!       CALL write_flow_boundaries(nbou_mod,nvel_mod,fbseg_mod,fbnds_mod)

      END PROGRAM convert_bou2sta