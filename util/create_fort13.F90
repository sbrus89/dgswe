      PROGRAM create_fort13

      USE globals, ONLY: rp      
      USE nodal_attributes_mod, ONLY: type13,read_fort13,write_fort13,apply_default_attributes, &
                                      create_neighbor_table,nearest_neighbor_interpolation, &
                                      determine_non_default_nodes,print_fort13_info
      USE grid_file_mod, ONLY: grid_type,read_grid
      USE kdtree2_module

      IMPLICIT NONE

      TYPE(grid_type) :: old14     
      TYPE(type13) :: old13
      
      TYPE(grid_type) :: prvi14
      TYPE(type13) :: prvi13

      TYPE(grid_type) :: new14
      TYPE(grid_type) :: manning14
      TYPE(grid_type) :: submerge14       
      TYPE(type13) :: new13
      
      TYPE(type13) :: it13

      TYPE(kdtree2), POINTER :: tree_xy      
      TYPE(kdtree2_result), ALLOCATABLE, DIMENSION(:) :: closest   
      REAL(rp) :: xy(2)       
      
      INTEGER :: i,j,k,l,attr
      INTEGER :: n,nattr
      INTEGER :: nd
      INTEGER :: clnd      
      INTEGER :: nattr_list,attr_list(7)
      INTEGER, DIMENSION(:), ALLOCATABLE :: neighbor_table
      
      old13%file_name = "/home/sbrus/data-drive/DWH_Oysters/inputs/fort.13/SL18TX33/sl18_tx_v2b.13.advectionstate.GULF"
      old14%grid_file = "/home/sbrus/data-drive/DWH_Oysters/inputs/fort.14/SL18TX33/sl18_tx_v2c_071311.grd"
      
      prvi13%file_name = "/home/sbrus/data-drive/DWH_Oysters/inputs/fort.13/PRVI/PRVI-2017.13.V2"
      prvi14%grid_file = "/home/sbrus/data-drive/DWH_Oysters/inputs/fort.14/PRVI/PRVI-2017.14"
      
      manning14%grid_file = "/home/sbrus/data-drive/DWH_Oysters/inputs/fort.13/SL18TX33+PRVI/SL18_GoM_PRVI_manning.grd"
      submerge14%grid_file = "/home/sbrus/data-drive/DWH_Oysters/inputs/fort.13/SL18TX33+PRVI/SL18_GoM_PRVI_submergence.grd"
      
      it13%file_name = "/home/sbrus/data-drive/DWH_Oysters/inputs/fort.13/SL18TX33+PRVI/fort.13_internal_tides"
      
      new14%grid_file = "/home/sbrus/data-drive/DWH_Oysters/inputs/fort.14/SL18TX33+PRVI/SL18_GoM_PRVI_natural.grd"
      new13%file_name = "/home/sbrus/data-drive/DWH_Oysters/inputs/fort.13/SL18TX33+PRVI/fort.13_fixed_sub"
      
      ! average_horizontal_eddy_viscosity_in_sea_water_wrt_depth
      ! sea_surface_height_above_geoid
      ! surface_submergence_state
      ! surface_canopy_coefficient
      ! primitive_weighting_in_continuity_equation
      ! surface_directional_effective_roughness_length
      ! mannings_n_at_sea_floor
      ! advection_state
      
!       CALL read_fort13(new13)
!       CALL print_fort13_info(new13)
      
!       STOP
      
      CALL read_fort13(old13)
      CALL read_grid(old14,1)
      CALL print_fort13_info(old13)      
      
      CALL read_fort13(prvi13)
      CALL read_grid(prvi14,3)
      CALL print_fort13_info(prvi13)      
      
      CALL read_grid(new14,3)      
      CALL read_grid(manning14,1)
      CALL read_grid(submerge14,1)
   
      CALL read_fort13(it13)
      CALL print_fort13_info(it13)       
      
!       PRINT*, MINVAL(new14%xy(1,:)),MAXVAL(new14%xy(1,:)),MINVAL(new14%xy(2,:)),MAXVAL(new14%xy(2,:))            
!        PRINT*, MINVAL(manning14%depth(:)),MAXVAL(manning14%depth(:))   
!        PRINT*, MINVAL(submerge14%depth(:)),MAXVAL(submerge14%depth(:))  
      
      attr_list = (/ 1, 2, 3, 5, 7, 3, 1 /)
      
!       CALL create_neighbor_table(old14%nn,old14%xy,new14%nn,new14%xy,neighbor_table)      
!       
!       DO k = 1,3
!         IF (k == 2) THEN
!           CYCLE
!         ENDIF
!       
!         j = attr_list(k)
!         DO i = 1,new13%nn
!           nd = neighbor_table(i)
!           
!           IF (abs(new13%defined_val(i,k)-old13%defined_val(nd,j)) > 1d-12) THEN
!             PRINT*, "matching error old13", k, j
!             STOP
!           ENDIF
!         ENDDO        
!       ENDDO
!       
!       k = 4
!       n = 0
!       DO i = 1,new13%nn
!           IF (abs(new13%defined_val(i,k)-manning14%depth(i)) > 1d-12) THEN
! !             PRINT*, new13%defined_val(i,k)
! !             PRINT*, manning14%depth(i)
! !             PRINT*, "matching error manning" 
! !             STOP
!               n = n + 1
!           ENDIF        
!       ENDDO
!       
!       IF (n > 0) THEN
!             PRINT*, "matching error manning", n      
!       ENDIF
!       
!       k = 5
!       n = 0
!       DO i = 1,new13%nn
!           IF (abs(new13%defined_val(i,k)-submerge14%depth(i)) > 1d-12) THEN
! !             PRINT*, "matching error submerge" 
! !             STOP
!               n = n + 1
!           ENDIF        
!       ENDDO      
!       
!       IF (n > 0) THEN
!             PRINT*, "matching error submerge", n      
!       ENDIF      
!       
!       k = 6
!       CALL create_neighbor_table(prvi14%nn,prvi14%xy,new14%nn,new14%xy,neighbor_table)
!       j = attr_list(k)
!       DO i = 1,new13%nn
!         nd = neighbor_table(i)
!           
!         IF (abs(new13%defined_val(i,k)-prvi13%defined_val(nd,j)) > 1d-12) THEN
!           PRINT*, "matching error advection"
!           STOP
!         ENDIF
!       ENDDO    
!       
!       STOP
!       
!       k = 7
!       DO i = 1,new13%nn
!           IF (abs(new13%defined_val(i,k) -it13%defined_val(i,1)) > 1d-12) THEN
!             PRINT*, "matching error internal tides"
!             STOP
!           ENDIF        
!       ENDDO        
!       
!       STOP
      
      ALLOCATE(new13%defined_val(new14%nn,7))
      ALLOCATE(new13%default_val(7))
      ALLOCATE(new13%ndns(new14%nn,7))
      ALLOCATE(new13%nndns(old13%nattr))
      ALLOCATE(new13%name(old13%nattr))
      ALLOCATE(new13%units(old13%nattr))
      ALLOCATE(new13%nvals(old13%nattr))
      ALLOCATE(new13%attr_index(7))
      
      nattr = 5
      new13%nn = new14%nn
      new13%nattr = nattr
      new13%grid_name = new14%grid_name    
      
      new13%n = 0
      DO k = 1,nattr
        j = attr_list(k)
        i = old13%attr_index(j)-1 
!         PRINT*, k,j,i
        new13%name(k)%line = old13%name(j)%line
        new13%units(k)%line = old13%units(j)%line
        new13%nvals(k) = old13%nvals(j)
        DO l = 1,new13%nvals(k)
          new13%default_val(new13%n+l) = old13%default_val(i+l) 
!           PRINT*, new13%default_val(new13%n+l)
        ENDDO
        new13%attr_index(k) = new13%n+1
        new13%n = new13%n + new13%nvals(k)
      ENDDO
      
      
      CALL apply_default_attributes(new13%nn,new13%nattr,new13%default_val,new13%defined_val)         
      
      CALL create_neighbor_table(old14%nn,old14%xy,new14%nn,new14%xy,neighbor_table)
      
      ! Eddy viscosity, Sea surface height above geoid, primative weighting from SL18TX33
      DO k = 1,4
        j = attr_list(k)
        i = old13%attr_index(j)
        CALL nearest_neighbor_interpolation(new13%nn,neighbor_table,old13%defined_val(:,i),old13%default_val(i), &
                                            new13%nndns(k),new13%ndns(:,k),new13%defined_val(:,k))
      ENDDO
      
      new13%nndns(2) = 0 ! make sea_surface_height_above_geoid constant
      
      ! Manning's n from Andika fort.14
      PRINT*, "Finding non-default manning n values"
      k = 5
      j = attr_list(k)
      i = old13%attr_index(j)      
      CALL determine_non_default_nodes(new14%nn,manning14%depth,old13%default_val(i), &
                                       new13%nndns(k),new13%ndns(:,k),new13%defined_val(:,k))                                                      
      PRINT("(2(A,I9,2x))"), "non-default = ", new13%nndns(k), "default = ",new14%nn - new13%nndns(k)                                       
      PRINT*, "" 
      
!       ! Surface submergence state from Andika fort.14
!       PRINT*, "Finding non-default surface submergence state values"
!       k = 5
!       j = attr_list(k)
!       i = old13%attr_index(j) 
!       CALL determine_non_default_nodes(new14%nn,submerge14%depth,old13%default_val(i), &
!                                        new13%nndns(k),new13%ndns(:,k),new13%defined_val(:,k))                                                      
!       PRINT("(2(A,I9,2x))"), "non-default = ", new13%nndns(k), "default = ",new14%nn - new13%nndns(k)                                       
!       PRINT*, ""       
      
      CALL create_neighbor_table(prvi14%nn,prvi14%xy,new14%nn,new14%xy,neighbor_table)
      
      ! Advection state from PRVI grid
      k = 6
      j = attr_list(k)
      i = prvi13%attr_index(j)
      
!       DO nd = 1,new13%nn
!         new13%defined_val(nd,k) = prvi13%default_val(i)
!       ENDDO
      
      CALL nearest_neighbor_interpolation(new13%nn,neighbor_table,prvi13%defined_val(:,i),prvi13%default_val(i), &
                                          new13%nndns(k),new13%ndns(:,k),new13%defined_val(:,k))
                   
      new13%n = new13%n + 1
      new13%nattr = new13%nattr + 1                   
      new13%name(k)%line = prvi13%name(j)%line
      new13%units(k)%line = prvi13%units(j)%line
      new13%nvals(k) = prvi13%nvals(j)
      new13%default_val(new13%n) = prvi13%default_val(i) 
      
      ! Internal tides
      k = 7
      
      new13%n = new13%n + 1
      new13%nattr = new13%nattr + 1                   
      new13%name(k)%line = it13%name(1)%line
      new13%units(k)%line = it13%units(1)%line
      new13%nvals(k) = it13%nvals(1)
      new13%default_val(new13%n) = it13%default_val(1)   
      new13%nndns(k) = it13%nndns(1)
      
      DO i = 1,it13%nndns(1)
        nd = it13%ndns(i,1)
        new13%ndns(i,k) = nd
        new13%defined_val(nd,k) = it13%defined_val(nd,1)
      ENDDO
      

      CALL write_fort13(new13)
      CALL print_fort13_info(new13)
      
!       PRINT*, new14%xy(1,new14%obnds(1,1)),new14%xy(2,new14%obnds(1,1))
!       PRINT*, prvi14%xy(1,prvi14%obnds(1,1)),prvi14%xy(2,prvi14%obnds(1,1))      

      END PROGRAM create_fort13