      SUBROUTINE curvilinear(mesh)

      USE globals, ONLY: rp,grid,nel_type,ctp, &
                         np,nnds,nverts
      USE curvilinear_nodes_mod                         

      IMPLICIT NONE
      
      TYPE(grid) :: mesh    
      REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: psiv
      REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: psic
      
      
      CALL shape_functions_linear_at_ctp(nel_type,np,psiv)             
      
      CALL eval_coordinates_curved(ctp,nnds,nverts,mesh%el_type,mesh%xy,mesh%ect,mesh%fbseg,mesh%fbnds, &
                                   mesh%nnfbed,mesh%nfbedn,mesh%nfbednn,mesh%ged2el,mesh%ged2led, &
                                   psiv,mesh%bndxy,mesh%elxy)     
                                   
                                   
      CALL shape_functions_eltypes_at_hbp(nel_type,np,psic)                                      
      
      CALL bathy_coordinates(mesh%ne,nnds,mesh%el_type,mesh%elxy,psic,mesh%elxyh)      
      
   
      
      
      RETURN
      END SUBROUTINE curvilinear
