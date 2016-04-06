      SUBROUTINE curvilinear(mesh)

      USE globals, ONLY: rp,grid,nel_type, &
                         nverts
      USE curvilinear_nodes_mod                         

      IMPLICIT NONE
      
      TYPE(grid) :: mesh    
      REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: psiv
      
      
      CALL shape_functions_linear_at_ctp(nel_type,mesh%np,psiv)             
      
      CALL eval_coordinates_curved(mesh%ctp,mesh%nnds,nverts,mesh%el_type,mesh%xy,mesh%ect,mesh%fbseg,mesh%fbnds, &
                                   mesh%nnfbed,mesh%nfbedn,mesh%nfbednn,mesh%ged2el,mesh%ged2led, &
                                   psiv,mesh%bndxy,mesh%elxy)                                                                               
      

      RETURN
      END SUBROUTINE curvilinear
