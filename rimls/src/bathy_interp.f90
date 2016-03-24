      SUBROUTINE bathymetry_eval_nodes(mesh)

      USE globals, ONLY: rp,grid,nel_type,np,nnds,nverts
      USE bathymetry_interp_mod               

      IMPLICIT NONE
      
      TYPE(grid) :: mesh
      REAL(rp), DIMENSION(:,:,:), ALLOCATABLE :: psic

      
      CALL shape_functions_eltypes_at_hbp(nel_type,np,psic) 
      
      CALL bathy_coordinates(mesh%ne,nnds,nverts,mesh%el_type,mesh%elxy,psic,mesh%elxyh, &
                             mesh%depth,mesh%ect,mesh%elhb)


      RETURN
      END SUBROUTINE bathymetry_eval_nodes
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!          
