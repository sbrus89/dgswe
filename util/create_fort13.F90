      PROGRAM create_fort13

      USE nodal_attributes_mod, ONLY: type13,read_fort13

      IMPLICIT NONE

      TYPE(type13) :: old_grid

      
      old_grid%file_name = "/home/sbrus/data-drive/DWH_Oysters/grid/sl18_tx_v2b.13.advectionstate.GULF"
      
      CALL read_fort13(old_grid)



      END PROGRAM create_fort13