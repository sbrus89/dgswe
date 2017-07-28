      PROGRAM essential2natural

      USE grid_file_mod

      IMPLICIT NONE
 
      TYPE(grid_type) :: mesh

      INTEGER :: bou
      INTEGER :: btype

      mesh%grid_file = "/home/sbrus/data-drive/DWH_Oysters/inputs/fort.14/SL18TX33+PRVI/SL18_GoM_PRVI_rivers.grd"

      CALL read_grid(mesh,3)

      DO bou = 1,mesh%nbou
        
        btype = mesh%fbseg(2,bou)
        
        IF (btype == 0 .OR. btype == 10) THEN
          btype = 20      
          PRINT*, "Change to ", btype, " on boundary ", bou
        ELSE IF (btype == 1 .OR. btype == 11) THEN     
          btype = 21        
          PRINT*, "Change to ", btype, " on boundary ", bou          
        ELSE IF (btype == 2 .OR. btype == 12) THEN     
          btype = 22        
          PRINT*, "Change to ", btype, " on boundary ", bou          
        ELSE IF (btype == 3 .OR. btype == 13) THEN
          btype = 23
          PRINT*, "Change to ", btype, " on boundary ", bou          
        ELSE IF (btype == 4 .OR. btype == 14) THEN
          btype = 24       
          PRINT*, "Change to ", btype, " on boundary ", bou          
        ELSE IF (btype == 5 .OR. btype == 15) THEN
          PRINT*, "Change to ", btype, " on boundary ", bou        
          btype = 25                
        ENDIF

        mesh%fbseg(2,bou) = btype
        
      ENDDO

      mesh%grid_file = "/home/sbrus/data-drive/DWH_Oysters/inputs/fort.14/SL18TX33+PRVI/SL18_GoM_PRVI_natural.grd"
      
      CALL write_grid(mesh)
      
      END PROGRAM essential2natural