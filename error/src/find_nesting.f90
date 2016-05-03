      SUBROUTINE find_nesting(coarse,fine,elf2elc)

      USE globals, ONLY: rp,solution,nverts,nel_type
      USE find_element, ONLY: find_element_init,in_element,find_element_final
      USE transformation, ONLY: element_transformation

      IMPLICIT NONE
      
      TYPE(solution) :: coarse,fine
      
      INTEGER :: elf,nd
      INTEGER :: etf,elin
      INTEGER :: n1,n2
      INTEGER :: el_found,found
      INTEGER :: elf2elc(fine%ne)
      
      REAL(rp) :: xf(2)   
      REAL(rp) :: rs(2)
      
      CALL find_element_init(nel_type,nverts,coarse%np,coarse%nnds,coarse%nn,coarse%xy,coarse%nepn,coarse%epn)  
      
      elf2elc(:) = 0
      el_found = 0 
      
      PRINT("(A)"), TRIM(coarse%sol_name) // " / " // TRIM(fine%sol_name) // " grid"
      
      ! Find coarse element each fine element is located in
elemf:DO elf = 1,fine%ne
           
        etf = fine%el_type(elf)
        
        ! Compute the (x,y) coordinates of the first quadrature point (from fine element)         
        CALL element_transformation(fine%nnds(etf),fine%elxy(:,elf,1),fine%elxy(:,elf,2),fine%l(:,1,etf),xf(1),xf(2))
                        
        CALL in_element(xf,coarse%el_type,coarse%elxy,elin,rs)          
       
        elf2elc(elf) = elin            
                                                              
      ENDDO elemf
      

      CALL find_element_final()
      
      
      RETURN
      END SUBROUTINE find_nesting