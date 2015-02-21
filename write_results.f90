      MODULE write_results

      USE globals, ONLY: grid,np,mninds

      CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE write_linear_nodes(mesh)

      USE globals, ONLY: Erad,phi0,lambda0
      
      IMPLICIT NONE
      
      INTEGER :: el,ed,i,j
      TYPE(grid) :: mesh
      
!       ALLOCATE(xyhw(mnnds,mesh%ned,3))      

!       CALL invcpp(ne,mninds,mnnds,xyhi,xyhw)

      PRINT "(A)", "Writing linearly interpolated nodes..."   
      PRINT*, ""
      
      OPEN(unit=9,file='interior_nodes.d')
      WRITE(9,*) mesh%ne,mninds
      DO el = 1,mesh%ne
        DO i = 1,mninds
!           WRITE(9,*) (xyhw(i,el,j), j = 1,3)
          WRITE(9,*) mesh%xyhi(i,el,1)/(Erad*cos(phi0))+lambda0, mesh%xyhi(i,el,2)/Erad, mesh%xyhi(i,el,3)
        ENDDO
      ENDDO
      
      CLOSE(9)
      
!       CALL invcpp(ned,np(3)-1,mnnds,xyhe,xyhw)
      
      OPEN(unit=10,file='edge_nodes.d')
      WRITE(10,*) mesh%ned,np(3)-1
      DO ed = 1,mesh%ned      
        DO i = 1,np(3)-1
!           WRITE(10,*) (xyhw(i,ed,j), j = 1,3)
          WRITE(10,*) mesh%xyhe(i,ed,1)/(Erad*cos(phi0))+lambda0, mesh%xyhe(i,ed,2)/Erad, mesh%xyhe(i,ed,3)          
        ENDDO
      ENDDO
      
      CLOSE(10)
      
      OPEN(unit=13,file='boundary_nodes.d')
      WRITE(13,*) mesh%ned,np(3)-1
      DO ed = 1,mesh%ned      
        DO i = 1,np(3)-1
          WRITE(13,*) mesh%bnd_flag(i,ed)        
        ENDDO
      ENDDO
      
      CLOSE(13)         

      RETURN
      END SUBROUTINE write_linear_nodes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE write_normals()

      USE globals, ONLY: base
      
      IMPLICIT NONE
      
      INTEGER :: el,i
      
      PRINT*, "" 
      PRINT "(A)", "Writing element normals..."             
      
      OPEN(unit=11,file='centers.d')
      WRITE(11,*) base%ne    
      DO el = 1,base%ne      
        WRITE(11,*) (base%xyhc(i,el), i = 1,3)
      ENDDO      
      
      CLOSE(11)
      
      OPEN(unit=12,file='normals.d')
      WRITE(12,*) base%ne
      DO el = 1,base%ne
        WRITE(12,*) (base%nhb(i,el),i=1,3)
      ENDDO
      
      CLOSE(12)
      
      RETURN
      END SUBROUTINE write_normals

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE write_rimls_nodes(mesh)

      IMPLICIT NONE
      
      INTEGER :: el,nd,ed,i,j
      TYPE(grid) :: mesh

      PRINT*, ""
      PRINT "(A)", "Writing rimls surface nodes..."  
      PRINT*, ""      
      
      OPEN(unit=9,file='rimls_interior_nodes.d')
      WRITE(9,*) mesh%ne,mninds
      DO el = 1,mesh%ne
        DO i = 1,mninds
          WRITE(9,*) (mesh%xyhi(i,el,j), j = 1,3)
        ENDDO
      ENDDO
      
      CLOSE(9)      
      
      
      OPEN(unit=10,file='rimls_edge_nodes.d')
      WRITE(10,*) mesh%ned,np(3)-1
      DO ed = 1,mesh%ned      
        DO i = 1,np(3)-1
          WRITE(10,*) (mesh%xyhe(i,ed,j), j = 1,3)
        ENDDO
      ENDDO
      
      CLOSE(10)
      
      OPEN(unit=11,file='rimls_vertex_nodes.d')
      WRITE(11,*) mesh%nn,1
      DO nd = 1,mesh%nn
          WRITE(11,*) (mesh%xyhv(1,nd,j), j = 1,3)
      ENDDO
      
      CLOSE(11)        

      RETURN
      END SUBROUTINE write_rimls_nodes
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE rewrite_fort14(mesh)

      IMPLICIT NONE
      
      INTEGER :: i,j
      TYPE(grid) :: mesh
      
      PRINT*, ""
      PRINT "(A)", "Writing fort.14 with rimls nodes..."  
      PRINT*, ""            
      
      OPEN(UNIT = 14, FILE = "fort.14_rimls")
      
      
      WRITE(14,*) mesh%grid_name      
      WRITE(14,*) mesh%ne,mesh%nn         
      
      DO i = 1,mesh%nn
        WRITE(14,*) i, (mesh%xyhv(1,i,j), j = 1,3)  ! write new rimls x,y, depth coordinates
      ENDDO
      
      DO i = 1,mesh%ne
        WRITE(14,*) i, mesh%nelnds(i), (mesh%ect(j,i), j = 1,mesh%nelnds(i))               
      ENDDO
      
      WRITE(14,*) mesh%nope
      WRITE(14,*) mesh%neta
      
      DO i = 1,mesh%nope
        WRITE(14,*) mesh%obseg(i)
        DO j = 1,mesh%obseg(i)
          WRITE(14,*) mesh%obnds(j,i)
        ENDDO
      ENDDO
      
      WRITE(14,*) mesh%nbou      
      WRITE(14,*) mesh%nvel
      
      DO i = 1,mesh%nbou
        WRITE(14,*) mesh%fbseg(1,i), mesh%fbseg(2,i)
        DO j = 1,mesh%fbseg(1,i)
          WRITE(14,*) mesh%fbnds(j,i)
        ENDDO
      ENDDO
      
      CLOSE(14)
      
      



      RETURN
      END SUBROUTINE rewrite_fort14

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      END MODULE write_results