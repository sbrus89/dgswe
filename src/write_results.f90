      MODULE write_results

      USE globals, ONLY: pres,grid,np,mninds,mnnds

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
      
      OPEN(unit=9,file='../output/interior_nodes.d')
      WRITE(9,*) mesh%ne,mninds
      DO el = 1,mesh%ne
        DO i = 1,mninds
!           WRITE(9,*) (xyhw(i,el,j), j = 1,3)
          WRITE(9,*) mesh%xyhi(i,el,1)/(Erad*cos(phi0))+lambda0, mesh%xyhi(i,el,2)/Erad, mesh%xyhi(i,el,3)
        ENDDO
      ENDDO
      
      CLOSE(9)
      
!       CALL invcpp(ned,np(3)-1,mnnds,xyhe,xyhw)
      
      OPEN(unit=10,file='../output/edge_nodes.d')
      WRITE(10,*) mesh%ned,np(3)-1
      DO ed = 1,mesh%ned      
        DO i = 1,np(3)-1
!           WRITE(10,*) (xyhw(i,ed,j), j = 1,3)
          WRITE(10,*) mesh%xyhe(i,ed,1)/(Erad*cos(phi0))+lambda0, mesh%xyhe(i,ed,2)/Erad, mesh%xyhe(i,ed,3)          
        ENDDO
      ENDDO
      
      CLOSE(10)
      
      OPEN(unit=13,file='../output/boundary_nodes.d')
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
      
      OPEN(unit=11,file='../output/centers.d')
      WRITE(11,*) base%ne    
      DO el = 1,base%ne      
        WRITE(11,*) (base%xyhc(i,el), i = 1,3)
      ENDDO      
      
      CLOSE(11)
      
      OPEN(unit=12,file='../output/normals.d')
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
      
      OPEN(unit=9,file='../output/rimls_interior_nodes.d')
      WRITE(9,*) mesh%ne,mninds
      DO el = 1,mesh%ne
        DO i = 1,mninds
          WRITE(9,*) (mesh%xyhi(i,el,j), j = 1,3)
        ENDDO
      ENDDO
      
      CLOSE(9)      
      
      
      OPEN(unit=10,file='../output/rimls_edge_nodes.d')
      WRITE(10,*) mesh%ned,np(3)-1
      DO ed = 1,mesh%ned      
        DO i = 1,np(3)-1
          WRITE(10,*) (mesh%xyhe(i,ed,j), j = 1,3)
        ENDDO
      ENDDO
      
      CLOSE(10)
      
      OPEN(unit=11,file='../output/rimls_vertex_nodes.d')
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
      
      PRINT "(A)", "Writing fort.14 with rimls nodes..."       
      
      OPEN(UNIT = 14, FILE = "../output/fort.14_rimls")
      
      
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

      SUBROUTINE write_elem_nodes(mesh)
      
      USE globals, ONLY: nverts,ctp,nnds

      IMPLICIT NONE
      
      INTEGER :: el,ed,et,nv,pt,i,pn,led,n,nd,v,nnd
      TYPE(grid) :: mesh
      REAL(pres) :: nodes(mnnds,mesh%ne)
      
      nodes = 0d0
      
      PRINT "(A)", "Writing rimls element nodes..."  
      PRINT*, ""     
      
      ! Verticies
      DO el = 1,mesh%ne
        et = mesh%el_type(el)
        nv = nverts(et)
        
        DO v = 1,nv
          i = (v-1)*ctp + 1
          nd = mesh%ect(v,el)
          nodes(i,el) = mesh%xyhv(1,nd,3)
        ENDDO
      ENDDO
      
      ! Edge points
      DO ed = 1,mesh%ned
        
        el = mesh%ged2el(1,ed)
        led = mesh%ged2led(1,ed)
        et = mesh%el_type(el)
        nv = nverts(et)       
          
        DO pt = 1,ctp-1          
          i = mod(led,nv)*ctp + pt + 1
          nodes(i,el) = mesh%xyhe(pt,ed,3)
        ENDDO          
          
        IF (mesh%bed_flag(ed) == 0) THEN
        
          el = mesh%ged2el(2,ed)
          led = mesh%ged2led(2,ed)
          et = mesh%el_type(el)
          nv = nverts(et)                 
          
          DO pt = 1,ctp-1          
            i = mod(led,nv)*ctp + ctp - pt + 1            
            nodes(i,el) = mesh%xyhe(pt,ed,3)
          ENDDO            

        ENDIF          
          
      ENDDO
      
      ! Interior points
      DO el = 1,mesh%ne
        et = mesh%el_type(el)
        n = nnds(et)
        nv = nverts(et)
        
        IF (mod(et,2) == 1) THEN   
          nnd = nnds(3)
        ELSE IF (mod(et,2) == 0) THEN
          nnd = nnds(4)          
        ENDIF    

        pt = 0
        DO i = nv*(ctp-1)+nv+1,nnd 
          pt = pt + 1
          nodes(i,el) = mesh%xyhi(pt,el,3) 
        ENDDO
      ENDDO
      
      
      OPEN(UNIT = 14, FILE = "../output/elem_nodes.d")      
      
      DO el = 1,mesh%ne
        et = mesh%el_type(el)
        
        IF (mod(et,2) == 1) THEN   
          nnd = nnds(3)
        ELSE IF (mod(et,2) == 0) THEN
          nnd = nnds(4)          
        ENDIF    

        WRITE(14,"(2(I7),6(e24.17,1x))") el,nnd,(nodes(i,el), i = 1,nnd)
      ENDDO
      
      CLOSE(14)


      RETURN
      END SUBROUTINE write_elem_nodes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      END MODULE write_results