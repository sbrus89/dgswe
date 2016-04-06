      MODULE write_results

      USE globals, ONLY: rp,grid,out_direc,nverts

      CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE write_linear_nodes(mesh)

      USE globals, ONLY: Erad,phi0,lambda0
      
      IMPLICIT NONE
      
      INTEGER :: el,ed,i,j
      TYPE(grid) :: mesh
    

      PRINT "(A)", "Writing linearly interpolated nodes..."   
      PRINT*, ""
      
      OPEN(unit=9,file=TRIM(out_direc) // 'interior_nodes.d')
      WRITE(9,*) mesh%tpts_interior
      DO el = 1,mesh%ne
        DO i = 1,mesh%npts_interior(el)
          WRITE(9,*) mesh%xyh_interior(i,el,1), mesh%xyh_interior(i,el,2), mesh%xyh_interior(i,el,3)
        ENDDO
      ENDDO
      
      CLOSE(9)
      
      OPEN(unit=10,file=TRIM(out_direc) // 'edge_nodes.d')
      WRITE(10,*) mesh%tpts_edge
      DO ed = 1,mesh%ned      
        DO i = 1,mesh%npts_edge(ed)
          WRITE(10,*) mesh%xyh_edge(i,ed,1), mesh%xyh_edge(i,ed,2), mesh%xyh_edge(i,ed,3)          
        ENDDO
      ENDDO
      
      CLOSE(10)
      
      OPEN(unit=13,file=TRIM(out_direc) // 'boundary_nodes.d')
      WRITE(13,*) mesh%tpts_edge
      DO ed = 1,mesh%ned      
        DO i = 1,mesh%npts_edge(ed)
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
      
      INTEGER :: el,i,j
      
      PRINT*, "" 
      PRINT "(A)", "Writing element normals..."       
      

      
      OPEN(unit=11,file=TRIM(out_direc) // 'centers.d')
      WRITE(11,*) base%tpts_interior    
      DO el = 1,base%ne      
        DO i = 1,base%npts_interior(el)
          WRITE(11,*) base%elxyh(i,el,1),base%elxyh(i,el,2),base%elhb(i,el)
        ENDDO
      ENDDO      
      
      CLOSE(11)
      
      OPEN(unit=12,file=TRIM(out_direc) // 'normals.d')
      WRITE(12,*) base%ne
      DO el = 1,base%ne
        DO i = 1,base%mninds
          WRITE(12,*) (base%nhb(i,el,j),j=1,3)
        ENDDO
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
      
      OPEN(unit=9,file=TRIM(out_direc) // 'rimls_interior_nodes.d')
      WRITE(9,*) mesh%tpts_interior
      DO el = 1,mesh%ne
        DO i = 1,mesh%npts_interior(el)
          WRITE(9,*) (mesh%xyh_interior(i,el,j), j = 1,3)
        ENDDO
      ENDDO
      
      CLOSE(9)      
      
      
      OPEN(unit=10,file=TRIM(out_direc) // 'rimls_edge_nodes.d')
      WRITE(10,*) mesh%tpts_edge
      DO ed = 1,mesh%ned      
        DO i = 1,mesh%npts_edge(ed)
          WRITE(10,*) (mesh%xyh_edge(i,ed,j), j = 1,3)
        ENDDO
      ENDDO
      
      CLOSE(10)
      
      OPEN(unit=11,file=TRIM(out_direc) // 'rimls_vertex_nodes.d')
      WRITE(11,*) mesh%nn
      DO nd = 1,mesh%nn
          WRITE(11,*) (mesh%xyh_vertex(1,nd,j), j = 1,3)
      ENDDO
      
      CLOSE(11)        

      RETURN
      END SUBROUTINE write_rimls_nodes
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE rewrite_fort14(mesh)

      IMPLICIT NONE
      
      INTEGER :: i,j
      INTEGER :: et,nv
      TYPE(grid) :: mesh
      
      PRINT "(A)", "Writing fort.14 with rimls nodes..."       
      
      OPEN(UNIT = 14, FILE = TRIM(out_direc) // "fort.14_rimls")
      
      
      WRITE(14,"(A)") mesh%grid_name      
      WRITE(14,*) mesh%ne,mesh%nn         
      
      DO i = 1,mesh%nn
        WRITE(14,"(I7,1x,3(e24.17,1x))") i, (mesh%xyh_vertex(1,i,j), j = 1,3)  ! write new rimls x,y, depth coordinates
      ENDDO
      
      DO i = 1,mesh%ne
        et = mesh%el_type(i)
        nv = nverts(et)
        WRITE(14,*) i, nv, (mesh%ect(j,i), j = 1,nv)               
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
      
      USE globals, ONLY: nverts

      IMPLICIT NONE
      
      INTEGER :: el,ed,et,nv,pt,i,pn,led,n,nd,v,nnd,hbp
      INTEGER :: sind,eind
      CHARACTER(100) :: name
      CHARACTER(1) :: hbp_char
      TYPE(grid) :: mesh
      REAL(rp) :: nodes(mesh%mnnds,mesh%ne)
      REAL(rp) :: xnodes(mesh%mnnds,mesh%ne)
      REAL(rp) :: ynodes(mesh%mnnds,mesh%ne)      
      
      nodes = 0d0
      
      PRINT "(A)", "Writing rimls element nodes..."  
      PRINT*, ""     
      
      hbp = mesh%hbp
      
      ! Verticies
      DO el = 1,mesh%ne
        et = mesh%el_type(el)
        nv = nverts(et)
        
        DO v = 1,nv
          i = (v-1)*hbp + 1
          nd = mesh%ect(v,el)
          nodes(i,el) = mesh%xyh_vertex(1,nd,3)
          xnodes(i,el) = mesh%xyh_vertex(1,nd,1)
          ynodes(i,el) = mesh%xyh_vertex(1,nd,2)
        ENDDO
      ENDDO
      
      ! Edge points
      DO ed = 1,mesh%ned
        
        el = mesh%ged2el(1,ed)
        led = mesh%ged2led(1,ed)
        et = mesh%el_type(el)
        nv = nverts(et)       
          
        DO pt = 1,hbp-1          
          i = mod(led,nv)*hbp + pt + 1
          nodes(i,el) = mesh%xyh_edge(pt,ed,3)
          xnodes(i,el) = mesh%xyh_edge(pt,ed,1)
          ynodes(i,el) = mesh%xyh_edge(pt,ed,2)          
          
!           PRINT*, el,led,nodes(i,el)
        ENDDO          
          
        IF (mesh%ed_type(ed) == 0) THEN
        
          el = mesh%ged2el(2,ed)
          led = mesh%ged2led(2,ed)
          et = mesh%el_type(el)
          nv = nverts(et)                 
          
          DO pt = 1,hbp-1          
            i = mod(led,nv)*hbp + hbp - pt + 1     
!             i = mod(led,nv)*hbp + pt + 1
            nodes(i,el) = mesh%xyh_edge(pt,ed,3)
            xnodes(i,el) = mesh%xyh_edge(pt,ed,1)
            ynodes(i,el) = mesh%xyh_edge(pt,ed,2)            
!             PRINT*, el,led,nodes(i,el)            
          ENDDO            

        ENDIF   
        
!         PRINT*, " "
          
      ENDDO
      
      ! Interior points
      DO el = 1,mesh%ne
        et = mesh%el_type(el)
        n = mesh%nnds(et)
        nv = nverts(et)
        
        IF (mod(et,2) == 1) THEN   
          nnd = mesh%nnds(3)
        ELSE IF (mod(et,2) == 0) THEN
          nnd = mesh%nnds(4)          
        ENDIF    

        pt = 0
        DO i = nv*(hbp-1)+nv+1,nnd 
          pt = pt + 1
          nodes(i,el) = mesh%xyh_interior(pt,el,3) 
          xnodes(i,el) = mesh%xyh_interior(pt,el,1)
          ynodes(i,el) = mesh%xyh_interior(pt,el,2)
        ENDDO
      ENDDO
      
      sind = INDEX(ADJUSTL(TRIM(mesh%grid_file)),"/",.true.)
      eind = INDEX(ADJUSTL(TRIM(mesh%grid_file)),".",.false.)
      
!       name = ADJUSTL(TRIM(mesh%grid_file(sind+1:eind-1)))     
!       OPEN(UNIT = 13, FILE = TRIM(out_direc) // ADJUSTL(TRIM(name)) // ".hb")    

      WRITE(hbp_char,"(I1)") hbp
      
      name = ADJUSTL(TRIM(mesh%grid_file(1:eind-1)))     
      OPEN(UNIT = 13, FILE = ADJUSTL(TRIM(name)) // "_hbp" // hbp_char // "_rimls.hb")            
      WRITE(13,"(2(I7,1x))") mesh%ne,hbp
      
      OPEN(UNIT = 14, FILE = TRIM(out_direc) // "elem_nodes.d")      
      WRITE(14,"(2(I7,1x))") mesh%ne,hbp
      
      OPEN(UNIT = 15, FILE = TRIM(out_direc) // "xelem_nodes.d")      
      WRITE(15,"(2(I7,1x))") mesh%ne,hbp
      
      OPEN(UNIT = 16, FILE = TRIM(out_direc) // "yelem_nodes.d")      
      WRITE(16,"(2(I7,1x))") mesh%ne,hbp      
      
      DO el = 1,mesh%ne
        et = mesh%el_type(el)
        
        IF (mod(et,2) == 1) THEN   
          nnd = mesh%nnds(3)
        ELSE IF (mod(et,2) == 0) THEN
          nnd = mesh%nnds(4)          
        ENDIF    

        WRITE(13,"(2(I7),1x,60(e24.17,1x))") el,nnd,(nodes(i,el), i = 1,nnd)        
        WRITE(14,"(2(I7),1x,60(e24.17,1x))") el,nnd,(nodes(i,el), i = 1,nnd)
        WRITE(15,"(2(I7),1x,60(e24.17,1x))") el,nnd,(xnodes(i,el), i = 1,nnd)
        WRITE(16,"(2(I7),1x,60(e24.17,1x))") el,nnd,(ynodes(i,el), i = 1,nnd)
      ENDDO
      
      CLOSE(14)
      CLOSE(15)
      CLOSE(16)


      RETURN
      END SUBROUTINE write_elem_nodes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      END MODULE write_results