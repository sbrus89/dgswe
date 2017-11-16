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

     SUBROUTINE write_curved_coordinates(mesh)
     
     IMPLICIT NONE
     
     TYPE(grid) :: mesh 
     INTEGER :: i,ed,pt
     INTEGER :: ged,el,et
     INTEGER :: nlist,nnodes
     INTEGER :: found
     INTEGER :: el_list(mesh%ne)    
     
     nlist = 0
     nnodes = 0
     DO ed = 1,mesh%nnfbed
       ged = mesh%nfbedn(ed)
       el = mesh%ged2el(1,ged)
       et = mesh%el_type(el)
       
       found = 0
search:DO i = 1,nlist
         IF (el == el_list(i)) THEN
           found = 1           
           EXIT search       
         ENDIF
       ENDDO search
       
       IF (found == 0) THEN
         nlist = nlist + 1     
         el_list(nlist) = el
         
         nnodes = nnodes + mesh%nnds(et)
       ENDIF
     ENDDO
     
     OPEN(unit=10,file=TRIM(out_direc) // 'curved_element_nodes.d')
     WRITE(10,*) nnodes
     DO i = 1,nlist
       el = el_list(i)
       et = mesh%el_type(el)
       DO pt = 1,mesh%nnds(et)
         WRITE(10,*) mesh%elxy(pt,el,1), mesh%elxy(pt,el,2)
       ENDDO
     ENDDO
     CLOSE(10)
     
     RETURN
     END SUBROUTINE write_curved_coordinates


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE write_normals()

      USE globals, ONLY: base
      
      IMPLICIT NONE
      
      INTEGER :: el,i,j
      
      PRINT*, "" 
      PRINT "(A)", "Writing element normals..."       
      

      
      OPEN(unit=11,file=TRIM(out_direc) // 'base_nodes.d')
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
      
      OPEN(unit=13,file=TRIM(out_direc) // 'base_centers.d')
      WRITE(13,*) base%ne    
      DO el = 1,base%ne      
        WRITE(13,*) base%elxy_center(1,el,1),base%elxy_center(1,el,2)
      ENDDO      
      
      CLOSE(13)      
      
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
      
      USE grid_file_mod, ONLY: write_header,write_coords,write_connectivity, &
                               write_open_boundaries,write_flow_boundaries,copy_footer    
      USE version, ONLY: version_information

      IMPLICIT NONE
      
      INTEGER :: i,j,nd,bnd
      INTEGER :: et,nv
      TYPE(grid) :: mesh
      INTEGER :: nbseg,btype
      
      PRINT "(A)", "Writing fort.14 with rimls nodes..."       
      
      DO bnd = 1,mesh%nbou
        btype = mesh%fbseg(2,bnd)
        
        IF (btype == 1 .OR. btype == 11 .OR. btype == 21) THEN
          mesh%fbseg(1,bnd) = mesh%fbseg(1,bnd) - 1
        ENDIF        
        
      ENDDO      
      
      CALL write_header(TRIM(out_direc)//"fort.14_rimls",mesh%grid_name,mesh%ne,mesh%nn)  
      
      CALL write_coords(mesh%nn,mesh%xy,mesh%xyh_vertex(1,:,3))
      
      CALL write_connectivity(mesh%ne,mesh%ect,mesh%el_type,nverts)
      
      CALL write_open_boundaries(mesh%nope,mesh%neta,mesh%obseg,mesh%obnds)
      
      CALL write_flow_boundaries(mesh%nbou,mesh%nvel,mesh%fbseg,mesh%fbnds)
      
      CALL copy_footer(mesh%grid_file,TRIM(out_direc)//"fort.14_rimls")         
      
      OPEN(UNIT=40, FILE=TRIM(out_direc)//"fort.14_rimls", POSITION="APPEND")
      CALL version_information(40)
      
      WRITE(40,"(A)") "-----------------------------------------------------------------------"           
      
      CALL write_input(40)         
      CLOSE(40)      
      
!       OPEN(UNIT = 14, FILE = TRIM(out_direc) // "fort.14_rimls")
!       
!       
!       WRITE(14,"(A)") mesh%grid_name      
!       WRITE(14,*) mesh%ne,mesh%nn         
!       
!       DO i = 1,mesh%nn
!         WRITE(14,"(I7,1x,3(e24.17,1x))") i, mesh%xy(1,i), mesh%xy(2,i), mesh%xyh_vertex(1,i,3)  ! write new rimls x,y, depth coordinates
!       ENDDO
!       
!       DO i = 1,mesh%ne
!         et = mesh%el_type(i)
!         nv = nverts(et)
!         WRITE(14,*) i, nv, (mesh%ect(j,i), j = 1,nv)               
!       ENDDO
!       
!       WRITE(14,"(I9,A)") mesh%nope, "   = Number of open boundaries"
!       WRITE(14,"(I9,A)") mesh%neta, "   = Total number of open boundary nodes"
!       DO bnd = 1,mesh%nope
!         WRITE(14,"(I9)") mesh%obseg(bnd)
!         DO nd = 1,mesh%obseg(bnd)
!           WRITE(14,"(I9)") mesh%obnds(nd,bnd)
!         ENDDO
!       ENDDO
!       
!       WRITE(14,"(I9,A)") mesh%nbou, "   = Number of land boundaries"
!       WRITE(14,"(I9,A)") mesh%nvel, "   = Total number of land boundary nodes"
!       DO bnd = 1,mesh%nbou
!         btype = mesh%fbseg(2,bnd)
!         nbseg = mesh%fbseg(1,bnd)
!         
!         IF (btype == 1 .OR. btype == 11 .OR. btype == 21) THEN
!           nbseg = nbseg - 1
!         ENDIF        
!         
!         WRITE(14,"(2(I9))") nbseg, btype
!         DO nd = 1,nbseg
!           WRITE(14,"(I9)") mesh%fbnds(nd,bnd)
!         ENDDO   
!       ENDDO
!       
!       CLOSE(14)
      
      



      RETURN
      END SUBROUTINE rewrite_fort14

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE write_elem_nodes(eval,base)
      
      USE globals, ONLY: nverts

      IMPLICIT NONE
      
      INTEGER :: el,et,i,nnd
      TYPE(grid) :: eval,base

      
      PRINT "(A)", "Writing rimls element nodes..."  
      PRINT*, ""     

      OPEN(UNIT = 14, FILE = TRIM(out_direc) // "elem_nodes.d")      
      WRITE(14,"(2(I7,1x))") eval%ne,eval%hbp
      
      OPEN(UNIT = 15, FILE = TRIM(out_direc) // "xelem_nodes.d")      
      WRITE(15,"(2(I7,1x))") eval%ne,eval%hbp
      
      OPEN(UNIT = 16, FILE = TRIM(out_direc) // "yelem_nodes.d")      
      WRITE(16,"(2(I7,1x))") eval%ne,eval%hbp      
      
      DO el = 1,eval%ne
        et = eval%el_type(el)
        
        IF (mod(et,2) == 1) THEN   
          nnd = eval%nnds(5)
        ELSE IF (mod(et,2) == 0) THEN
          nnd = eval%nnds(6)          
        ENDIF    
     
        WRITE(14,"(2(I7),1x,60(e24.17,1x))") el,nnd,(eval%elhb(i,el), i = 1,nnd)
        WRITE(15,"(2(I7),1x,60(e24.17,1x))") el,nnd,(eval%elxyh(i,el,1), i = 1,nnd)
        WRITE(16,"(2(I7),1x,60(e24.17,1x))") el,nnd,(eval%elxyh(i,el,2), i = 1,nnd)
      ENDDO
      
      CLOSE(14)
      CLOSE(15)
      CLOSE(16)


      RETURN
      END SUBROUTINE write_elem_nodes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE write_hb_file(mesh)
      
      USE version, ONLY: version_information
      
      IMPLICIT NONE
      
      TYPE(grid) :: mesh
      INTEGER :: el,et,i,nnd
      INTEGER :: sind,eind
      CHARACTER(100) :: name
      CHARACTER(1) :: hbp_char      
      
      
      eind = INDEX(ADJUSTL(TRIM(mesh%grid_file)),".",.false.)
      name = ADJUSTL(TRIM(mesh%grid_file(1:eind-1)))           
      WRITE(hbp_char,"(I1)") mesh%hbp
      
      OPEN(UNIT = 13, FILE = ADJUSTL(TRIM(name)) // "_hbp" // hbp_char // "_rimls.hb") 
      
      CALL version_information(13)
      
      WRITE(13,"(A)") "-----------------------------------------------------------------------"           
      
      CALL write_input(13)           
      
      WRITE(13,"(A)") "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"         

      WRITE(13,"(2(I7,1x),A,I8)") mesh%ne,mesh%hbp      
      
      DO el = 1,mesh%ne
        et = mesh%el_type(el)
        
        IF (mod(et,2) == 1) THEN   
          nnd = mesh%nnds(5)
        ELSE IF (mod(et,2) == 0) THEN
          nnd = mesh%nnds(6)          
        ENDIF    

        WRITE(13,"(2(I7),1x,60(e24.17,1x))") el,nnd,(mesh%elhb(i,el), i = 1,nnd)        
      ENDDO      
      
      CLOSE(13)
      
      RETURN
      END SUBROUTINE write_hb_file


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      END MODULE write_results